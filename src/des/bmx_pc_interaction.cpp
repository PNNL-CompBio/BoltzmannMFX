//
//     Copyright (c) 2013 Battelle Memorial Institute
//     Licensed under modified BSD License. A copy of this license can be found
//     in the LICENSE file in the top level directory of this distribution.
//
#include <bmx_pc.H>
#include <bmx_dem_parms.H>
#include <bmx_bc_parms.H>
#include <bmx_cell_interaction_K.H>

using namespace amrex;

void BMXParticleContainer::EvolveParticles (Real dt,
                                            const Vector<MultiFab*> cost,
                                            std::string& knapsack_weight_type,
                                            int& nsubsteps)
{
    BL_PROFILE_REGION_START("bmx_dem::EvolveParticles()");
    BL_PROFILE("bmx_dem::EvolveParticles()");

    Real eps = std::numeric_limits<Real>::epsilon();

    for (int lev = 0; lev <= finest_level; lev++)
    {

    int n_at_lev = this->NumberOfParticlesAtLevel(lev);
    amrex::Print() << "In Evolve Particles with " << n_at_lev << " particles at level " << lev << std::endl;

    if (n_at_lev == 0) continue;

    amrex::Print() << "Evolving particles on level: " << lev
                   << " ... with fluid dt " << dt << std::endl;

    BMXCellInteraction *interaction = BMXCellInteraction::instance();
    std::vector<Real> fpar_vec = interaction->getForceParams();
    Real *fpar = &fpar_vec[0];
    /****************************************************************************
     * DEBUG flag toggles:                                                      *
     *   -> Print number of collisions                                          *
     *   -> Print max (over substeps) particle velocity at each time step       *
     *   -> Print max particle-wall and particle-particle forces                *
     ***************************************************************************/

    // Debug level controls the detail of debug output:
    //   -> debug_level = 0 : no debug output
    //   -> debug_level = 1 : debug output for every fluid step
    //   -> debug_level = 2 : debug output for every substep
    const int debug_level = 0;

    /****************************************************************************
     * Init substeps                                                            *
     ***************************************************************************/

    Real subdt;
    // TODO: This needs to be initialized somewhere else using a better rule
    DEM::dtsolid = dt;
    // des_init_time_loop(&dt, &nsubsteps, &subdt);
    if ( dt >= DEM::dtsolid )
    {
       nsubsteps = static_cast<int>(amrex::Math::ceil (  dt / DEM::dtsolid ));
       subdt     =  dt / nsubsteps;
    } else {
       nsubsteps = 1;
       subdt     = dt;
    }

    /****************************************************************************
     * Init temporary storage:                                                  *
     *   -> particle-particle, and particle-wall forces                         *
     *   -> particle-particle, and particle-wall torques                        *
     ***************************************************************************/
    // total force=particle+wall, particle force, wall force, torque?
    std::map<PairIndex, amrex::Gpu::DeviceVector<Real>> fc, pfor, wfor, torq;

    std::map<PairIndex, bool> tile_has_walls;

    for (BMXParIter pti(*this, lev); pti.isValid(); ++pti)
    {
        PairIndex index(pti.index(), pti.LocalTileIndex());
        fc[index]   = amrex::Gpu::DeviceVector<Real>();
        pfor[index] = amrex::Gpu::DeviceVector<Real>();
        wfor[index] = amrex::Gpu::DeviceVector<Real>();
        torq[index] = amrex::Gpu::DeviceVector<Real>();
    }

    /****************************************************************************
     * Iterate over sub-steps                                                   *
     ***************************************************************************/

    int n = 0; // Counts sub-steps

    while (n < nsubsteps)
    {
        // Redistribute particles every so often BUT always update the neighbour
        // list (Note that this fills the neighbour list after every
        // redistribute operation)
        if (n % 25 == 0) {
            clearNeighbors();
            Redistribute(0, 0, 0, 1);
            fillNeighbors();
            // send in "false" for sort_neighbor_list option

            buildNeighborList(BMXCheckPair(DEM::neighborhood,false), false);
        } else {
            updateNeighbors();
        }

        /********************************************************************
         * Particles routines                                               *
         *******************************************************************/
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for (BMXParIter pti(*this, lev); pti.isValid(); ++pti)
        {
            // Timer used for load-balancing
            Real wt = ParallelDescriptor::second();

            //const Box& bx = pti.tilebox(); // UNUSED_VARIABLE
            PairIndex index(pti.index(), pti.LocalTileIndex());

            auto& plev = GetParticles(lev);
            auto& ptile = plev[index];
            auto& aos   = ptile.GetArrayOfStructs();
            ParticleType* pstruct = aos().dataPtr();

            const int nrp = GetParticles(lev)[index].numRealParticles();

            // Number of particles including neighbor particles
            int ntot = nrp;

            // Particle-particle (and particle-wall) forces and torques. We need
            // these to be zero every time we start a new batch (i.e tile and
            // substep) of particles.
            fc[index].clear();
            fc[index].resize(6*ntot, 0.0);

            Real* fc_ptr = fc[index].dataPtr();

            // For debugging: keep track of particle-particle (pfor),
            // particle-wall (wfor) forces and torques (torq)
            pfor[index].clear();
            wfor[index].clear();
            torq[index].clear();
            pfor[index].resize(3*ntot, 0.0);
            wfor[index].resize(3*ntot, 0.0);
            torq[index].resize(3*ntot, 0.0);

            /********************************************************************
             * Particle-Particle collision forces (and torques)                 *
             *******************************************************************/

            BL_PROFILE_VAR("calc_particle_collisions()", calc_particle_collisions);

            auto nbor_data = m_neighbor_list[lev][index].data();

            constexpr Real small_number = 1.0e-15;

#if 0
            Real l_bndry_width   = BMXCellInteraction::p_bndry_width;
            Real l_stiffness     = BMXCellInteraction::p_stiffness;
            Real l_z_bndry_width = BMXCellInteraction::p_z_bndry_width;
            Real l_z_stiffness   = BMXCellInteraction::p_z_stiffness;
            Real l_z_wall        = BMXCellInteraction::p_z_wall;
            Real l_z_gravity     = BMXCellInteraction::p_z_gravity;
#endif

            // now we loop over the neighbor list and compute the forces
            int me = ParallelDescriptor::MyProc();
            amrex::ParallelFor(nrp,
                [nrp,pstruct,fc_ptr,nbor_data,subdt,ntot,fpar,me,n]
              AMREX_GPU_DEVICE (int i) noexcept
              {
                  auto& particle = pstruct[i];

                  RealVect pos1(particle.pos());
                  // clean up flags
                  particle.idata(intIdx::fuse_flag) = 0;
                  particle.idata(intIdx::split_flag) = 0;
                  particle.idata(intIdx::new_flag) = 0;
//                  printf("p[%d] ID: %d CPU: %d RX: %e RY: %e RZ: %e\n",me,
//                      static_cast<int>(particle.id()),static_cast<int>(particle.cpu()),
//                      pos1[0],pos1[1],pos1[2]);

                  const auto neighbs = nbor_data.getNeighbors(i);
                  for (auto mit = neighbs.begin(); mit != neighbs.end(); ++mit)
                  {
                      auto p2 = *mit;
                      const int j = mit.index();

                      Real dist_x = pos1[0] - p2.pos(0);
                      Real dist_y = pos1[1] - p2.pos(1);
                      Real dist_z = pos1[2] - p2.pos(2);

                      Real r2 = dist_x*dist_x +
                                dist_y*dist_y +
                                dist_z*dist_z;

//                      printf("SEPARATION: %f\n",sqrt(r2));
                      RealVect diff(dist_x,dist_y,dist_z);

                      Real r_lm = maxInteractionDistance(&particle.rdata(0),&p2.rdata(0),
                                                         &particle.idata(0),&p2.idata(0),
                                                         fpar[0]);
//                      printf("R_LM: %f\n",r_lm);

                      AMREX_ASSERT_WITH_MESSAGE(
                          not (particle.id() == p2.id() and
                               particle.cpu() == p2.cpu()),
                        "A particle should not be its own neighbor!");

                      if ( r2 <= (r_lm - small_number)*(r_lm - small_number) )
                      {
                       //   Real dist_mag = sqrt(r2);
                       //   Real dist_mag_inv = 1.e0/dist_mag;
                       //   RealVect normal(0.);
                       //   normal[0] = dist_x * dist_mag_inv;
                       //   normal[1] = dist_y * dist_mag_inv;
                       //   normal[2] = dist_z * dist_mag_inv;

                          RealVect v1(0.);
                          RealVect v2(0.);
                          RealVect rot1(0.);
                          RealVect rot2(0.);

                          if (p2.idata(intIdx::fuse_flag) != 0) {
                            printf("p[%d] Found FUSING particle j in force"
                                " loop id: %d cpu: %d\n",
                                me,p2.idata(intIdx::id),p2.idata(intIdx::cpu));
                          }
                          evaluateForce(&diff[0],&particle.rdata(0),
                              &p2.rdata(0), &particle.idata(0),
                              &p2.idata(0), &v1[0], &v2[0], &rot1[0],
                              &rot2[0], fpar, me, i, j, n);

//                          printf("w1x: %e w1y: %e w1z: %e w2x: %e w2y: %e w2z: %e\n",
//                              rot1[0],rot1[1],rot1[2],rot2[0],rot2[1],rot2[2]);
//                          printf("v1x: %e r1x: %e v2x: %e r2x: %e\n",
//                              v1[0],pos1[0],v2[0],p2.pos(0));
#ifdef _OPENMP
#pragma omp critical
                          {
#endif
                            amrex::Gpu::Atomic::Add(&fc_ptr[i         ], v1[0]);
                            amrex::Gpu::Atomic::Add(&fc_ptr[i + ntot  ], v1[1]);
                            amrex::Gpu::Atomic::Add(&fc_ptr[i + 2*ntot], v1[2]);
                            amrex::Gpu::Atomic::Add(&fc_ptr[i + 3*ntot], rot1[0]);
                            amrex::Gpu::Atomic::Add(&fc_ptr[i + 4*ntot], rot1[1]);
                            amrex::Gpu::Atomic::Add(&fc_ptr[i + 5*ntot], rot1[2]);

#if 0
                            if (j < nrp)
                            {
                              amrex::Gpu::Atomic::Add(&fc_ptr[j         ], v2[0]);
                              amrex::Gpu::Atomic::Add(&fc_ptr[j + ntot  ], v2[1]);
                              amrex::Gpu::Atomic::Add(&fc_ptr[j + 2*ntot], v2[2]);
                              amrex::Gpu::Atomic::Add(&fc_ptr[j + 3*ntot], rot2[0]);
                              amrex::Gpu::Atomic::Add(&fc_ptr[j + 4*ntot], rot2[1]);
                              amrex::Gpu::Atomic::Add(&fc_ptr[j + 5*ntot], rot2[2]);
                            }
#endif
#ifdef _OPENMP
                          }
#endif
                          // TODO: Do we need an OPENMP pragma here?

                      }
                  } // end of neighbor loop
                  RealVect vcom(0.);
                  RealVect vrot(0.);
                  evaluateSurfaceForce(&pos1[0],&particle.rdata(0),
                      &particle.idata(0),&vcom[0],&vrot[0],fpar);
                  amrex::Gpu::Atomic::Add(&fc_ptr[i         ], vcom[0]);
                  amrex::Gpu::Atomic::Add(&fc_ptr[i + ntot  ], vcom[1]);
                  amrex::Gpu::Atomic::Add(&fc_ptr[i + 2*ntot], vcom[2]);
                  amrex::Gpu::Atomic::Add(&fc_ptr[i + 3*ntot], vrot[0]);
                  amrex::Gpu::Atomic::Add(&fc_ptr[i + 4*ntot], vrot[1]);
                  amrex::Gpu::Atomic::Add(&fc_ptr[i + 5*ntot], vrot[2]);
              }); // end of loop over particles

            amrex::Gpu::Device::synchronize();

            // Debugging: copy data from the fc (all forces) vector to the wfor
            // (wall forces) vector. Note that since fc already contains the
            // wall forces, these need to be subtracted here.
            if (debug_level > 0)
            {
                for (size_t i = 0; i < pfor[index].size(); i++ ) {
                    pfor[index][i] = fc[index][i] - wfor[index][i];
                }
            }

            BL_PROFILE_VAR_STOP(calc_particle_collisions);

            BL_PROFILE_VAR("des::update_particle_velocity_and_position()", des_time_march);
            /********************************************************************
             * Move particles based on collision forces and torques             *
             *******************************************************************/

            const auto p_lo = Geom(lev).ProbLoArray();
            const auto p_hi = Geom(lev).ProbHiArray();

            int x_lo_bc = BC::domain_bc[0];
            int x_hi_bc = BC::domain_bc[1];
            int y_lo_bc = BC::domain_bc[2];
            int y_hi_bc = BC::domain_bc[3];
            int z_lo_bc = BC::domain_bc[4];
            int z_hi_bc = BC::domain_bc[5];

            bool verbose = p_verbose;
            amrex::ParallelFor(nrp,
              [pstruct,subdt,fc_ptr,ntot,eps,p_hi,p_lo,
               x_lo_bc,x_hi_bc,y_lo_bc,y_hi_bc,z_lo_bc,z_hi_bc,
              verbose]
              AMREX_GPU_DEVICE (int i) noexcept
              {
                auto& particle = pstruct[i];

                RealVect ppos(particle.pos());

                particle.rdata(realIdx::velx) = fc_ptr[i];
                particle.rdata(realIdx::vely) = fc_ptr[i+ntot];
                particle.rdata(realIdx::velz) = fc_ptr[i+2*ntot];

                particle.rdata(realIdx::wx) = fc_ptr[i+3*ntot];
                particle.rdata(realIdx::wy) = fc_ptr[i+4*ntot];
                particle.rdata(realIdx::wz) = fc_ptr[i+5*ntot];

                ppos[0] += subdt * particle.rdata(realIdx::velx);
                ppos[1] += subdt * particle.rdata(realIdx::vely);
                ppos[2] += subdt * particle.rdata(realIdx::velz);

                // Modify orientation based on angular momentum
                if (particle.idata(intIdx::cell_type) == cellType::FUNGI) {
                  Real theta = particle.rdata(realIdx::theta);
                  Real phi = particle.rdata(realIdx::phi);
                  Real cp = cos(phi);
                  Real sp = sin(phi);
                  Real ct = cos(theta);
                  Real st = sin(theta);
                  // Real x = sp*st;
                  // Real y = cp*st;
                  // Real z = ct;
                  // construct matrix to rotate x-axis to segment orientation
                  Real rot[3][3];
                  rot[0][0] = cp*st;
                  rot[0][1] = -sp;
                  rot[0][2] = -cp*ct;
                  rot[1][0] = sp*st;
                  rot[1][1] = cp;
                  rot[1][2] = -sp*ct;
                  rot[2][0] = ct;
                  rot[2][1] = 0.0;
                  rot[2][2] = st;
                  // get angular momentum
                  Real lm[3];
                  lm[0] = particle.rdata(realIdx::wx);
                  lm[1] = particle.rdata(realIdx::wy);
                  lm[2] = particle.rdata(realIdx::wz);
                  // transform angular momentum using inverse rotation
                  // and calculate rotational velocity
                  Real om[3];
                  for (int ii=0; ii<3; ii++) {
                    om[ii] = 0.0;
                    for (int jj=0; jj<3; jj++) {
                      om[ii] += rot[jj][ii]*lm[jj];
                    }
                  }
                  Real clen = particle.rdata(realIdx::c_length);
                  om[0] = 0.0;
                  om[1] *= 2.0/(clen*clen);
                  om[2] *= 2.0/(clen*clen);
                  // Find rotation angle of angular momentum with respect to
                  // z-axis
                  Real on = sqrt(om[1]*om[1]+om[2]*om[2]);
                  if (on > 0.0) {
                    om[1] /= on;
                    om[2] /= on;
                  }
                  // Construct matrix to rotate system about x-axis so that rotational
                  // velocity is along z-axis. Cosine of the angle theta with
                  // the z-axis is just om[2] and sine of theta is om[1]
                  ct = om[2];
                  st = om[1];
                  Real orot[3][3];
                  orot[0][0] = 1.0;
                  orot[0][1] = 0.0;
                  orot[0][2] = 0.0;
                  orot[1][0] = 0.0;
                  orot[1][1] = ct;
                  orot[1][2] = -st;
                  orot[2][0] = 0.0;
                  orot[2][1] = st;
                  orot[2][2] = ct;
                  // Cylinder is currently aligned along x-axis, so this
                  // rotation has no effect on it. Calculate how much system
                  // rotates about the z-axis in 1 time step and then calculate
                  // how much x-axis is rotated
                  Real dtheta = on*subdt;
                  Real ndir[3]; 
                  ndir[0] = cos(dtheta);
                  ndir[1] = sin(dtheta);
                  ndir[2] = 0.0;
                  // apply inverse of orot to ndir (inverse is equal to
                  // transpose)
                  for (int ii=0; ii<3; ii++) {
                    om[ii] = 0.0;
                    for (int jj=0; jj<3; jj++) {
                      om[ii] += orot[jj][ii]*ndir[jj];
                    }
                  }

                  // apply rot to new direction to recover final orientation
                  for (int ii=0; ii<3; ii++) {
                    ndir[ii] = 0.0;
                    for (int jj=0; jj<3; jj++) {
                      ndir[ii] += rot[ii][jj]*om[jj];
                    }
                  }
                  // get orientation angles
                  theta = acos(ndir[2]);
                  on = sqrt(ndir[0]*ndir[0]+ndir[1]*ndir[1]);
                  if (on > 0.0) {
                    ndir[0] /= on;
                  }
                  phi = acos(ndir[0]);
                  if (ndir[1] < 0.0) {
                    phi = 2.0*M_PI-phi;
                  }
                  particle.rdata(realIdx::theta) = theta;
                  particle.rdata(realIdx::phi) = phi;
                }

                particle.pos(0) = ppos[0];
                particle.pos(1) = ppos[1];
                particle.pos(2) = ppos[2];

#if !defined(AMREX_USE_GPU)
                if (verbose) {
                  char sbuf[128];
                  int p_id = particle.id();
                  sprintf(sbuf,"particle: %d position: %14.8f %14.8f %14.8f"
                      " volume: %16.8e length: %14.8f",p_id,particle.pos(0),
                      particle.pos(1),particle.pos(2), particle.rdata(realIdx::vol),
                      particle.rdata(realIdx::c_length));
                  std::cout << sbuf << std::endl;
                }
#endif
              });

            BL_PROFILE_VAR_STOP(des_time_march);

            amrex::Gpu::synchronize();


            /********************************************************************
             * Update runtime cost (used in load-balancing)                     *
             *******************************************************************/

            if (cost[lev])
            {
                // Runtime cost is either (weighted by tile box size):
                //   * time spent
                //   * number of particles
                const Box& tbx = pti.tilebox();
                if (knapsack_weight_type == "RunTimeCosts")
                {
                    wt = (ParallelDescriptor::second() - wt) / tbx.d_numPts();
                }
                else if (knapsack_weight_type == "NumParticles")
                {
                    wt = nrp / tbx.d_numPts();
                }
                (*cost[lev])[pti].plus<RunOn::Device>(wt, tbx);
            }
        }

        // Update substep count
        n += 1;

    } // end of loop over substeps

    // Redistribute particles at the end of all substeps (note that the particle
    // neighbour list needs to be reset when redistributing).
    clearNeighbors();
    Redistribute(0, 0, 0, 1);

    } // lev

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif

    BL_PROFILE_REGION_STOP("bmx_dem::EvolveParticles()");
}

/*******************************************************************************
 *    *
 ******************************************************************************/
void BMXParticleContainer::InitBonds (const Vector<MultiFab*> cost,
                                              std::string& knapsack_weight_type)
{
  printf("Calling InitBonds\n");
    BL_PROFILE_REGION_START("bmx_dem::InitBonds()");
    BL_PROFILE("bmx_dem::InitBonds()");

    Real eps = std::numeric_limits<Real>::epsilon();

    Real rlim = 1.0e-7;

    BMXCellInteraction *interaction = BMXCellInteraction::instance();
    std::vector<Real> fpar_vec = interaction->getForceParams();
    Real *fpar = &fpar_vec[0];

    for (int lev = 0; lev <= finest_level; lev++)
    {

    int n_at_lev = this->NumberOfParticlesAtLevel(lev);

    if (n_at_lev == 0) continue;


    // Debug level controls the detail of debug output:
    //   -> debug_level = 0 : no debug output
    //   -> debug_level = 1 : debug output for every fluid step
    //   -> debug_level = 2 : debug output for every substep
    const int debug_level = 0;

    // update the neighbour list with fresh data
    clearNeighbors();
    Redistribute(0, 0, 0, 1);
    fillNeighbors();
    // send in "false" for sort_neighbor_list option

    buildNeighborList(BMXCheckPair(DEM::neighborhood, false), false);

    /********************************************************************
     * Particles routines                                               *
     *******************************************************************/
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (BMXParIter pti(*this, lev); pti.isValid(); ++pti)
    {
      // Timer used for load-balancing
      Real wt = ParallelDescriptor::second();

      //const Box& bx = pti.tilebox(); // UNUSED_VARIABLE
      PairIndex index(pti.index(), pti.LocalTileIndex());

      auto& plev = GetParticles(lev);
      auto& ptile = plev[index];
      auto& aos   = ptile.GetArrayOfStructs();
      ParticleType* pstruct = aos().dataPtr();

      const int nrp = GetParticles(lev)[index].numRealParticles();

      // Number of particles including neighbor particles
      int ntot = nrp;

      /********************************************************************
       * Particle-Particle collision forces (and torques)                 *
       *******************************************************************/

      BL_PROFILE_VAR("set_bonds()", set_bonds);

      auto nbor_data = m_neighbor_list[lev][index].data();

      constexpr Real small_number = 1.0e-15;

      // now we loop over the neighbor list and compute the forces
      int me = ParallelDescriptor::MyProc();
      amrex::ParallelForRNG(nrp,
          [nrp,pstruct,nbor_data,rlim,fpar,me]
          AMREX_GPU_DEVICE (int i, amrex::RandomEngine const& engine) noexcept
//            [=] AMREX_GPU_DEVICE (int i, amrex::RandomEngine const& engine) noexcept
          {
          auto& particle = pstruct[i];

          /*
          if (particle.idata(intIdx::position) == siteLocation::TIP &&
              particle.idata(intIdx::n_bnds) > 2) {
          int *idata = &particle.idata(0);
          printf("particle id: %d cpu: %d is Tip with %d bonds\n",
              idata[intIdx::id],idata[intIdx::cpu],idata[intIdx::n_bnds]);
          }
          */
          RealVect pos1(particle.pos());

          const auto neighbs = nbor_data.getNeighbors(i);
          for (auto mit = neighbs.begin(); mit != neighbs.end(); ++mit)
          {
          auto p2 = *mit;
          const int j = mit.index();

          Real dist_x = pos1[0] - p2.pos(0);
          Real dist_y = pos1[1] - p2.pos(1);
          Real dist_z = pos1[2] - p2.pos(2);

          Real r2 = dist_x*dist_x +
            dist_y*dist_y +
            dist_z*dist_z;

          printf("R: %e\n",sqrt(r2));
          RealVect diff(dist_x,dist_y,dist_z);


          AMREX_ASSERT_WITH_MESSAGE(
              not (particle.id() == p2.id() and
                particle.cpu() == p2.cpu()),
              "A particle should not be its own neighbor!");

          Real r_lm = maxInteractionDistance(&particle.rdata(0),&p2.rdata(0),
              &particle.idata(0),&p2.idata(0), fpar[0]);

          if ( r2 <= (r_lm - small_number)*(r_lm - small_number) )
          {

            // create bond if fungi segment ends are close enough to
            // each other
            fixBonds(&diff[0], &particle.rdata(0), &particle.idata(0),
                &p2.rdata(0), &p2.idata(0), rlim);

            // TODO: Do we need an OPENMP pragma here?

          }
          } // end of neighbor loop
          }); // end of loop over particles

      amrex::Gpu::Device::synchronize();

      BL_PROFILE_VAR_STOP(set_bonds);

      /********************************************************************
       * Update runtime cost (used in load-balancing)                     *
       *******************************************************************/

      if (cost[lev])
      {
        // Runtime cost is either (weighted by tile box size):
        //   * time spent
        //   * number of particles
        const Box& tbx = pti.tilebox();
        if (knapsack_weight_type == "RunTimeCosts")
        {
          wt = (ParallelDescriptor::second() - wt) / tbx.d_numPts();
        }
        else if (knapsack_weight_type == "NumParticles")
        {
          wt = nrp / tbx.d_numPts();
        }
        (*cost[lev])[pti].plus<RunOn::Device>(wt, tbx);
      }
    }

    // Redistribute particles at the end of all substeps (note that the particle
    // neighbour list needs to be reset when redistributing).
    //clearNeighbors();
    //Redistribute(0, 0, 0, 1);
    //updateNeighbors();

    } // lev

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif

    BL_PROFILE_REGION_STOP("bmx_dem::InitBonds()");
}
