#include <bmx_des_K.H>

#include <bmx_dem_parms.H>
#include <bmx_bc_parms.H>

using namespace amrex;

int  BMXParticleContainer::domain_bc[6] {0};

BMXParticleContainer::BMXParticleContainer (AmrCore* amr_core)
    : NeighborParticleContainer<realData::count,intData::count>
      (amr_core->GetParGDB(), 1)
{
    ReadStaticParameters();

    this->SetVerbose(0);

    nlev = amr_core->maxLevel() + 1;
}

void BMXParticleContainer::AllocData ()
{
    reserveData();
    resizeData();
}

void BMXParticleContainer::PrintParticleCounts ()
{
  const int lev = 0;
  amrex::AllPrintToFile("load_balance") << "Particles on each box: \n";
  long local_count = 0;
  for (BMXParIter pti(*this, lev); pti.isValid(); ++pti)
    {
      long np = pti.numParticles();
      local_count += np;
      amrex::AllPrintToFile("load_balance") << "Box:" << pti.index() << ", count: " << np << std::endl;
    }
  amrex::AllPrintToFile("load_balance") << "Total for this process: " << local_count << std::endl << std::endl;
}

void BMXParticleContainer::printParticles ()
{
    const int lev = 0;
    auto& plevel = GetParticles(lev);

    for (auto& kv : plevel)
    {
       const auto& particles = kv.second.GetArrayOfStructs();
       auto& soa = kv.second.GetStructOfArrays();
       auto p_realarray = soa.realarray();
       auto p_intarray = soa.intarray();

       for (int i = 0; i < particles.numParticles(); ++i)
       {
          std::cout << "Particle ID  = " << i << " " << std::endl;
          std::cout << "X            = " << particles[i].pos(0) << " " << std::endl;
          std::cout << "Y            = " << particles[i].pos(1) << " " << std::endl;
          std::cout << "Z            = " << particles[i].pos(2) << " " << std::endl;
          std::cout << "state        = " << p_intarray[intData::state][i] << " " << std::endl;
          std::cout << "phase        = " << p_intarray[intData::phase][i] << " " << std::endl;
          std::cout << "Real properties = " << std::endl;

          for (int j = 0; j < realData::count; j++)
            std::cout << "property " << j << "  = " << p_realarray[j][i] << " " << std::endl;

          std::cout << std::endl;
       }
    }
}

void BMXParticleContainer::ReadStaticParameters ()
{
    static bool initialized = false;

    if (!initialized)
        initialized = true;
}

void BMXParticleContainer::EvolveParticles (int lev,
                                             int nstep,
                                             Real dt,
                                             Real time,
                                             MultiFab* cost,
                                             std::string& knapsack_weight_type,
                                             int& nsubsteps)
{
    BL_PROFILE_REGION_START("bmx_dem::EvolveParticles()");
    BL_PROFILE("bmx_dem::EvolveParticles()");

    Real eps = std::numeric_limits<Real>::epsilon();

    amrex::Print() << "Evolving particles on level: " << lev
                   << " ... with fluid dt " << dt << std::endl;

    /****************************************************************************
     * Geometry                                                                 *
     ***************************************************************************/

    const Real* dx = Geom(lev).CellSize();

    /****************************************************************************
     * Init substeps                                                            *
     ***************************************************************************/

    Real subdt;
    // des_init_time_loop(&dt, &nsubsteps, &subdt);
    // if ( dt >= DEM::dtsolid )
    // {
    //    nsubsteps = amrex::Math::ceil (  dt / DEM::dtsolid );
    //    subdt     =  dt / nsubsteps;
    // } else 
    {
       nsubsteps = 1;
       subdt     = dt;
    }

    /****************************************************************************
     * Init temporary storage:                                                  *
     *   -> particle-particle, and particle-wall forces                         *
     *   -> particle-particle, and particle-wall torques                        *
     ***************************************************************************/
    std::map<PairIndex, Gpu::DeviceVector<Real>> tow;
    std::map<PairIndex, Gpu::DeviceVector<Real>> fc, pfor, wfor;

    std::map<PairIndex, bool> tile_has_walls;

    for (BMXParIter pti(*this, lev); pti.isValid(); ++pti)
    {
        const Box& bx = pti.tilebox();
        PairIndex index(pti.index(), pti.LocalTileIndex());
        tow[index]  = Gpu::DeviceVector<Real>();
        fc[index]   = Gpu::DeviceVector<Real>();
        pfor[index] = Gpu::DeviceVector<Real>();
        wfor[index] = Gpu::DeviceVector<Real>();
    }

    /****************************************************************************
     * Iterate over sub-steps                                                   *
     ***************************************************************************/

    int n = 0; // Counts sub-steps

    while (n < nsubsteps)
    {
        // Redistribute particles ever so often BUT always update the neighbour
        // list (Note that this fills the neighbour list after every
        // redistribute operation)
        if (n % 25 == 0) {
            clearNeighbors();
            Redistribute(0, 0, 0, 1);
            fillNeighbors();
            // send in "false" for sort_neighbor_list option

            buildNeighborList(BMXCheckPair(DEM::neighborhood), false);
        } else {
            updateNeighbors();
        }

#if 0
        /********************************************************************
         * Particles routines                                               *
         *******************************************************************/
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
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
            amrex::Print() << "THERE ARE " << nrp << " PARTICLES IN PLAY " << std::endl;

            auto& soa   = pti.GetStructOfArrays();
            auto p_realarray = soa.realarray();
            auto p_intarray = soa.intarray();

            // Number of particles including neighbor particles
            int ntot = nrp;

            // Particle-particle (and particle-wall) forces and torques. We need
            // these to be zero every time we start a new batch (i.e tile and
            // substep) of particles.
            tow[index].clear();
            fc[index].clear();
            tow[index].resize(3*ntot, 0.0);
            fc[index].resize(3*ntot, 0.0);

            Real* fc_ptr = fc[index].dataPtr();
            Real* tow_ptr = tow[index].dataPtr();

            // For debugging: keep track of particle-particle (pfor) and
            // particle-wall (wfor) forces
            pfor[index].clear();
            wfor[index].clear();
            pfor[index].resize(3*ntot, 0.0);
            wfor[index].resize(3*ntot, 0.0);

            /********************************************************************
             * Particle-Particle collision forces (and torques)                 *
             *******************************************************************/

            BL_PROFILE_VAR("calc_particle_collisions()", calc_particle_collisions);

            auto nbor_data = m_neighbor_list[lev][index].data();

            constexpr Real small_number = 1.0e-15;

            // now we loop over the neighbor list and compute the forces
            amrex::ParallelFor(nrp,
                [nrp,pstruct,p_realarray,p_intarray,fc_ptr,tow_ptr,nbor_data,
#if defined(AMREX_DEBUG) || defined(AMREX_USE_ASSERTION)
                 eps,
#endif
                 subdt,ntot]
              AMREX_GPU_DEVICE (int i) noexcept
              {
                  auto particle = pstruct[i];

                  RealVect pos1(particle.pos());

                  const auto neighbs = nbor_data.getNeighbors(i);
                  for (auto mit = neighbs.begin(); mit != neighbs.end(); ++mit)
                  {
                      const auto p2 = *mit;
                      const int j = mit.index();

                      Real dist_x = p2.pos(0) - pos1[0];
                      Real dist_y = p2.pos(1) - pos1[1];
                      Real dist_z = p2.pos(2) - pos1[2];

                      Real r2 = dist_x*dist_x +
                                dist_y*dist_y +
                                dist_z*dist_z;

                      const Real p1radius = p_realarray[realData::radius][i];
                      const Real p2radius = p_realarray[realData::radius][j];

                      Real r_lm = p1radius + p2radius;

                      AMREX_ASSERT_WITH_MESSAGE(
                          not (particle.id() == p2.id() and
                               particle.cpu() == p2.cpu()),
                        "A particle should not be its own neighbor!");

                      if ( r2 <= (r_lm - small_number)*(r_lm - small_number) )
                      {
#if 0
                          Real dist_mag = sqrt(r2);

                          AMREX_ASSERT(dist_mag >= eps);

                          Real dist_mag_inv = 1.e0/dist_mag;

                          RealVect normal(0.);
                          normal[0] = dist_x * dist_mag_inv;
                          normal[1] = dist_y * dist_mag_inv;
                          normal[2] = dist_z * dist_mag_inv;

                          Real overlap_n = r_lm - dist_mag;
                          Real vrel_trans_norm;
                          RealVect vrel_t(0.);

                          RealVect p1vel(p_realarray[realData::velx][i],
                                         p_realarray[realData::vely][i],
                                         p_realarray[realData::velz][i]);

                          RealVect p2vel(p_realarray[realData::velx][j],
                                         p_realarray[realData::vely][j],
                                         p_realarray[realData::velz][j]);

                          cfrelvel(p1vel, p2vel, p1radius, p2radius, 
                                   vrel_trans_norm, vrel_t, normal, dist_mag);

                          int phase1 = p_intarray[intData::phase][i];
                          int phase2 = p_intarray[intData::phase][j];

                          Real kn_des = local_kn;
                          Real etan_des = local_etan(phase1-1,phase2-1);

                          // NOTE - we don't use the tangential components right now,
                          // but we might in the future
                          // Real kt_des = DEM::kt;
                          // Real etat_des = DEM::etat[phase1-1][phase2-1];

                          RealVect fn(0.);
                          RealVect ft(0.);
                          RealVect overlap_t(0.);
                          Real mag_overlap_t(0.);

                          // calculate the normal contact force
                          fn[0] = -(kn_des*overlap_n*normal[0]
                                  + etan_des*vrel_trans_norm*normal[0]);
                          fn[1] = -(kn_des*overlap_n*normal[1]
                                  + etan_des*vrel_trans_norm*normal[1]);
                          fn[2] = -(kn_des*overlap_n*normal[2]
                                  + etan_des*vrel_trans_norm*normal[2]);

#ifdef _OPENMP
#pragma omp critical
                          {
#endif
                            Gpu::Atomic::Add(&fc_ptr[i         ], fn[0] + ft[0]);
                            Gpu::Atomic::Add(&fc_ptr[i + ntot  ], fn[1] + ft[1]);
                            Gpu::Atomic::Add(&fc_ptr[i + 2*ntot], fn[2] + ft[2]);

                            if (j < nrp)
                            {
                              Gpu::Atomic::Add(&fc_ptr[j         ], -(fn[0] + ft[0]));
                              Gpu::Atomic::Add(&fc_ptr[j + ntot  ], -(fn[1] + ft[1]));
                              Gpu::Atomic::Add(&fc_ptr[j + 2*ntot], -(fn[2] + ft[2]));
                            }
#ifdef _OPENMP
                          }
#endif

#endif // end of "if 0"
                      }
                  }
              });

            Gpu::Device::synchronize();

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

            amrex::ParallelFor(nrp,
              [pstruct,subdt,fc_ptr,ntot,tow_ptr,eps,p_hi,p_lo,
               x_lo_bc,x_hi_bc,y_lo_bc,y_hi_bc,z_lo_bc,z_hi_bc]
              AMREX_GPU_DEVICE (int i) noexcept
              {
                auto& p= pstruct[i];

                p.rdata(realData::velx) += subdt * fc_ptr[i];
                p.rdata(realData::vely) += subdt * fc_ptr[i+ntot];
                p.rdata(realData::velz) += subdt * fc_ptr[i+2*ntot];

                p.pos(0) += subdt * p.rdata(realData::velx);
                p.pos(1) += subdt * p.rdata(realData::vely);
                p.pos(2) += subdt * p.rdata(realData::velz);

                if (x_lo_bc and p.pos(0) < p_lo[0])
                {
                  p.pos(0) = p_lo[0] + eps;
                  p.rdata(realData::velx) = -p.rdata(realData::velx);
                }
                else if (x_hi_bc and p.pos(0) > p_hi[0])
                {
                  p.pos(0) = p_hi[0] - eps;
                  p.rdata(realData::velx) = -p.rdata(realData::velx);
                }
                else if (y_lo_bc and p.pos(1) < p_lo[1])
                {
                  p.pos(1) = p_lo[1] + eps;
                  p.rdata(realData::vely) = -p.rdata(realData::vely);
                }
                else if (y_hi_bc and p.pos(1) > p_hi[1])
                {
                  p.pos(1) = p_hi[1] - eps;
                  p.rdata(realData::vely) = -p.rdata(realData::vely);
                }
                else if (z_lo_bc and p.pos(2) < p_lo[2])
                {
                  p.pos(2) = p_lo[2] + eps;
                  p.rdata(realData::velz) = -p.rdata(realData::velz);
                }
                else if (z_hi_bc and p.pos(2) > p_hi[2])
                {
                  p.pos(2) = p_hi[2] - eps;
                  p.rdata(realData::velz) = -p.rdata(realData::velz);
                }
              });

            BL_PROFILE_VAR_STOP(des_time_march);

            Gpu::synchronize();

            /********************************************************************
             * Update runtime cost (used in load-balancing)                     *
             *******************************************************************/

            if (cost)
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
                (*cost)[pti].plus<RunOn::Device>(wt, tbx);
            }
        }
#endif

        // Update substep count
        n += 1;

    } // end of loop over substeps

    // Redistribute particles at the end of all substeps (note that the particle
    // neighbour list needs to be reset when redistributing).
    clearNeighbors();
    Redistribute(0, 0, 0, 1);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (BMXParIter pti(*this, lev); pti.isValid(); ++pti)
    {
        const int nrp   = NumberOfParticles(pti);
        void* particles = pti.GetArrayOfStructs().data();
    }

    amrex::Print() << "done. \n";

    BL_PROFILE_REGION_STOP("bmx_dem::EvolveParticles()");
}

void BMXParticleContainer::writeAllAtLevel (int lev)
{
    // Not threaded because its print to terminal
    for (BMXParIter pti(*this, lev); pti.isValid(); ++pti)
    {
        auto& particles = pti.GetArrayOfStructs();
        int np = pti.numParticles();
        Gpu::HostVector<ParticleType> host_particles(np);
        Gpu::copy(Gpu::deviceToHost, particles.begin(), particles.end(), host_particles.begin());

        for (const auto& p: host_particles)
        {
           const IntVect& iv = Index(p, lev);

           RealVect xyz(p.pos(0), p.pos(1), p.pos(2));
           std::cout << " id " << p.id()
                << " index " << iv
                << " position " << xyz << std::endl;
       }
    }
}

