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
     * Geometry                                                                 *
     ***************************************************************************/

    const Real* dx = Geom(lev).CellSize();

    /****************************************************************************
     * Init substeps                                                            *
     ***************************************************************************/

    Real subdt;
    // des_init_time_loop(&dt, &nsubsteps, &subdt);
    if ( dt >= DEM::dtsolid )
    {
       nsubsteps = amrex::Math::ceil (  dt / DEM::dtsolid );
       subdt     =  dt / nsubsteps;
    } else {
       nsubsteps = 1;
       subdt     = dt;
    }

    /****************************************************************************
     * Get particle EB geometric info
     ***************************************************************************/
    const FabArray<EBCellFlagFab>* flags = &(ebfactory->getMultiEBCellFlagFab());

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

    int ncoll_total = 0;  // Counts total number of collisions
    loc_maxvel  = RealVect(0., 0., 0.);  // Tracks max (absolute) velocity
    loc_maxpfor = RealVect(0., 0., 0.);  // Tracks max particle-particle force
    loc_maxwfor = RealVect(0., 0., 0.);  // Tracks max particle-wall force
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

        /********************************************************************
         * Compute number of Particle-Particle collisions
         *******************************************************************/
        int ncoll = 0;  // Counts number of collisions (over sub-steps)

        if (debug_level > 0) 
        {
#ifdef AMREX_USE_GPU
          if (Gpu::inLaunchRegion())
          {
            // Reduce sum operation for ncoll
            ReduceOps<ReduceOpSum> reduce_op;
            ReduceData<int> reduce_data(reduce_op);
            using ReduceTuple = typename decltype(reduce_data)::Type;

            for (BMXParIter pti(*this, lev); pti.isValid(); ++pti)
            {
              PairIndex index(pti.index(), pti.LocalTileIndex());

              const int nrp = GetParticles(lev)[index].numRealParticles();

              auto& plev = GetParticles(lev);
              auto& ptile = plev[index];
              auto& aos   = ptile.GetArrayOfStructs();
              ParticleType* pstruct = aos().dataPtr();

              auto& soa = ptile.GetStructOfArrays();
              auto p_realarray = soa.realarray();

              auto nbor_data = m_neighbor_list[lev][index].data();

              constexpr Real small_number = 1.0e-15;

              reduce_op.eval(nrp, reduce_data, [pstruct,p_realarray,
                  nbor_data,small_number]
                AMREX_GPU_DEVICE (int i) -> ReduceTuple
              {
                int l_ncoll(0);

                ParticleType p1 = pstruct[i];
                const RealVect pos1 = p1.pos();
                const Real radius1 = p_realarray[SoArealData::radius][i];

                const auto neighbs = nbor_data.getNeighbors(i);
                for (auto mit = neighbs.begin(); mit != neighbs.end(); ++mit)
                {
                  const auto p2 = *mit;
                  const int j = mit.index();

                  const RealVect pos2 = p2.pos();
                  const Real radius2 = p_realarray[SoArealData::radius][j];

                  Real r2 = (pos1 - pos2).radSquared();

                  Real r_lm = radius1 + radius2;

                  if (r2 <= (r_lm-small_number)*(r_lm-small_number))
                    l_ncoll = 1;
                }

                return {l_ncoll};
              });
            }

            ReduceTuple host_tuple = reduce_data.value();
            ncoll += amrex::get<0>(host_tuple);
          }
          else
#endif
          {
#ifdef _OPENMP
#pragma omp parallel reduction(+:ncoll) if (Gpu::notInLaunchRegion())
#endif
            for (BMXParIter pti(*this, lev); pti.isValid(); ++pti)
            {
              PairIndex index(pti.index(), pti.LocalTileIndex());

              const int nrp = GetParticles(lev)[index].numRealParticles();

              auto& plev = GetParticles(lev);
              auto& ptile = plev[index];
              auto& aos   = ptile.GetArrayOfStructs();
              ParticleType* pstruct = aos().dataPtr();

              auto& soa = ptile.GetStructOfArrays();
              auto p_realarray = soa.realarray();

              auto nbor_data = m_neighbor_list[lev][index].data();

              constexpr Real small_number = 1.0e-15;

              for(int i(0); i < nrp; ++i)
              {
                ParticleType p1 = pstruct[i];
                const RealVect pos1 = p1.pos();
                const Real radius1 = p_realarray[SoArealData::radius][i];

                const auto neighbs = nbor_data.getNeighbors(i);
                for (auto mit = neighbs.begin(); mit != neighbs.end(); ++mit)
                {
                  const auto p2 = *mit;
                  const int j = mit.index();

                  const RealVect pos2 = p2.pos();
                  const Real radius2 = p_realarray[SoArealData::radius][j];

                  Real r2 = (pos1 - pos2).radSquared();

                  Real r_lm = radius1 + radius2;

                  if (r2 <= (r_lm-small_number)*(r_lm-small_number))
                  {
                    ncoll += 1;
                  }
                }
              }
            }
          } 
        } // end if (debug_level > 0)

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
                 subdt,ntot,local_mew=DEM::mew,local_kn=DEM::kn,
                 local_etan=DEM::etan]
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

                      const Real p1radius = p_realarray[SoArealData::radius][i];
                      const Real p2radius = p_realarray[SoArealData::radius][j];

                      Real r_lm = p1radius + p2radius;

                      AMREX_ASSERT_WITH_MESSAGE(
                          not (particle.id() == p2.id() and
                               particle.cpu() == p2.cpu()),
                        "A particle should not be its own neighbor!");

                      if ( r2 <= (r_lm - small_number)*(r_lm - small_number) )
                      {
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

                          RealVect p1vel(p_realarray[SoArealData::velx][i],
                                         p_realarray[SoArealData::vely][i],
                                         p_realarray[SoArealData::velz][i]);

                          RealVect p2vel(p_realarray[SoArealData::velx][j],
                                         p_realarray[SoArealData::vely][j],
                                         p_realarray[SoArealData::velz][j]);

                          RealVect p1omega(p_realarray[SoArealData::omegax][i],
                                           p_realarray[SoArealData::omegay][i],
                                           p_realarray[SoArealData::omegaz][i]);

                          RealVect p2omega(p_realarray[SoArealData::omegax][j],
                                           p_realarray[SoArealData::omegay][j],
                                           p_realarray[SoArealData::omegaz][j]);

                          cfrelvel(p1vel, p2vel, p1radius, p2radius, p1omega,
                              p2omega, vrel_trans_norm, vrel_t, normal, dist_mag);

                          int phase1 = p_intarray[SoAintData::phase][i];
                          int phase2 = p_intarray[SoAintData::phase][j];

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

                          // calculate the tangential overlap
                          overlap_t[0] = subdt*vrel_t[0];
                          overlap_t[1] = subdt*vrel_t[1];
                          overlap_t[2] = subdt*vrel_t[2];
                          mag_overlap_t = sqrt(dot_product(overlap_t, overlap_t));

                          if (mag_overlap_t > 0.0) {
                              Real fnmd = local_mew * sqrt(dot_product(fn, fn));
                              RealVect tangent(0.);
                              tangent[0] = overlap_t[0]/mag_overlap_t;
                              tangent[1] = overlap_t[1]/mag_overlap_t;
                              tangent[2] = overlap_t[2]/mag_overlap_t;
                              ft[0] = -fnmd * tangent[0];
                              ft[1] = -fnmd * tangent[1];
                              ft[2] = -fnmd * tangent[2];
                          } else {
                              ft[0] = 0.0;
                              ft[1] = 0.0;
                              ft[2] = 0.0;
                          }


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

                          Real dist_cl1 = 0.5 * (dist_mag + (p1radius*p1radius - p2radius*p2radius) * dist_mag_inv);
                          dist_cl1 = dist_mag - dist_cl1;

                          Real dist_cl2 = 0.5 * (dist_mag + (p2radius*p2radius - p1radius*p1radius) * dist_mag_inv);
                          dist_cl2 = dist_mag - dist_cl2;

                          RealVect tow_force(0.);

                          cross_product(normal, ft, tow_force);

#ifdef _OPENMP
#pragma omp critical
                          {
#endif
                            Gpu::Atomic::Add(&tow_ptr[i         ], dist_cl1*tow_force[0]);
                            Gpu::Atomic::Add(&tow_ptr[i + ntot  ], dist_cl1*tow_force[1]);
                            Gpu::Atomic::Add(&tow_ptr[i + 2*ntot], dist_cl1*tow_force[2]);

                            if (j < nrp)
                            {
                                Gpu::Atomic::Add(&tow_ptr[j         ], dist_cl2*tow_force[0]);
                                Gpu::Atomic::Add(&tow_ptr[j + ntot  ], dist_cl2*tow_force[1]);
                                Gpu::Atomic::Add(&tow_ptr[j + 2*ntot], dist_cl2*tow_force[2]);
                            }
#ifdef _OPENMP
                          }
#endif
                      }
                  }
              });

            Gpu::Device::synchronize();

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

            amrex::ParallelFor(nrp,
              [pstruct,p_realarray,subdt,fc_ptr,ntot,gravity,tow_ptr,eps,p_hi,p_lo,
               x_lo_bc,x_hi_bc,y_lo_bc,y_hi_bc,z_lo_bc,z_hi_bc]
              AMREX_GPU_DEVICE (int i) noexcept
              {
                auto& particle = pstruct[i];

                RealVect ppos(particle.pos());

                Real mass = p_realarray[SoArealData::mass][i];

                p_realarray[SoArealData::velx][i] += subdt * (
                    (p_realarray[SoArealData::dragx][i] + fc_ptr[i]) /
                     mass + gravity[0]
                );
                p_realarray[SoArealData::vely][i] += subdt * (
                    (p_realarray[SoArealData::dragy][i] + fc_ptr[i+ntot]) /
                     mass + gravity[1]
                );
                p_realarray[SoArealData::velz][i] += subdt * (
                    (p_realarray[SoArealData::dragz][i] + fc_ptr[i+2*ntot]) /
                     mass + gravity[2]
                );

                p_realarray[SoArealData::omegax][i] +=
                  subdt * p_realarray[SoArealData::oneOverI][i] * tow_ptr[i];
                p_realarray[SoArealData::omegay][i] +=
                  subdt * p_realarray[SoArealData::oneOverI][i] * tow_ptr[i+ntot];
                p_realarray[SoArealData::omegaz][i] +=
                  subdt * p_realarray[SoArealData::oneOverI][i] * tow_ptr[i+2*ntot];

                ppos[0] += subdt * p_realarray[SoArealData::velx][i];
                ppos[1] += subdt * p_realarray[SoArealData::vely][i];
                ppos[2] += subdt * p_realarray[SoArealData::velz][i];

                particle.pos(0) = ppos[0];
                particle.pos(1) = ppos[1];
                particle.pos(2) = ppos[2];
              });

            BL_PROFILE_VAR_STOP(des_time_march);

            Gpu::synchronize();

            usr2_des(nrp, ptile);

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

        // Update substep count
        n += 1;

        /************************************************************************
         * DEBUG: output the number of collisions in current substep            *
         *        output the max velocity (and forces) in current substep       *
         *        update max velocities and forces                              *
         ***********************************************************************/

        if (debug_level > 0) ncoll_total += ncoll;

        if (debug_level > 1) {
            ParallelDescriptor::ReduceIntSum(ncoll, ParallelDescriptor::IOProcessorNumber());
            Print() << "Number of collisions: " << ncoll << " at step " << n << std::endl;
        }

        if (debug_level > 0){
            UpdateMaxVelocity();
            UpdateMaxForces(pfor, wfor);
        }

        if (debug_level > 1) {
            RealVect max_vel = GetMaxVelocity();
            Vector<RealVect> max_forces = GetMaxForces();

            const Real * dx_crse = Geom(0).CellSize();
            amrex::Print() << "Maximum distance traveled:"
                           << std::endl
                           <<  "x= " << max_vel[0] * dt
                           << " y= " << max_vel[1] * dt
                           << " z= " << max_vel[2] * dt
                           << " and note that "
                           << " dx= " << dx_crse[0] << std::endl;

            amrex::Print() << "Maximum particle-particle (pp) and particle-wall (pw) forces:"
                           << std::endl
                           <<  "ppx= " << max_forces[0][0]
                           << " ppy= " << max_forces[0][1]
                           << " ppz= " << max_forces[0][2] << std::endl
                           <<  "pwx= " << max_forces[1][0]
                           << " pwy= " << max_forces[1][1]
                           << " pwz= " << max_forces[1][2] << std::endl;

        }

    } // end of loop over substeps

    // Redistribute particles at the end of all substeps (note that the particle
    // neighbour list needs to be reset when redistributing).
    clearNeighbors();
    Redistribute(0, 0, 0, 1);

    /****************************************************************************
     * DEBUG: output the total number of collisions over all substeps           *
     *        output the maximum velocity and forces over all substeps          *
     ***************************************************************************/
    if (debug_level > 0) {
        ParallelDescriptor::ReduceIntSum(ncoll_total, ParallelDescriptor::IOProcessorNumber());
        amrex::Print() << "Number of collisions: " << ncoll_total << " in "
                       << nsubsteps << " substeps " << std::endl;
    }

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (BMXParIter pti(*this, lev); pti.isValid(); ++pti)
    {
        const int nrp   = NumberOfParticles(pti);
        void* particles = pti.GetArrayOfStructs().data();

        usr3_des(nrp,particles);
    }

    if (debug_level > 0) {
        RealVect max_vel = GetMaxVelocity();
        Vector<RealVect> max_forces = GetMaxForces();

        const Real * dx_crse = Geom(0).CellSize();
        amrex::Print() << "Maximum possible distance traveled:" << std::endl
                       <<  "x= " << max_vel[0] * dt
                       << " y= " << max_vel[1] * dt
                       << " z= " << max_vel[2] * dt
                       << " and note that "
                       << " dx= " << dx_crse[0] << std::endl;

        amrex::Print() << "Maximum particle-particle (pp) and particle-wall (pw) forces:" << std::endl
                       <<  "ppx= " << max_forces[0][0]
                       << " ppy= " << max_forces[0][1]
                       << " ppz= " << max_forces[0][2] << std::endl
                       <<  "pwx= " << max_forces[1][0]
                       << " pwy= " << max_forces[1][1]
                       << " pwz= " << max_forces[1][2] << std::endl;
    }

    amrex::Print() << "done. \n";

    BL_PROFILE_REGION_STOP("bmx_dem::EvolveParticles()");
}
