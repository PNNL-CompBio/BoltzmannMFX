//
//     Copyright (c) 2013 Battelle Memorial Institute
//     Licensed under modified BSD License. A copy of this license can be found
//     in the LICENSE file in the top level directory of this distribution.
//
#include <bmx_pc.H>
#include <bmx_dem_parms.H>
#include <bmx_bc_parms.H>
#include <bmx_chem_K.H>

using namespace amrex;

void BMXParticleContainer::ParticleExchange (Real dt,
                                            const Vector<MultiFab*> cost,
                                            std::string& knapsack_weight_type,
                                            int& nsubsteps)
{
    BL_PROFILE_REGION_START("bmx_dem::ParticleExchange()");
    BL_PROFILE("bmx_dem::ParticleExchange()");

    Real eps = std::numeric_limits<Real>::epsilon();

    for (int lev = 0; lev <= finest_level; lev++)
    {

    int n_at_lev = this->NumberOfParticlesAtLevel(lev);
    amrex::Print() << "In Particle Exchange with " << n_at_lev << " particles at level " << lev << std::endl;

    if (n_at_lev == 0) continue;

    amrex::Print() << " particles on level: " << lev
                   << " ... with fluid dt " << dt << std::endl;

    BMXChemistry *chemistry = BMXChemistry::instance();
    amrex::Gpu::DeviceVector<Real> xpar_vec = chemistry->getExchangeParameters();
    //Real *xpar = &xpar_vec[0];
    auto xpar = xpar_vec.data();

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

            /********************************************************************
             * Particle-Particle collision forces (and torques)                 *
             *******************************************************************/

            BL_PROFILE_VAR("calc_particle_exchanges()", calc_particle_exchanges);

            auto nbor_data = m_neighbor_list[lev][index].data();


            // now we loop over the neighbor list and compute exchanges
            int me = ParallelDescriptor::MyProc();
            amrex::ParallelFor(nrp,
                [nrp,pstruct,nbor_data,subdt,ntot,xpar]
              AMREX_GPU_DEVICE (int i) noexcept
              {
              auto& particle = pstruct[i];
              if (particle.idata(intIdx::position) == siteLocation::TIP &&
                  particle.idata(intIdx::n_bnds) > 2) {
                int *idata = &particle.idata(0);
                std::printf("particle id: %d cpu: %d is Tip with %d bonds\n",
                    idata[intIdx::id],idata[intIdx::cpu],idata[intIdx::n_bnds]);
              }
#if 0
              std::printf("xpar[0]: %f\n",xpar[0]);
              std::printf("xpar[1]: %f\n",xpar[1]);
              std::printf("xpar[2]: %f\n",xpar[2]);
#endif


                  // initialize particle before evaluating exchange
                  initExchange(&particle.rdata(0),&particle.idata(0));

                  // Initialize increments to zero
                  const auto neighbs = nbor_data.getNeighbors(i);
                  for (auto mit = neighbs.begin(); mit != neighbs.end(); ++mit)
                  {
                      auto p2 = *mit;
                      const int j = mit.index();

                      AMREX_ASSERT_WITH_MESSAGE(
                          not (particle.id() == p2.id() and
                            particle.cpu() == p2.cpu()),
                          "A particle should not be its own neighbor!");

                      evaluateExchange(&particle.rdata(0),&p2.rdata(0),
                          &particle.idata(0), &p2.idata(0), xpar, subdt);

                      // TODO: Do we need an OPENMP pragma here?

                  } // end of neighbor loop
              }); // end of loop over particles

            amrex::Gpu::Device::synchronize();

            BL_PROFILE_VAR_STOP(calc_particle_exchanges);


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

    for (int lev = 0; lev <= finest_level; lev++)
    {

    int n_at_lev = this->NumberOfParticlesAtLevel(lev);
    if (n_at_lev == 0) continue;

//    BMXChemistry *chemistry = BMXChemistry::instance();
//    std::vector<Real> xpar_vec = chemistry->getExchangeParameters();
//    Real *xpar = &xpar_vec[0];

    // Debug level controls the detail of debug output:
    //   -> debug_level = 0 : no debug output
    //   -> debug_level = 1 : debug output for every fluid step
    //   -> debug_level = 2 : debug output for every substep
    const int debug_level = 0;

    /********************************************************************
     * Particles routines                                               *
     *******************************************************************/
     BL_PROFILE_VAR("des::update_particle_concentrations()", des_time_march);
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


      /********************************************************************
       * Update particles based on calculated chemical increments         *
       *******************************************************************/

      // now we loop over the neighbor list and compute exchanges
      amrex::ParallelFor(nrp,
          [nrp,pstruct]
          AMREX_GPU_DEVICE (int i) noexcept
          {
          auto& particle = pstruct[i];

          applyIncrements(&particle.rdata(0),&particle.idata(0));
          }); // end of loop over particles

      amrex::Gpu::Device::synchronize();
    }
    BL_PROFILE_VAR_STOP(des_time_march);
    }


    BL_PROFILE_REGION_STOP("bmx_dem::ParticleExchange()");
}
