//
//     Copyright (c) 2013 Battelle Memorial Institute
//     Licensed under modified BSD License. A copy of this license can be found
//     in the LICENSE file in the top level directory of this distribution.
//
#include <bmx_pc.H>
#include <bmx_dem_parms.H>
#include <bmx_bc_parms.H>
#include <bmx_cell_interaction_K.H>
#include <bmx_chem_K.H>

using namespace amrex;

/*******************************************************************************
 *  Check all tips to see if they fuse with a neighbor.                        *
 ******************************************************************************/
bool BMXParticleContainer::EvaluateTipFusion (const Vector<MultiFab*> cost,
                                              std::string& knapsack_weight_type)
{
    bool ret = false;
    BL_PROFILE_REGION_START("bmx_dem::EvaluateTipFusion()");
    BL_PROFILE("bmx_dem::EvaluateTipFusion()");

    int global_fused = 0;
    for (int lev = 0; lev <= finest_level; lev++)
    {

    int n_at_lev = this->NumberOfParticlesAtLevel(lev);

    if (n_at_lev == 0) continue;

    BMXChemistry *chemistry = BMXChemistry::instance();
    std::vector<Real> fpar_vec = chemistry->getFusionParameters();
    Real *fpar = &fpar_vec[0];
    BMXCellInteraction *interaction = BMXCellInteraction::instance();
    std::vector<Real> xpar_vec = interaction->getForceParams();
    Real *xpar = &xpar_vec[0];
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
    // const int debug_level = 0;

    // Don't redistribute particles since that gets done elsewhere but update
    // the neighbour list with fresh data
    //if (n % 25 == 0) {
      clearNeighbors();
      Redistribute(0, 0, 0, 1);
      fillNeighbors();
      // send in "false" for sort_neighbor_list option

      buildNeighborList(BMXCheckPair(DEM::neighborhood, false), false);
    // } else {
    //   updateNeighbors();
    // }

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

      Gpu::DeviceScalar<int> fused_gpu(0);
      int* fused = fused_gpu.dataPtr();

      const int nrp = GetParticles(lev)[index].numRealParticles();
      // Number of particles including neighbor particles

      /********************************************************************
       * Particle-Particle collision forces (and torques)                 *
       *******************************************************************/

      BL_PROFILE_VAR("calc_tip_fusions()", calc_tip_fusions);

      auto nbor_data = m_neighbor_list[lev][index].data();

      constexpr Real small_number = 1.0e-15;

      // now we loop over the neighbor list and compute the forces
      int me = ParallelDescriptor::MyProc();
      amrex::ParallelForRNG(nrp,
          [nrp,pstruct,nbor_data,fpar,xpar,fused,me]
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
          int fusing;

          Real dist_x = pos1[0] - p2.pos(0);
          Real dist_y = pos1[1] - p2.pos(1);
          Real dist_z = pos1[2] - p2.pos(2);

          Real r2 = dist_x*dist_x +
            dist_y*dist_y +
            dist_z*dist_z;

          RealVect diff(dist_x,dist_y,dist_z);

          Real r_lm = maxInteractionDistance(&particle.rdata(0),&p2.rdata(0),
              &particle.idata(0),&p2.idata(0),xpar[0]);

          AMREX_ASSERT_WITH_MESSAGE(
              not (particle.id() == p2.id() and
                particle.cpu() == p2.cpu()),
              "A particle should not be its own neighbor!");

          if ( r2 <= (r_lm - small_number)*(r_lm - small_number) )
          {

            // Check to see if particle is a TIP. If it is, then decide
            // whether it is fusing to p2. If fusion occurs, mark particle
            // as being fused to p2.
            checkTipFusion(&diff[0], &particle.rdata(0), &particle.idata(0),
                &p2.rdata(0), &p2.idata(0), fpar, me, &fusing, engine);

            Gpu::Atomic::Max(fused, fusing);
            // TODO: Do we need an OPENMP pragma here?

          }
          } // end of neighbor loop
          }); // end of loop over particles

      amrex::Gpu::Device::synchronize();

      global_fused = std::max(global_fused, fused_gpu.dataValue());


      BL_PROFILE_VAR_STOP(calc_tip_fusions);

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
    ParallelDescriptor::ReduceIntMax(global_fused);

    // Redistribute particles at the end of all substeps (note that the particle
    // neighbour list needs to be reset when redistributing).
    //clearNeighbors();
    //Redistribute(0, 0, 0, 1);
    //updateNeighbors();

    } // lev
    if (global_fused == 1) ret = true;

    if (!ret) Redistribute(0, 0, 0, 1);
    BL_PROFILE_REGION_STOP("bmx_dem::EvaluateTipFusion()");
    return ret;
}

/*******************************************************************************
 *  Check all interior segments to see if they have fused to a tip. Split the  *
 *  segment if fusion has occured                                              *
 ******************************************************************************/
void BMXParticleContainer::EvaluateInteriorFusion (const Vector<MultiFab*> cost,
                                              std::string& knapsack_weight_type)
{
  BL_PROFILE_REGION_START("bmx_dem::EvaluateInteriorFusion()");
  BL_PROFILE("bmx_dem::EvaluateInteriorFusion()");

  for (int lev = 0; lev <= finest_level; lev++)
  {

    int n_at_lev = this->NumberOfParticlesAtLevel(lev);

    if (n_at_lev == 0) continue;

    BMXChemistry *chemistry = BMXChemistry::instance();
    std::vector<Real> fpar_vec = chemistry->getFusionParameters();
    Real *fpar = &fpar_vec[0];
    BMXCellInteraction *interaction = BMXCellInteraction::instance();
    std::vector<Real> xpar_vec = interaction->getForceParams();
    Real *xpar = &xpar_vec[0];
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
    // const int debug_level = 0;

    // Don't redistribute particles since that gets done elsewhere but update
    // the neighbour list with fresh data
#if 1
      clearNeighbors();
      //Redistribute(0, 0, 0, 1);
      fillNeighbors();
      // send in "false" for sort_neighbor_list option

      buildNeighborList(BMXCheckPair(DEM::neighborhood, false), false);
#else
      updateNeighbors();
#endif

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
      BL_PROFILE_VAR("calc_interior_fusions()", calc_interior_fusions);

      //const Box& bx = pti.tilebox(); // UNUSED_VARIABLE
      PairIndex index(pti.index(), pti.LocalTileIndex());

      auto& plev = GetParticles(lev);
      auto& ptile = plev[index];
      auto& particles  = ptile.GetArrayOfStructs();
      ParticleType* pstruct = particles().dataPtr();

      const int grid = pti.index();
      const int tile = pti.LocalTileIndex();
      auto& particle_tile = this->GetParticles(lev)[std::make_pair(grid,tile)];

      const int nrp = GetParticles(lev)[index].numRealParticles();
      const int num_total = GetParticles(lev)[index].numTotalParticles();

      const int np = particles.size();

      /********************************************************************
       * Particle-Particle collision forces (and torques)                 *
       *******************************************************************/


      auto nbor_data = m_neighbor_list[lev][index].data();

      constexpr Real small_number = 1.0e-15;

      // now we loop over the neighbor list and compute the forces
      int me = ParallelDescriptor::MyProc();
      Gpu::DeviceVector<unsigned int> do_split(nrp+1, 0);
      auto do_split_p = do_split.data();
      bool did_fusion = false;
      amrex::ParallelFor(nrp,
          [nrp,pstruct,nbor_data,fpar,xpar,me,do_split_p,did_fusion]
          AMREX_GPU_DEVICE (int i) noexcept
          {
          auto& particle = pstruct[i];

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

          RealVect diff(dist_x,dist_y,dist_z);

          Real r_lm = maxInteractionDistance(&particle.rdata(0),&p2.rdata(0),
              &particle.idata(0),&p2.idata(0),xpar[0]);
          AMREX_ASSERT_WITH_MESSAGE(
              not (particle.id() == p2.id() and
                particle.cpu() == p2.cpu()),
              "A particle should not be its own neighbor!");

          int split_flag = 0;
          if ( r2 <= (r_lm - small_number)*(r_lm - small_number) )
          {

            // Check to see if p2 is fused to particle. If it is,
            // then mark particle as being fused to p2 and set
            // split_flag to 1.
            checkInteriorFusion(&particle.rdata(0), &particle.idata(0),
                &p2.rdata(0), &p2.idata(0), &split_flag, me);
            if (split_flag == 1) do_split_p[i]++;
            if (split_flag == 1) {
              printf("p[%d] Splitting particle id: %d cpu: %d\n",me,
                  particle.idata(intIdx::id),particle.idata(intIdx::cpu));
            }
            

            // TODO: Do we need an OPENMP pragma here?

          }
          } // end of neighbor loop
          }); // end of loop over particles

      amrex::Gpu::Device::synchronize();

      // Prefix sum to count total number of new particles to create
      Gpu::DeviceVector<unsigned int> offsets(nrp+1);
      Gpu::exclusive_scan(do_split.begin(), do_split.end(), offsets.begin());
      unsigned int num_split;
#ifdef AMREX_USE_GPU
      Gpu::dtoh_memcpy(&num_split,offsets.dataPtr()+nrp,sizeof(unsigned
            int));
#else
      std::memcpy(&num_split,offsets.dataPtr()+nrp,sizeof(unsigned
            int));
#endif

      // make room for new particles - invalidates iterators, so get the
      // ptr again
      particle_tile.resize(num_total+num_split);
//      particle_tile.setNumNeighbors(0);
      pstruct = particles().dataPtr();
      // Update NextID to include particles created in this function
      Long next_pid;
#ifdef AMREX_OPENMP
#pragma omp critical (add_plasma_nextid)
#endif
      {
        next_pid = ParticleType::NextID();
        ParticleType::NextID(next_pid+num_split);
      }
      // Fill new particle data. If particle pid is split, the new particle
      // is at index np + poffsets[pid]


      int l_num_reals = BMXChemistry::p_num_reals;
      int l_num_ints  = BMXChemistry::p_num_ints;

      auto poffsets = offsets.data();
      int my_proc = amrex::ParallelDescriptor::MyProc();
      amrex::ParallelFor( nrp, [=] AMREX_GPU_DEVICE (int pid) noexcept
          {
          BMXParticleContainer::ParticleType& p_orig = pstruct[pid];
          // Check to see if particle is splitting
          // into two new particles
          if (do_split_p[pid] == 1) {
            ParticleType& p = pstruct[nrp+poffsets[pid]];
            p.id()  = next_pid + poffsets[pid];
            p.cpu() = my_proc;

            Real *pos_orig = &p_orig.pos(0);
            Real *pos_new  = &p.pos(0);

            Real *par_orig = &p_orig.rdata(0);
            Real *par_new  = &p.rdata(0);

            int *ipar_orig = &p_orig.idata(0);
            int *ipar_new  = &p.idata(0);

            // Set parameters on new particled base on values from
            // original particle
            setSplitSegment(pos_orig, pos_new, par_orig,
                par_new, ipar_orig, ipar_new,
                l_num_reals, l_num_ints, p.id(), p.cpu());
            ipar_new[intIdx::id] = p.id();
            ipar_new[intIdx::cpu] = p.cpu();
            ipar_orig[intIdx::fuse_flag] = 0;
            ipar_orig[intIdx::split_flag] = 1;
            ipar_orig[intIdx::new_flag] = 0;
            ipar_orig[intIdx::fuse_id] = -1;
            ipar_orig[intIdx::fuse_cpu] = -1;
            ipar_new[intIdx::fuse_flag] = 0;
            ipar_new[intIdx::split_flag] = 0;
            ipar_new[intIdx::new_flag] = 1;
            ipar_new[intIdx::fuse_id] = -1;
            ipar_new[intIdx::fuse_cpu] = -1;
          } else if (do_split_p[pid] > 1) {
            //TODO: Simultaneous fusion happened. We don't know how to
            //handle this.
#ifndef AMREX_USE_GPU
            std::cout<<"Simultaneous fusion event happened. We cannot"
              " handle this situation"<<std::endl;
#endif
          } else {
            int *ipar_orig = &p_orig.idata(0);
            ipar_orig[intIdx::split_flag] = 0;
          } // if test
          }); // pid
          BL_PROFILE_VAR_STOP(calc_interior_fusions);

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
    } // pti


  } // lev

  // Redistribute particles at the end of all substeps (note that the particle
  // neighbour list needs to be reset when redistributing).
  clearNeighbors();
  Redistribute(0, 0, 0, 1);

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif

  BL_PROFILE_REGION_STOP("bmx_dem::EvaluateInteriorFusion()");
}

/*******************************************************************************
 *  Clean up bonding information if an interior split has occured. Additional  *
 *  cleanup occurs in force loop                                               *
 ******************************************************************************/
void BMXParticleContainer::CleanupFusion (const Vector<MultiFab*> cost,
                                              std::string& knapsack_weight_type)
{
#if 1
  BL_PROFILE_REGION_START("bmx_dem::CleanupFusion()");
  BL_PROFILE("bmx_dem::CleanupFusion()");

  /*
  BMXCellInteraction *interaction = BMXCellInteraction::instance();
  std::vector<Real> xpar_vec = interaction->getForceParams();
  Real *xpar = &xpar_vec[0];
  */
  BMXChemistry *bmxchem = BMXChemistry::instance();
  std::vector<Real> chempar_vec;
  bmxchem->getChemParams(chempar_vec);
  Real *chempar = &chempar_vec[0];
  std::vector<Real> fpar_vec = bmxchem->getFusionParameters();
  Real *fpar = &fpar_vec[0];

  for (int lev = 0; lev <= finest_level; lev++)
  {

    int n_at_lev = this->NumberOfParticlesAtLevel(lev);

    if (n_at_lev == 0) continue;

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
    // const int debug_level = 0;

    // Don't redistribute particles since that gets done elsewhere but update
    // the neighbour list with fresh data
#if 1
      clearNeighbors();
      //Redistribute(0, 0, 0, 1);
      //Redistribute(0, finest_level, 0, 1);
      Redistribute();
      fillNeighbors();
      // send in "false" for sort_neighbor_list option

      buildNeighborList(BMXCheckPair(DEM::neighborhood, false), false);
#else
      updateNeighbors();
#endif

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
      BL_PROFILE_VAR("cleanup_fusions()", cleanup_fusions);

      //const Box& bx = pti.tilebox(); // UNUSED_VARIABLE
      PairIndex index(pti.index(), pti.LocalTileIndex());

      auto& plev = GetParticles(lev);
      auto& ptile = plev[index];
      auto& particles  = ptile.GetArrayOfStructs();
      ParticleType* pstruct = particles().dataPtr();

      const int grid = pti.index();
      const int tile = pti.LocalTileIndex();
      auto& particle_tile = this->GetParticles(lev)[std::make_pair(grid,tile)];

      const int nrp = GetParticles(lev)[index].numRealParticles();
      const int num_total = GetParticles(lev)[index].numTotalParticles();

       // Number of particles including neighbor particles
       int ntot = nrp;

      /********************************************************************
       * Particle-Particle collision forces (and torques)                 *
       *******************************************************************/

      auto nbor_data = m_neighbor_list[lev][index].data();

      constexpr Real small_number = 1.0e-15;

      // now we loop over the neighbor list and look for invalid connections
      int me = ParallelDescriptor::MyProc();
      amrex::ParallelFor(nrp,
          [nrp,pstruct,nbor_data,chempar,fpar,ntot,me]
          AMREX_GPU_DEVICE (int i) noexcept
          {
          auto& particle = pstruct[i];

          // increment bond scaling parameter
          incrementBondScale(&particle.rdata(0),&particle.idata(0),fpar);
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

          RealVect diff(dist_x,dist_y,dist_z);

          /*
          Real r_lm = maxInteractionDistance(&particle.rdata(0),&p2.rdata(0),
              &particle.idata(0),&p2.idata(0),xpar[0]);
              */
          Real r_lm = 1.5*chempar[9];
          AMREX_ASSERT_WITH_MESSAGE(
              not (particle.id() == p2.id() and
                particle.cpu() == p2.cpu()),
              "A particle should not be its own neighbor!");

          // int split_flag = 0;
          if ( r2 <= (r_lm - small_number)*(r_lm - small_number) )
          {
            // Check to see if particle is bound to p2. If it is
            // and p2 has split flag set to 1 and if site on particle is
            // 1 then get rid of bond on particle
            cleanupFusionBond(&particle.idata(0),&p2.idata(0),me);

          }
          } // end of neighbor loop
          }); // end of loop over particles

      amrex::Gpu::Device::synchronize();

          BL_PROFILE_VAR_STOP(cleanup_fusions);

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
    } // pti


  } // lev
  // Redistribute particles at the end of all substeps (note that the particle
  // neighbour list needs to be reset when redistributing).
  clearNeighbors();
  Redistribute(0, 0, 0, 1);

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif

  BL_PROFILE_REGION_STOP("bmx_dem::CleanupFusion()");
#endif
}

/*******************************************************************************
 *  Print out connectivity information for all segments. Only use this for
 *  debugging
 ******************************************************************************/
void BMXParticleContainer::PrintConnectivity (const Vector<MultiFab*> cost,
                                              std::string& knapsack_weight_type)
{
  BL_PROFILE_REGION_START("bmx_dem::PrintConnectivity()");
  BL_PROFILE("bmx_dem::PrintConnectivity()");

  for (int lev = 0; lev <= finest_level; lev++)
  {

    int n_at_lev = this->NumberOfParticlesAtLevel(lev);

    if (n_at_lev == 0) continue;

    /********************************************************************
     * Particle routines                                                *
     *******************************************************************/
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (BMXParIter pti(*this, lev); pti.isValid(); ++pti)
    {
      // Timer used for load-balancing
      Real wt = ParallelDescriptor::second();
      BL_PROFILE_VAR("print_connectivity()", print_connectivity);

      //const Box& bx = pti.tilebox(); // UNUSED_VARIABLE
      PairIndex index(pti.index(), pti.LocalTileIndex());

      auto& plev = GetParticles(lev);
      auto& ptile = plev[index];
      auto& particles  = ptile.GetArrayOfStructs();
      ParticleType* pstruct = particles().dataPtr();

      const int grid = pti.index();
      const int tile = pti.LocalTileIndex();
      auto& particle_tile = this->GetParticles(lev)[std::make_pair(grid,tile)];

      const int nrp = GetParticles(lev)[index].numRealParticles();
      const int num_total = GetParticles(lev)[index].numTotalParticles();

      /********************************************************************
       * Particle-Particle collision forces (and torques)                 *
       *******************************************************************/

      // now we loop over the neighbor list and compute the forces
      int me = ParallelDescriptor::MyProc();
      amrex::Gpu::Device::synchronize();

      pstruct = particles().dataPtr();
      amrex::ParallelFor( nrp, [=] AMREX_GPU_DEVICE (int pid) noexcept
          {
          BMXParticleContainer::ParticleType& p_orig = pstruct[pid];

            int *ipar = &p_orig.idata(0);
            if (ipar[intIdx::n_bnds] == 1) {
            printf("particle id: %d cpu: %d nbnds: %d site1: %d id1: %d cpu1: %d\n",
                ipar[intIdx::id],ipar[intIdx::cpu],ipar[intIdx::n_bnds],
                ipar[intIdx::site1],ipar[intIdx::seg1_id1],ipar[intIdx::seg1_id2]);
            } else if (ipar[intIdx::n_bnds] == 2) {
            printf("particle id: %d cpu: %d nbnds: %d site1: %d id1: %d cpu1: %d"
                " site2: %d id2: %d cpu2: %d\n",
                ipar[intIdx::id],ipar[intIdx::cpu],ipar[intIdx::n_bnds],
                ipar[intIdx::site1],ipar[intIdx::seg1_id1],ipar[intIdx::seg1_id2],
                ipar[intIdx::site2],ipar[intIdx::seg2_id1],ipar[intIdx::seg2_id2]);
            } else if (ipar[intIdx::n_bnds] == 3) {
            printf("particle id: %d cpu: %d nbnds: %d site1: %d id1: %d cpu1: %d"
                " site2: %d id2: %d cpu2: %d site3: %d id3: %d cpu3: %d\n",
                ipar[intIdx::id],ipar[intIdx::cpu],ipar[intIdx::n_bnds],
                ipar[intIdx::site1],ipar[intIdx::seg1_id1],ipar[intIdx::seg1_id2],
                ipar[intIdx::site2],ipar[intIdx::seg2_id1],ipar[intIdx::seg2_id2],
                ipar[intIdx::site3],ipar[intIdx::seg3_id1],ipar[intIdx::seg3_id2]);
            } else if (ipar[intIdx::n_bnds] == 4) {
            printf("particle id: %d cpu: %d nbnds: %d site1: %d id1: %d cpu1: %d"
                " site2: %d id2: %d cpu2: %d site3: %d id3: %d cpu3: %d"
                " site4: %d id4: %d cpu4: %d\n",
                ipar[intIdx::id],ipar[intIdx::cpu],ipar[intIdx::n_bnds],
                ipar[intIdx::site1],ipar[intIdx::seg1_id1],ipar[intIdx::seg1_id2],
                ipar[intIdx::site2],ipar[intIdx::seg2_id1],ipar[intIdx::seg2_id2],
                ipar[intIdx::site3],ipar[intIdx::seg3_id1],ipar[intIdx::seg3_id2],
                ipar[intIdx::site4],ipar[intIdx::seg4_id1],ipar[intIdx::seg4_id2]);
            }


          }); // pid
          BL_PROFILE_VAR_STOP(print_connectivity);

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
    } // pti


  } // lev
  // Redistribute particles at the end of all substeps (note that the particle
  // neighbour list needs to be reset when redistributing).
  clearNeighbors();
  Redistribute(0, 0, 0, 1);

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif

  BL_PROFILE_REGION_STOP("bmx_dem::PrintConnectivity()");
}
