#include <bmx.H>
#include <bmx_des_K.H>

/**
 * @brief this function splits particles if some criterion is met
 */
void
BMXParticleContainer::split_particles ()
{
  BMXChemistry *bmxchem = BMXChemistry::instance();

  for (int lev = 0; lev <= finest_level; lev++) {

      for (BMXParIter pti(*this, lev); pti.isValid(); ++pti)
      {
        auto& particles = pti.GetArrayOfStructs();
        BMXParticleContainer::ParticleType* pstruct = particles().dataPtr();

        const int grid = pti.index();
        const int tile = pti.LocalTileIndex();
        auto& particle_tile = this->GetParticles(lev)[std::make_pair(grid,tile)];

        const int np = particles.size();

        // count how many particles will be split
        Gpu::DeviceVector<unsigned int> counts(np+1, 0);
        auto pcounts = counts.data();
        amrex::ParallelFor( np, [=] AMREX_GPU_DEVICE (int i) noexcept
        {
            BMXParticleContainer::ParticleType& p = pstruct[i];
            if (bmxchem->checkSplit(&p.rdata(0), &p.rdata(realIdx::first_data)))
            {
                pcounts[i] = 1;
            }
        });
        Gpu::Device::synchronize();

        // Prefix sum to count total number of new particles to create
        Gpu::DeviceVector<unsigned int> offsets(np);
        Gpu::exclusive_scan(counts.begin(), counts.end(), offsets.begin());
        unsigned int num_split;
#ifdef AMREX_USE_GPU
        Gpu::dtoh_memcpy(&num_split,offsets.dataPtr()+np,sizeof(unsigned int));
#else
        std::memcpy(&num_split,offsets.dataPtr()+np,sizeof(unsigned int));
#endif

        // make room for new particles
        particle_tile.resize(np+num_split);

        // Update NextID to include particles created in this function
        Long next_pid;
#ifdef AMREX_USE_OMP
#pragma omp critical (add_plasma_nextid)
#endif
        {
            next_pid = ParticleType::NextID();
            ParticleType::NextID(next_pid+num_split);
        }

        // Fill new particle data. If particle pid is split, the new particle
        // is at index np + poffsets[pid];
        auto poffsets = offsets.data();
        pstruct = particles().dataPtr();
        amrex::ParallelFor( np, [=] AMREX_GPU_DEVICE (int pid) noexcept
        {
            BMXParticleContainer::ParticleType& p_orig = pstruct[pid];

            // Check to see if particle satisfies some criteria for splitting
            // into two new particles
            Real* p_par =  &p_orig.rdata(0);
            //            std::printf("TESTING PARTICLE WITH VOL = %d %f \n", (Long) p_orig.id(), p_par[realIdx::vol]);
            if (pcounts[pid])
            {
                std::printf("MAKING NEW PARTICLE FROM VOL = %f \n", p_par[realIdx::vol]);
                ParticleType& p = pstruct[np+poffsets[pid]];
                p.id()  = next_pid + poffsets[pid];
                p.cpu() = amrex::ParallelDescriptor::MyProc();

                Real *pos_orig = &p_orig.pos(0);
                Real *pos_new  = &p.pos(0);

                Real *par_orig = &p_orig.rdata(0);
                Real *par_new  = &p.rdata(0);

                int *ipar_orig = &p_orig.idata(0);
                int *ipar_new  = &p.idata(0);

                // Set parameters on new particle base on values from
                // original particle
                bmxchem->setNewCell(pos_orig, pos_new, par_orig,
                                    par_new, ipar_orig, ipar_new);

                //                std::printf("OLD VOLUME OUT = %f \n", par_orig[realIdx::vol]);
                //                std::printf("NEW VOLUME OUT = %f \n", par_new[realIdx::vol]);
            } // if test
        }); // pid
      } // pti
  } // lev

  // HACK HACK -- THIS LOOP IS JUST FOR DEBUGGING!
  for (int lev = 0; lev <= finest_level; lev++) {

      for (BMXParIter pti(*this, lev); pti.isValid(); ++pti)
      {
        auto& particles = pti.GetArrayOfStructs();
        BMXParticleContainer::ParticleType* pstruct = particles().dataPtr();

        const int grid = pti.index();
        const int tile = pti.LocalTileIndex(); 
        auto& particle_tile = this->GetParticles(lev)[std::make_pair(grid,tile)];

        const int np = particles.size();

        for (int pid = 0; pid < np; ++pid)     
        {
              BMXParticleContainer::ParticleType& p_orig = pstruct[pid];

              Real* p_par =  &p_orig.rdata(0);
            //  std::cout << "PRE_REDIST:TESTING PARTICLE WITH VOL = " << p_orig.id() << " " << p_par[realIdx::vol] <<  std::endl;
        } // pid
      } // pti
  } // lev

  // Redistribute the particles so the new particles end up on the right processor
  Redistribute();

  // HACK HACK -- THIS LOOP IS JUST FOR DEBUGGING!
  for (int lev = 0; lev <= finest_level; lev++) {

      for (BMXParIter pti(*this, lev); pti.isValid(); ++pti)
      {
        auto& particles = pti.GetArrayOfStructs();
        BMXParticleContainer::ParticleType* pstruct = particles().dataPtr();

        const int grid = pti.index();
        const int tile = pti.LocalTileIndex(); 
        auto& particle_tile = this->GetParticles(lev)[std::make_pair(grid,tile)];

        const int np = particles.size();

        for (int pid = 0; pid < np; ++pid)     
        {
              BMXParticleContainer::ParticleType& p_orig = pstruct[pid];

              Real* p_par =  &p_orig.rdata(0);
            //  std::cout << "POST_REDIST:TESTING PARTICLE WITH VOL = " << p_orig.id() << " " << p_par[realIdx::vol] <<  std::endl;
        } // pid
      } // pti
  } // lev
}
