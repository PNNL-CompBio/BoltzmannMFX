#include <bmx.H>
#include <bmx_chem.H>
#include <bmx_chem_K.H>
#include <bmx_des_K.H>
#include <bmx_chem_species_parms.H>

/**
 * @brief this function splits particles if some criterion is met
 */
void
BMXParticleContainer::split_particles ()
{
  BMXChemistry *bmxchem = BMXChemistry::instance();

  Real l_max_vol  = SPECIES::max_vol; 
  Real l_overlap  = BMXChemistry::p_overlap;
  int l_num_reals = BMXChemistry::p_num_reals;
  int l_num_ints  = BMXChemistry::p_num_ints;

  for (int lev = 0; lev <= finest_level; lev++) {

      for (BMXParIter pti(*this, lev); pti.isValid(); ++pti)
      {
        auto& particles = pti.GetArrayOfStructs();
        BMXParticleContainer::ParticleType* pstruct = particles().dataPtr();

        const int grid = pti.index();
        const int tile = pti.LocalTileIndex();
        auto& particle_tile = this->GetParticles(lev)[std::make_pair(grid,tile)];

        const int np = particles.size();

        // whether each particle will be split
        Gpu::DeviceVector<unsigned int> do_split(np+1, 0);
        auto do_split_p = do_split.data();
        amrex::ParallelFor( np, [=] AMREX_GPU_DEVICE (int i) noexcept
        {
            BMXParticleContainer::ParticleType& p = pstruct[i];
            if (checkSplit(&p.rdata(0), &p.rdata(realIdx::first_data), l_max_vol))
            {
                do_split_p[i] = 1;
            }
        });
        Gpu::Device::synchronize();

        // Prefix sum to count total number of new particles to create
        Gpu::DeviceVector<unsigned int> offsets(np);
        Gpu::exclusive_scan(do_split.begin(), do_split.end(), offsets.begin());
        unsigned int num_split;
#ifdef AMREX_USE_GPU
        Gpu::dtoh_memcpy(&num_split,offsets.dataPtr()+np,sizeof(unsigned int));
#else
        std::memcpy(&num_split,offsets.dataPtr()+np,sizeof(unsigned int));
#endif

        // make room for new particles - invalidates iterators, so get the ptr again
        particle_tile.resize(np+num_split);
        pstruct = particles().dataPtr();

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
        int my_proc = amrex::ParallelDescriptor::MyProc();
        amrex::ParallelFor( np, [=] AMREX_GPU_DEVICE (int pid) noexcept
        {
            BMXParticleContainer::ParticleType& p_orig = pstruct[pid];

            // Check to see if particle satisfies some criteria for splitting
            // into two new particles
            Real* p_par =  &p_orig.rdata(0);
            if (do_split_p[pid])
            {
                // std::printf("MAKING NEW PARTICLE FROM VOL = %f \n", p_par[realIdx::vol]);
                ParticleType& p = pstruct[np+poffsets[pid]];
                p.id()  = next_pid + poffsets[pid];
                p.cpu() = my_proc;

                Real *pos_orig = &p_orig.pos(0);
                Real *pos_new  = &p.pos(0);

                Real *par_orig = &p_orig.rdata(0);
                Real *par_new  = &p.rdata(0);

                int *ipar_orig = &p_orig.idata(0);
                int *ipar_new  = &p.idata(0);

                // Set parameters on new particle base on values from
                // original particle
                setNewCell(pos_orig, pos_new, par_orig,
                           par_new, ipar_orig, ipar_new,
                           l_overlap, l_num_reals, l_num_ints);

                //                std::printf("OLD VOLUME OUT = %f \n", par_orig[realIdx::vol]);
                //                std::printf("NEW VOLUME OUT = %f \n", par_new[realIdx::vol]);
            } // if test
        }); // pid
      } // pti
  } // lev

  // Redistribute the particles so the new particles end up on the right processor
  Redistribute();
}
