#include <bmx.H>
#include <bmx_des_K.H>
#include <bmx_interp_K.H>
#include <bmx_filcc.H>

#include <AMReX_BC_TYPES.H>
#include <AMReX_Box.H>
#include <AMReX_FillPatchUtil.H>
#include <bmx_mf_helpers.H>
#include <bmx_dem_parms.H>
#include <bmx_fluid_parms.H>
#include <bmx_chem.H>

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

        //
        // Note: this will happen on CPU only -- we will need to add something to count how many
        //       new particles, resize the particle_tile appropriately, then fill the data for the 
        //       new particles on the GPU
        //
        for (int pid = 0; pid < np; ++pid)     
        {
              BMXParticleContainer::ParticleType& p_orig = pstruct[pid];

              // Check to see if particle satisfies some criteria for splitting
              // into two new  particles
              if (bmxchem->checkSplit(&p_orig.rdata(0), &p_orig.rdata(realIdx::first_data)))
              {
                   ParticleType p;
                   p.id()  = ParticleType::NextID();
                   p.cpu() = amrex::ParallelDescriptor::MyProc();
                   Real *pos_orig = &p_orig.pos(0);
                   Real *pos_new = &p.pos(0);
                   Real *par_orig = &p_orig.rdata(0);
                   Real *par_new = &p.rdata(0);
                   int *ipar_orig = &p_orig.idata(0);
                   int *ipar_new = &p.idata(0);
                   // Set parameters on new particle base on values from
                   // original particle
                   bmxchem->setNewCell(pos_orig, pos_new, par_orig,
                       par_new, ipar_orig, ipar_new);

                   particle_tile.push_back(p);
              } // if test
            } // pid
      } // pti
  } // lev
}