#include <bmx.H>
#include <bmx_des_K.H>
#include <bmx_interp_K.H>
#include <bmx_filcc.H>

#include <AMReX_BC_TYPES.H>
#include <AMReX_Box.H>
#include <AMReX_FillPatchUtil.H>
#include <bmx_mf_helpers.H>
#include <bmx_dem_parms.H>

/**
 * @brief this function splits particles if some criterion is met
 */
void 
BMXParticleContainer::split_particles ()
{
  for (int lev = 0; lev < nlev; lev++) {

      for (BMXParIter pti(*this, lev); pti.isValid(); ++pti)
      {
        auto& particles = pti.GetArrayOfStructs();
        BMXParticleContainer::ParticleType* pstruct = particles().dataPtr();

        const int grid = pti.index();
        const int tile = pti.LocalTileIndex(); 
        auto& particle_tile = this->GetParticles(lev)[std::make_pair(grid,tile)];

        const int np = particles.size();

        Box bx = pti.tilebox ();

        // This is just a made-up threshold
        Real max_vol = 0.5;

        //
        // Note: this will happen on CPU only -- we will need to add something to count how many
        //       new particles, resize the particle_tile appropriately, then fill the data for the 
        //       new particles on the GPU
        //
        for (int pid = 0; pid < np; ++pid)     
        {
              BMXParticleContainer::ParticleType& p_orig = pstruct[pid];

              // This is just a made-up test on particle volume to trip splitting
              if (p_orig.rdata(realIdx::vol) > max_vol)
              {
                   ParticleType p;
                   p.id()  = ParticleType::NextID();
                   p.cpu() = amrex::ParallelDescriptor::MyProc();
                   p.pos(0) = p_orig.pos(0);
                   p.pos(1) = p_orig.pos(1);
                   p.pos(2) = p_orig.pos(2);

                   // THIS IS JUST AN EXAMPLE
                   // Copy values from original particle -- probably also want to 
                   //      split some of these in half
                   p_orig.rdata(0) *= 0.5;
                   p.rdata(0)        = p_orig.rdata(0);

                   // Need to set all the rest of the data as well

                   particle_tile.push_back(p);
              } // if test
            } // pid
      } // pti
  } // lev
}
