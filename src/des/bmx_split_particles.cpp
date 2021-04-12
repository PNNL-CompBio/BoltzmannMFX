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

/**
 * @brief this function splits particles if some criterion is met
 */
void 
BMXParticleContainer::split_particles ()
{
  BMXChemistry *bmxchem = BMXChemistry::instance();
  Real PI = 4.0*atan(1.0);

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

        //
        // Note: this will happen on CPU only -- we will need to add something to count how many
        //       new particles, resize the particle_tile appropriately, then fill the data for the 
        //       new particles on the GPU
        //
        for (int pid = 0; pid < np; ++pid)     
        {
              BMXParticleContainer::ParticleType& p_orig = pstruct[pid];

              // This is just a made-up test on particle volume to trip splitting
              if (p_orig.rdata(realIdx::vol) > FLUID::max_vol)
              {
                   ParticleType p;
                   p.id()  = ParticleType::NextID();
                   p.cpu() = amrex::ParallelDescriptor::MyProc();
                   Real x,y,z;
                   x = p_orig.pos(0);
                   y = p_orig.pos(1);
                   z = p_orig.pos(2);

                   // Copy values from original particle and modify some values
                   // as appropriate
                   bmxchem->setChildParameters(&p_orig.rdata(0), &p_orig.idata(0),
                       &p.rdata(0),&p.idata(0));
                   // Find new locations for split particles
                   Real radius = p.rdata(realIdx::a_size);
                   Real theta = PI * amrex::Random();
                   Real phi = 2.0 * PI * amrex::Random();
                   Real nx = cos(theta)*cos(phi);
                   Real ny = cos(theta)*sin(phi);
                   Real nz = sin(theta);
                   p.pos(0) = x + 0.5*nx*radius;
                   p.pos(1) = y + 0.5*ny*radius;
                   p.pos(2) = z + 0.5*nz*radius;
                   p_orig.pos(0) = x - 0.5*nx*radius;
                   p_orig.pos(1) = y - 0.5*ny*radius;
                   p_orig.pos(2) = z - 0.5*nz*radius;

                   particle_tile.push_back(p);
              } // if test
            } // pid
      } // pti
  } // lev
}
