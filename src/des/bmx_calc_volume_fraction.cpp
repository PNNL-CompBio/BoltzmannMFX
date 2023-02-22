//
//     Copyright (c) 2013 Battelle Memorial Institute
//     Licensed under modified BSD License. A copy of this license can be found
//     in the LICENSE file in the top level directory of this distribution.
//
#include <bmx.H>
#include <bmx_deposition_K.H>
#include <bmx_fluid_parms.H>
#include <AMReX_AmrParticles.H>

// Calculate the fraction of each grid cell not occupied by biological cells -- this
//   1) defines vf_n using the current particle locations
//   2) updates X_k on the grid to allow for the change in vf
void bmx::bmx_calc_volume_fraction (bool adjust_X)
{
  int start_part_comp = realIdx::vol;
  int start_mesh_comp = 0;
  int        num_comp = 1;
  const int nchem_species = FLUID::nchem_species;

  // Initialize to zero because the deposition routine will only change values
  // where there are particles (note this is the default) 
  bool zero_out_input = true;

  // Divide the quantity on the mesh after deposition so it is relative to
  // the volume of the cell (note this is the default) 
  bool vol_weight     = true;

  if (bmx::m_vf_deposition_scheme == DepositionScheme::one_to_one) 
  {
      ParticleToMesh(*pc,get_vf_new(),0,finest_level,
                     OneToOneDeposition{start_part_comp,start_mesh_comp,num_comp},
                     zero_out_input, vol_weight);

  } else if (bmx::m_vf_deposition_scheme == DepositionScheme::trilinear) {

      ParticleToMesh(*pc,get_vf_new(),0,finest_level,
                     TrilinearDeposition{start_part_comp,start_mesh_comp,num_comp},
                     zero_out_input, vol_weight);
#if 0
  } else if (bmx::m_vf_deposition_scheme == DepositionScheme::square_dpvm) {

    ParticleToMesh(*pc,get_vf(),0,finest_level,DPVMSquareDeposition());

  } else if (bmx::m_vf_deposition_scheme == DepositionScheme::true_dpvm) {

    ParticleToMesh(*pc,get_vf(),0,finest_level,TrueDPVMDeposition());

  } else if (bmx::m_vf_deposition_scheme == DepositionScheme::centroid) {

    ParticleToMesh(*pc,get_vf(),0,finest_level,CentroidDeposition());
#endif

  } else {
    amrex::Abort("Don't know this deposition_scheme for volume fraction!");
  }

  for (int lev = 0; lev <= finest_level; lev++)  
  {
    MultiFab& vf = *get_vf_new()[lev];
    Geometry& gm = Geom(lev); 

    // Now define this mf = (1 - particle_vol)
    vf.mult(-1.0, vf.nGrow());
    vf.plus( 1.0, vf.nGrow());

#if 1
    // Fix up volume fraction so that the minimum value is bounded
    auto& ld = *m_leveldata[lev];

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*ld.vf_n,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
      Box const& bx = mfi.tilebox();

      Array4<Real> const& vf_n     = ld.vf_n->array(mfi);

      ParallelFor(bx, [vf_n]
          AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          {
            vf_n(i,j,k) = amrex::max(vf_n(i,j,k),0.1);
          });
    } // mfi
#endif
 
    if (vf.min(0,0) < 0.0)
    {
       amrex::Print() << "MIN OF VF AT LEVEL " << lev  << " is " << vf.min(0,0) << std::endl;
       amrex::Abort("This is unphysical");
    }

    // Move any volume deposited outside the domain back into the domain
    // when BC is either a pressure inlet or mass inflow.
    bmx_deposition_bcs(lev, vf);

    // Fill the boundaries so we calculate the correct average
    // solids volume fraction for periodic boundaries.
    vf.FillBoundary(gm.periodicity());
  }

  // Now adjust the concentrations to account for the change in vf
  if (adjust_X)
  {
      for (int lev = 0; lev <= finest_level; lev++)
      {
          auto& ld = *m_leveldata[lev];

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
          for (MFIter mfi(*ld.X_k,TilingIfNotGPU()); mfi.isValid(); ++mfi)
          {
              Box const& bx = mfi.tilebox();
              Array4<Real      > const& X_k_n    = ld.X_k->array(mfi);
              Array4<Real const> const& vf_o     = ld.vf_o->const_array(mfi);
              Array4<Real const> const& vf_n     = ld.vf_n->const_array(mfi);

              ParallelFor(bx, nchem_species, [=]
                AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
              {
                  if (vf_n(i,j,k) != vf_o(i,j,k))
                  {
                      X_k_n(i,j,k,n) *= vf_o(i,j,k) / vf_n(i,j,k);
                  }
#if 0
                  if (vf_n(i,j,k) < 0.0) {
                    printf("Volume fraction at (%d,%d,%d) is %f\n",i,j,k,vf_n(i,j,k));
                  }
#endif
            });
        } // mfi
      } // lev
  } // if (adjust_X)
}
