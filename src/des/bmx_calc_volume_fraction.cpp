#include <bmx.H>
#include <bmx_pc.H>
#include <bmx_chem.H>
#include <bmx_deposition_K.H>
#include <AMReX_AmrParticles.H>

// This re-calculates the volume fraction within the domain
// but does not change the values outside the domain
void bmx::bmx_calc_volume_fraction ()
{
  int start_part_comp = realIdx::vol;
  int start_mesh_comp = 0;
  int        num_comp = 1;

  // Initialize to zero because the deposition routine will only change values
  // where there are particles (note this is the default) 
  bool zero_out_input = true;

  // Divide the quantity on the mesh after deposition so it is relative to
  // the volume of the cell (note this is the default) 
  bool vol_weight     = true;

  if (bmx::m_deposition_scheme == DepositionScheme::trilinear) 
  {
      ParticleToMesh(*pc,get_vf(),0,finest_level,
                     TrilinearDeposition{start_part_comp,start_mesh_comp,num_comp},
                     zero_out_input, vol_weight);
#if 0
  } else if (bmx::m_deposition_scheme == DepositionScheme::square_dpvm) {

    ParticleToMesh(*pc,get_vf(),0,finest_level,DPVMSquareDeposition());

  } else if (bmx::m_deposition_scheme == DepositionScheme::true_dpvm) {

    ParticleToMesh(*pc,get_vf(),0,finest_level,TrueDPVMDeposition());

  } else if (bmx::m_deposition_scheme == DepositionScheme::centroid) {

    ParticleToMesh(*pc,get_vf(),0,finest_level,CentroidDeposition());
#endif

  } else {
    amrex::Abort("Don't know this deposition_scheme!");
  }

  for (int lev = 0; lev <= finest_level; lev++)  
  {
    MultiFab& vf = *get_vf()[lev];
    Geometry& gm = Geom(lev); 

    // Now define this mf = (1 - particle_vol)
    vf.mult(-1.0, vf.nGrow());
    vf.plus( 1.0, vf.nGrow());

    // Move any volume deposited outside the domain back into the domain
    // when BC is either a pressure inlet or mass inflow.
    bmx_deposition_bcs(lev, vf);

    // Fill the boundaries so we calculate the correct average
    // solids volume fraction for periodic boundaries.
    vf.FillBoundary(gm.periodicity());
  }
}
