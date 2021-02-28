#include <bmx.H>
#include <bmx_filcc.H>

#include <AMReX_FillPatchUtil.H>
#include <bmx_pc.H>
#include <bmx_mf_helpers.H>

#include <bmx_dem_parms.H>
#include <bmx_diffusion_op.H>

void bmx::bmx_calc_volume_fraction (MultiFab& ep_g)
{
  BL_PROFILE("bmx::bmx_calc_volume_fraction()");

  // Start the timers ...
  const Real strttime = ParallelDescriptor::second();

  if (DEM::solve)
  {
    // This re-calculates the volume fraction within the domain
    // but does not change the values outside the domain

    const Geometry& gm  = Geom(0);

    // Deposit particle volume to the grid
    int lev = 0;
    pc->SolidsVolumeDeposition(lev, ep_g);

    // Move any volume deposited outside the domain back into the domain
    // when BC is either a pressure inlet or mass inflow.
    bmx_deposition_bcs(0, ep_g);

    // Sum grid boundaries to capture any material that was deposited into
    // your grid from an adjacent grid.
    ep_g.SumBoundary(gm.periodicity());

    // Fill the boundaries so we calculate the correct average
    // solids volume fraction for periodic boundaries.
    ep_g.FillBoundary(gm.periodicity());

  } else {
      ep_g.setVal(0.);
  }

  if (m_verbose > 1) {
    Real stoptime = ParallelDescriptor::second() - strttime;

    ParallelDescriptor::ReduceRealMax(stoptime,ParallelDescriptor::IOProcessorNumber());

    amrex::Print() << "BMXParticleContainer::PICDeposition time: " << stoptime << '\n';
  }
}
