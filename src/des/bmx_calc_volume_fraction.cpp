#include <bmx.H>
#include <bmx_pc.H>

// This re-calculates the volume fraction within the domain
// but does not change the values outside the domain
void bmx::bmx_calc_volume_fraction (int lev, MultiFab& vf)
{
    const Geometry& gm  = Geom(lev);

    // Deposit particle volume to the grid
    pc->SolidsVolumeDeposition(lev, vf);

    // Move any volume deposited outside the domain back into the domain
    // when BC is either a pressure inlet or mass inflow.
    bmx_deposition_bcs(lev, vf);

    // Sum grid boundaries to capture any material that was deposited into
    // your grid from an adjacent grid.
    vf.SumBoundary(gm.periodicity());

    // Now define this mf = (1 - particle_vol)
    vf.mult(-1.0, vf.nGrow());
    vf.plus( 1.0, vf.nGrow());

    // Fill the boundaries so we calculate the correct average
    // solids volume fraction for periodic boundaries.
    vf.FillBoundary(gm.periodicity());
}
