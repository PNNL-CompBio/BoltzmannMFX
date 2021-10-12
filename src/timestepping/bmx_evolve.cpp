#include <bmx.H>
#include <bmx_fluid_parms.H>
#include <bmx_dem_parms.H>

// This subroutine is the driver for time stepping the whole system
// (fluid + particles )
void
bmx::Evolve (int nstep, Real & dt, Real & prev_dt, Real time, Real stop_time)
{
    BL_PROFILE_REGION_START("bmx::Evolve");

    Real coupling_timing(0.);
    Real drag_timing(0.);

    /****************************************************************************
     *                                                                          *
     * Evolve Fluid and Update Chemistry inside Particles                       *
     *                                                                          *
     ***************************************************************************/
    Real start_fluid = ParallelDescriptor::second();
    BL_PROFILE_VAR("FLUID SOLVE",fluidSolve);
    if (FLUID::solve)
    {
       EvolveFluid(nstep,dt,prev_dt,time,stop_time,drag_timing);
       prev_dt = dt;
    }
    BL_PROFILE_VAR_STOP(fluidSolve);

    Real end_fluid = ParallelDescriptor::second() - start_fluid - drag_timing;
    ParallelDescriptor::ReduceRealMax(end_fluid, ParallelDescriptor::IOProcessorNumber());

#if 0
    const int nchem_species = FLUID::nchem_species;
    for (int lev = 0; lev <= finest_level; lev++)
    {
      auto& ld = *m_leveldata[lev];

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(*ld.vf_n,TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {
        Box const& bx = mfi.tilebox();

        Array4<Real const> const& vf_n     = ld.vf_n->const_array(mfi);

        ParallelFor(bx, nchem_species, [=]
            AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
              if (vf_n(i,j,k) != 1.0)
              {
                IntVect cell(i,j,k);
                print_state(*ld.X_k, cell, -1);
              }
            });
      } // mfi
    } // lev
#endif


    /****************************************************************************
     *                                                                          *
     * Evolve Particles (Using Particle MD)                                     *
     *                                                                          *
     ***************************************************************************/

    Real start_particles = ParallelDescriptor::second();

    BL_PROFILE_VAR("PARTICLES SOLVE", particlesSolve);

    int nsubsteps;

    if (DEM::solve)
    {
        pc->EvolveParticles(dt, particle_cost, knapsack_weight_type, nsubsteps);
        pc->split_particles();
    }

    BL_PROFILE_VAR_STOP(particlesSolve);

    Real end_particles = ParallelDescriptor::second() - start_particles;
    ParallelDescriptor::ReduceRealMax(end_particles, ParallelDescriptor::IOProcessorNumber());

    ComputeAndPrintSums();

    if (ParallelDescriptor::IOProcessor()) {
      if(FLUID::solve)
        std::cout << "   Time per fluid step      " << end_fluid << std::endl;

      if(DEM::solve)
        std::cout << "   Time per " << nsubsteps
                  << " particle steps " << end_particles << std::endl;

      if((DEM::solve) and FLUID::solve)
        std::cout << "   Coupling time per step   " << coupling_timing << std::endl;
    }

    BL_PROFILE_REGION_STOP("bmx::Evolve");
}
