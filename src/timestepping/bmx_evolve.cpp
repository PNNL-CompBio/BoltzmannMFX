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
     * Evolve Fluid                                                             *
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

    int lev = 0;
    const auto p_lo = Geom(lev).ProbLoArray();
    const auto p_hi = Geom(lev).ProbHiArray();

    Real dx = geom[lev].CellSize(0);
    Real dy = geom[lev].CellSize(1);
    Real dz = geom[lev].CellSize(2);

    Real domain_vol = (p_hi[2]-p_lo[2])*(p_hi[1]-p_lo[1])*(p_hi[0]-p_lo[0]);

    Real fluid_vol = volSum(0,*(m_leveldata[lev]->vf_n),false) * dx * dy * dz;

    Real particle_vol = pc->computeParticleVolume();

    Real A_in_fluid     = volWgtSum(0,*(m_leveldata[lev]->vf_n), *(m_leveldata[lev]->X_k), 0, false) * dx * dy * dz;
    Real A_in_particles = pc->computeParticleContent(22);

    Real B_in_fluid     = volWgtSum(0,*(m_leveldata[lev]->vf_n), *(m_leveldata[lev]->X_k), 1, false) * dx * dy * dz;
    Real B_in_particles = pc->computeParticleContent(23);

    Real C_in_fluid     = volWgtSum(0,*(m_leveldata[lev]->vf_n), *(m_leveldata[lev]->X_k), 2, false) * dx * dy * dz;
    Real C_in_particles = pc->computeParticleContent(24);

    amrex::Print() << "Domain   volume : " << domain_vol << std::endl;
    amrex::Print() << "Fluid    volume : " << fluid_vol  << std::endl;
    amrex::Print() << "Particle volume : " << particle_vol << std::endl;
    amrex::Print() << "Particle + Fluid: " << fluid_vol+particle_vol << std::endl;

    amrex::Print() << " A in fluid        : " << A_in_fluid << std::endl;
    amrex::Print() << " B in fluid        : " << B_in_fluid << std::endl;
    amrex::Print() << " C in fluid        : " << C_in_fluid << std::endl;
    amrex::Print() << " A in particles    : " << A_in_particles << std::endl;
    amrex::Print() << " B in particles    : " << B_in_particles << std::endl;
    amrex::Print() << " C in particles    : " << C_in_particles << std::endl;
    amrex::Print() << " A+C   in fluid    : " << A_in_fluid + C_in_fluid << std::endl;
    amrex::Print() << " A+C   in particles: " << A_in_particles + C_in_particles << std::endl;
    amrex::Print() << " Total A + C       : " << A_in_fluid + A_in_particles +
                                                 C_in_fluid + C_in_particles << std::endl;

    if (std::abs(domain_vol - (fluid_vol+particle_vol)) > 1.e-12 * domain_vol)
       amrex::Abort("Volumes don't match!");

    BL_PROFILE_VAR_STOP(particlesSolve);

    Real end_particles = ParallelDescriptor::second() - start_particles;
    ParallelDescriptor::ReduceRealMax(end_particles, ParallelDescriptor::IOProcessorNumber());

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
