#include <AMReX_MultiFabUtil.H>
#include <bmx_diffusion_op.H>
#include <bmx_fluid_parms.H>

using namespace amrex;

//
// Implicit solve for species mass fraction
//
void DiffusionOp::diffuse_species (      Vector< MultiFab* >    X_gk_in,
                                   const Vector< MultiFab* >    D_gk_in,
                                   Real theta, Real dt)
{
    BL_PROFILE("DiffusionOp::diffuse_species");

    int finest_level = amrcore->finestLevel();

    // Update the coefficients of the matrix going into the solve based on the current state of the
    // simulation. Recall that the relevant matrix is
    //
    //      alpha a - beta div ( b grad )   <--->   rho - dt div ( mu_s grad )
    //
    // So the constants and variable coefficients are:
    //
    //      alpha: 1
    //      beta: theta * dt
    //      a: 1
    //      b: D_gk

    if(verbose > 0)
      amrex::Print() << "Diffusing species mass fractions ..." << std::endl;

    // Set alpha and beta
    species_matrix->setScalars(1.0, theta*dt);

    // Number of fluid species
    const int nspecies_g = FLUID::nspecies;

    Vector<BCRec> bcs_X; 
    bcs_X.resize(3*nspecies_g);

    for(int lev = 0; lev <= finest_level; lev++)
    {
        average_cellcenter_to_face( GetArrOfPtrs(species_b[lev]), *D_gk_in[lev], geom[lev], nspecies_g );

        // This sets the coefficients
        species_matrix->setACoeffs (lev, 1.);
        species_matrix->setBCoeffs (lev, GetArrOfConstPtrs(species_b[lev]));

        // Zero these out just to have a clean start because they have 3 components
        //      (due to re-use with velocity solve)
        species_phi[lev]->setVal(0.0);
        species_rhs[lev]->setVal(0.0);

        // Set rhs equal to X_gk and
        // Multiply rhs by (D_gk) -- we are solving
        //
        //      X_star = X_old + theta * dt (div (D_gk grad X_gk) )
        //
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(*X_gk_in[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
          Box const& bx = mfi.growntilebox(IntVect(0));

          if (bx.ok())
          {
            Array4<Real const> const& X_gk_arr    = X_gk_in[lev]->const_array(mfi);
            Array4<Real      > const& rhs_arr     = species_rhs[lev]->array(mfi);

            amrex::ParallelFor(bx, [X_gk_arr,rhs_arr,nspecies_g]
              AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
              for (int n(0); n < nspecies_g; ++n)
                rhs_arr(i,j,k,n) = X_gk_arr(i,j,k,n);
            });
          }
        }

        MultiFab::Copy(*species_phi[lev], *X_gk_in[lev], 0, 0, nspecies_g, 1);
        species_matrix->setLevelBC(lev, GetVecOfConstPtrs(species_phi)[lev]);
    }

    MLMG solver(*species_matrix);
    setSolverSettings(solver);

    // This ensures that ghost cells of sol are correctly filled when returned from the solver
    solver.setFinalFillBC(true);

    solver.solve(GetVecOfPtrs(species_phi), GetVecOfConstPtrs(species_rhs), mg_rtol, mg_atol);

    for(int lev = 0; lev <= finest_level; lev++)
    {
        species_phi[lev]->FillBoundary(geom[lev].periodicity());
        MultiFab::Copy(*X_gk_in[lev], *species_phi[lev], 0, 0, nspecies_g, 1);
    }
}
