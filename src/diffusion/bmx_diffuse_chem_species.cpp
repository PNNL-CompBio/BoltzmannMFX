#include <AMReX_MultiFabUtil.H>
#include <bmx_diffusion_op.H>
#include <bmx_fluid_parms.H>

using namespace amrex;

//
// Implicit solve for chem_species mass fraction
//
void DiffusionOp::diffuse_chem_species (      Vector< MultiFab*      > X_k_in,
                                        const Vector< MultiFab const*> D_k_in,
                                        Real theta, Real dt)
{
    BL_PROFILE("DiffusionOp::diffuse_chem_species");

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
    //      b: D_k

    if(verbose > 0)
      amrex::Print() << "Diffusing chem_species mass fractions ..." << std::endl;

    // Set alpha and beta
    chem_species_matrix->setScalars(1.0, theta*dt);

    // Number of fluid chem_species
    const int nchem_species = FLUID::nchem_species;

    define_coeffs_on_faces(D_k_in);

    for(int lev = 0; lev <= finest_level; lev++)
    {
        // This sets the coefficients
        chem_species_matrix->setACoeffs (lev, 1.);
        chem_species_matrix->setBCoeffs (lev, GetArrOfConstPtrs(chem_species_b[lev]));

        // Zero these out just to have a clean start because they have 3 components
        //      (due to re-use with velocity solve)
        chem_species_phi[lev]->setVal(0.0);
        chem_species_rhs[lev]->setVal(0.0);

        // Set rhs equal to X_k and
        // Multiply rhs by (D_k) -- we are solving
        //
        //      X_star = X_old + theta * dt (div (D_k grad X_k) )
        //
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(*X_k_in[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
          Box const& bx = mfi.growntilebox(IntVect(0));

          if (bx.ok())
          {
            Array4<Real const> const& X_k_arr    = X_k_in[lev]->const_array(mfi);
            Array4<Real      > const& rhs_arr     = chem_species_rhs[lev]->array(mfi);

            amrex::ParallelFor(bx, [X_k_arr,rhs_arr,nchem_species]
              AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
              for (int n(0); n < nchem_species; ++n)
                rhs_arr(i,j,k,n) = X_k_arr(i,j,k,n);
            });
          }
        }

        MultiFab::Copy(*chem_species_phi[lev], *X_k_in[lev], 0, 0, nchem_species, 1);
        chem_species_matrix->setLevelBC(lev, GetVecOfConstPtrs(chem_species_phi)[lev]);
    }

    MLMG solver(*chem_species_matrix);
    setSolverSettings(solver);

    // This ensures that ghost cells of sol are correctly filled when returned from the solver
    solver.setFinalFillBC(true);

    solver.solve(GetVecOfPtrs(chem_species_phi), GetVecOfConstPtrs(chem_species_rhs), mg_rtol, mg_atol);

    for(int lev = 0; lev <= finest_level; lev++)
    {
        chem_species_phi[lev]->FillBoundary(geom[lev].periodicity());
        MultiFab::Copy(*X_k_in[lev], *chem_species_phi[lev], 0, 0, nchem_species, 1);
    }
}
