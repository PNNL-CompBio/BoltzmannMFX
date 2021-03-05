#include <AMReX_MultiFabUtil.H>
#include <bmx_diffusion_op.H>
#include <bmx_fluid_parms.H>

using namespace amrex;

//
// Implicit solve for chem_species mass fraction
//
void DiffusionOp::diffuse_chem_species (      Vector< MultiFab* >    X_gk_in,
                                        const Vector< MultiFab* >    D_gk_in,
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
    //      b: D_gk

    if(verbose > 0)
      amrex::Print() << "Diffusing chem_species mass fractions ..." << std::endl;

    // Set alpha and beta
    chem_species_matrix->setScalars(1.0, theta*dt);

    // Number of fluid chem_species
    const int nchem_species_g = FLUID::nchem_species;

    Vector<BCRec> bcs_X; 
    bcs_X.resize(3*nchem_species_g);

    for(int lev = 0; lev <= finest_level; lev++)
    {
        average_cellcenter_to_face( GetArrOfPtrs(chem_species_b[lev]), *D_gk_in[lev], geom[lev], nchem_species_g );

        chem_species_b[lev][0]->FillBoundary(geom[lev].periodicity());
        chem_species_b[lev][1]->FillBoundary(geom[lev].periodicity());
        chem_species_b[lev][2]->FillBoundary(geom[lev].periodicity());

        // Zero out the coefficients in the high-z part
        const Box& domain = geom[lev].Domain();
        const int zhi = domain.bigEnd()[2] / 2;

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(*X_gk_in[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
          Box const& bx = mfi.tilebox();
          if (bx.ok())
          {
            Array4<Real> const& bx_arr = chem_species_b[lev][0]->array(mfi);
            Array4<Real> const& by_arr = chem_species_b[lev][1]->array(mfi);
            Array4<Real> const& bz_arr = chem_species_b[lev][2]->array(mfi);

            Box const& xbx = mfi.growntilebox(IntVect(1,0,0));
            amrex::ParallelFor(xbx, nchem_species_g, [=]
              AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
                if (k > zhi/2)
                    bx_arr(i,j,k,n) = 0.;
            });

            Box const& ybx = mfi.growntilebox(IntVect(0,1,0));
            amrex::ParallelFor(ybx, nchem_species_g, [=]
              AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
                if (k > zhi/2)
                    by_arr(i,j,k,n) = 0.;
            });

            Box const& zbx = mfi.growntilebox(IntVect(0,0,1));
            amrex::ParallelFor(zbx, nchem_species_g, [=]
              AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
                if (k > zhi/2)
                    bz_arr(i,j,k,n) = 0.;
            });
          }
        }

        // This sets the coefficients
        chem_species_matrix->setACoeffs (lev, 1.);
        chem_species_matrix->setBCoeffs (lev, GetArrOfConstPtrs(chem_species_b[lev]));

        // Zero these out just to have a clean start because they have 3 components
        //      (due to re-use with velocity solve)
        chem_species_phi[lev]->setVal(0.0);
        chem_species_rhs[lev]->setVal(0.0);

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
            Array4<Real      > const& rhs_arr     = chem_species_rhs[lev]->array(mfi);

            amrex::ParallelFor(bx, [X_gk_arr,rhs_arr,nchem_species_g]
              AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
              for (int n(0); n < nchem_species_g; ++n)
                rhs_arr(i,j,k,n) = X_gk_arr(i,j,k,n);
            });
          }
        }

        MultiFab::Copy(*chem_species_phi[lev], *X_gk_in[lev], 0, 0, nchem_species_g, 1);
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
        MultiFab::Copy(*X_gk_in[lev], *chem_species_phi[lev], 0, 0, nchem_species_g, 1);
    }
}
