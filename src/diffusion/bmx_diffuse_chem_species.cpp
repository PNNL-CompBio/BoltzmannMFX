//
//     Copyright (c) 2013 Battelle Memorial Institute
//     Licensed under modified BSD License. A copy of this license can be found
//     in the LICENSE file in the top level directory of this distribution.
//
#include <AMReX_MultiFabUtil.H>
#include <bmx_diffusion_op.H>
#include <bmx_fluid_parms.H>

using namespace amrex;

//
// Implicit solve for chem_species mass fraction
//
void DiffusionOp::diffuse_chem_species (      Vector< MultiFab*      > X_k_in,
                                        const Vector< MultiFab const*> D_k_in,
                                        const Vector< MultiFab const*>  vf_in,
                                        Real omt, Real dt)
{
    BL_PROFILE("DiffusionOp::diffuse_chem_species");

    int finest_level = amrcore->finestLevel();

    // Update the coefficients of the matrix going into the solve based on the current state of the
    // simulation. Recall that the relevant matrix is
    //
    //      alpha a - beta div ( b grad )   <--->   vf - (1-theta) dt div ( D_k grad )
    //
    // So the constants and variable coefficients are:
    //
    //      alpha: 1
    //      beta: omt * dt
    //      a: vf
    //      b: D_k

    if(verbose > 0)
      amrex::Print() << "Diffusing chem_species mass fractions ..." << std::endl;

    // Set alpha and beta
    chem_species_matrix->setScalars(1.0, omt*dt);

    // Number of fluid chem_species
    const int nchem_species = FLUID::nchem_species;

    define_coeffs_on_faces(D_k_in, vf_in);

    for(int lev = 0; lev <= finest_level; lev++)
    {
        // This sets the coefficients
        chem_species_matrix->setACoeffs (lev, (*vf_in[lev]));
        chem_species_matrix->setBCoeffs (lev, GetArrOfConstPtrs(chem_species_b[lev]));

        // Zero these out just to have a clean start because they have 3 components
        //      (due to re-use with velocity solve)
        chem_species_phi[lev]->setVal(0.0);
        chem_species_rhs[lev]->setVal(0.0);

        // Set rhs equal to X_k and
        // Multiply rhs by (D_k) -- we are solving
        //
        //      vf X_star = vf X_old + (1-theta) * dt (div (D_k grad X_k) )
        //
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(*X_k_in[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
          Box const& bx = mfi.growntilebox(IntVect(0));

          if (bx.ok())
          {
            Array4<Real const> const& X_k_arr     = X_k_in[lev]->const_array(mfi);
            Array4<Real      > const& rhs_arr     = chem_species_rhs[lev]->array(mfi);
            Array4<Real const> const&  vf_arr     = vf_in[lev]->const_array(mfi);

            amrex::ParallelFor(bx, [X_k_arr,rhs_arr,vf_arr,nchem_species]
              AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
              for (int n(0); n < nchem_species; ++n)
              {
                rhs_arr(i,j,k,n) = X_k_arr(i,j,k,n) * vf_arr(i,j,k);
              }
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
        chem_species_phi[lev]->FillBoundary(amrcore->Geom(lev).periodicity());
        MultiFab::Copy(*X_k_in[lev], *chem_species_phi[lev], 0, 0, nchem_species, 1);
    }
}
