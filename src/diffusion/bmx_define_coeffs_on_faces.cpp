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
void DiffusionOp::define_coeffs_on_faces ( const Vector< MultiFab const* > D_k_in,
                                           const Vector< MultiFab const* > vf_in )
{
    BL_PROFILE("DiffusionOp::define_coeffs_on_faces");

    int finest_level = amrcore->finestLevel();

    // Number of fluid chem_species
    const int nchem_species = FLUID::nchem_species;

    // Make a copy of diffusion coefficients
    Vector<MultiFab> D_k_adj(D_k_in.size());

    for (int lev = 0; lev <= finest_level; lev++) 
    {
      const BoxArray &ba = D_k_in[lev]->boxArray();
      const DistributionMapping &dm = D_k_in[lev]->DistributionMap();
      int ncomp = D_k_in[lev]->nComp();
      int ngrow = D_k_in[lev]->nGrow();
      D_k_adj[lev].define(ba, dm, ncomp, ngrow);
    //  D_k_adj[lev] = *D_k_in[lev];

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(*D_k_in[lev],TilingIfNotGPU());
          mfi.isValid(); ++mfi)
      {
        Array4<Real const> const&   d0_arr = D_k_in[lev]->array(mfi);
        Array4<Real      > const& dadj_arr = D_k_adj[lev].array(mfi);
        Array4<Real const> const&   vf_arr = vf_in[lev]->array(mfi);

        Box const& bx = mfi.growntilebox(1);

        // Evaluate diffusion coefficient in grid cell modified by value of
        // volume fraction.
        // These formulas assume that volume fraction = fraction of fluid in grid cell.
        amrex::ParallelFor(bx, nchem_species, [=]
          AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
          {
            if (vf_arr(i,j,k) > 0.0) {
              dadj_arr(i,j,k,n) = vf_arr(i,j,k) * d0_arr(i,j,k,n) / (1.0-0.45*log(vf_arr(i,j,k)));
            } else {
              dadj_arr(i,j,k,n) = 0.1*d0_arr(i,j,k,n);
            }
            if (dadj_arr(i,j,k,n) < 0.1*d0_arr(i,j,k,n)) dadj_arr(i,j,k,n) = 0.1*d0_arr(i,j,k,n);
          });

      } // MFIter
    } // lev

    for (int lev = 0; lev <= finest_level; lev++)
    {
        bool use_harmonic_averaging = true;
        average_cellcenter_to_face( GetArrOfPtrs(chem_species_b[lev]), D_k_adj[lev], 
                                                 amrcore->Geom(lev), nchem_species,
                                                 use_harmonic_averaging );

        auto& geom = amrcore->Geom(lev);

        Real dz = geom.CellSize(2);
        Real zhi = FLUID::surface_location;

        //
        // The code below -- currently commented out --  is if we want to over-ride 
        // the harmonic average with some problem-specific information
        //
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(D_k_adj[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
          Array4<Real> const& bx_arr = chem_species_b[lev][0]->array(mfi);
          Array4<Real> const& by_arr = chem_species_b[lev][1]->array(mfi);
          Array4<Real> const& bz_arr = chem_species_b[lev][2]->array(mfi);


          Array4<Real const> const& vf_arr = vf_in[lev]->array(mfi);

          Box const& xbx = mfi.nodaltilebox(0);

          amrex::ParallelFor(xbx, nchem_species, [=]
            AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
          {
              Real z = (static_cast<double>(k)+.5) * dz; 

              if (z > zhi)
                 if (vf_arr(i,j,k) == 1. || vf_arr(i-1,j,k) == 1.)
                 {
               //  printf("bx[%d,%d,%d,%d] is zero\n",i,j,k,n);
                     bx_arr(i,j,k,n) = 0.;
                     }
          });

          Box const& ybx = mfi.nodaltilebox(1);

          amrex::ParallelFor(ybx, nchem_species, [=]
            AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
          {
              Real z = (static_cast<double>(k)+.5) * dz; 

              if (z > zhi)
                 if (vf_arr(i,j,k) == 1. || vf_arr(i,j-1,k) == 1.)
                 {
              //   printf("by[%d,%d,%d,%d] is zero\n",i,j,k,n);
                     by_arr(i,j,k,n) = 0.;
                     }
          });

          Box const& zbx = mfi.nodaltilebox(2);

          amrex::ParallelFor(zbx, nchem_species, [=]
            AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
          {
              Real zk   = (static_cast<double>(k)+.5) * dz; 
              Real zkm1 = (static_cast<double>(k)-.5) * dz; 

              if ( (zk   > zhi && vf_arr(i,j,k  ) == 1.) || // zk is above surface and is empty

                   (zkm1 > zhi && vf_arr(i,j,k-1) == 1.) )  // zk is above surface and might have particle
                                                            //   but cell below is empty
                 {
               //  printf("bz[%d,%d,%d,%d] is zero\n",i,j,k,n);
                       bz_arr(i,j,k,n) = 0.;
                       }
          });
        } // mfi
    } // lev
}
