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

    for (int lev = 0; lev <= finest_level; lev++)
    {
        bool use_harmonic_averaging = true;
        average_cellcenter_to_face( GetArrOfPtrs(chem_species_b[lev]), *D_k_in[lev], 
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
        for (MFIter mfi(*D_k_in[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
          Array4<Real> const& bx_arr = chem_species_b[lev][0]->array(mfi);
          Array4<Real> const& by_arr = chem_species_b[lev][1]->array(mfi);
          Array4<Real> const& bz_arr = chem_species_b[lev][2]->array(mfi);

          Array4<Real const> const& vf_arr = vf_in[lev]->array(mfi);

          Box const& xbx = mfi.growntilebox(IntVect(1,0,0));
          amrex::ParallelFor(xbx, nchem_species, [=]
            AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
          {
              Real z = (k+.5) * dz; 

              if (z > zhi)
                 if (vf_arr(i,j,k) == 1. || vf_arr(i-1,j,k) == 1.)
                     bx_arr(i,j,k,n) = 0.;
          });

          Box const& ybx = mfi.growntilebox(IntVect(0,1,0));
          amrex::ParallelFor(ybx, nchem_species, [=]
            AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
          {
              Real z = (k+.5) * dz; 

              if (z > zhi)
                 if (vf_arr(i,j,k) == 1. || vf_arr(i,j-1,k) == 1.)
                     by_arr(i,j,k,n) = 0.;
          });

          Box const& zbx = mfi.growntilebox(IntVect(0,0,1));
          amrex::ParallelFor(zbx, nchem_species, [=]
            AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
          {
              Real zk   = (k+.5) * dz; 
              Real zkm1 = (k-.5) * dz; 

              if ( (zk   > zhi && vf_arr(i,j,k  ) == 1.) || 
                   (zkm1 > zhi && vf_arr(i,j,k-1) == 1.) )
                       bz_arr(i,j,k,n) = 0.;
          });
        } // mfi
    } // lev
}
