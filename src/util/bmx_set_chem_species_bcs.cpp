//
//     Copyright (c) 2013 Battelle Memorial Institute
//     Licensed under modified BSD License. A copy of this license can be found
//     in the LICENSE file in the top level directory of this distribution.
//
#include <bmx.H>

#include <bmx_fluid_parms.H>
#include <bmx_chem_species_parms.H>

using namespace amrex;

//
// Set the BCs 
//
void
bmx::bmx_set_chem_species_bcs (Real time,
                               Vector< MultiFab* > const& X_k_in,
                               Vector< MultiFab* > const& D_k_in)
{
  BL_PROFILE("bmx::bmx_set_chem_species_bcs()");

  for (int lev = 0; lev <= finest_level; lev++)
  {
     Box domain(geom[lev].Domain());

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
     for (MFIter mfi(*(m_leveldata[lev]->X_k), false); mfi.isValid(); ++mfi)
     {
        set_neumann_bcs(time, lev, (*X_k_in[lev])[mfi], geom[lev].data());
        set_neumann_bcs(time, lev, (*D_k_in[lev])[mfi], geom[lev].data());
     }

     X_k_in[lev]->FillBoundary(geom[lev].periodicity());
     D_k_in[lev]->FillBoundary(geom[lev].periodicity());
  }
}

void 
bmx::set_neumann_bcs (Real /*time*/,
                      const int /*lev*/,
                      FArrayBox& scal_fab,
                      const GeometryData& geom_data)
{
  const Box& domain = geom_data.Domain();

  IntVect dom_lo(domain.loVect());
  IntVect dom_hi(domain.hiVect());

  Array4<Real> const& scal_arr = scal_fab.array();

  const int nchem_species = FLUID::nchem_species;

  Gpu::DeviceVector< Real > D_k0_d(nchem_species);
  Gpu::copyAsync(Gpu::hostToDevice, FLUID::D_k0.begin(), FLUID::D_k0.end(), D_k0_d.begin());

  IntVect scal_lo(scal_fab.loVect());
  IntVect scal_hi(scal_fab.hiVect());

  const int nlft = geom_data.isPeriodic(0) ? 0 : std::max(0, dom_lo[0]-scal_lo[0]);
  const int nbot = geom_data.isPeriodic(1) ? 0 : std::max(0, dom_lo[1]-scal_lo[1]);
  const int ndwn = geom_data.isPeriodic(2) ? 0 : std::max(0, dom_lo[2]-scal_lo[2]);

  const int nrgt = geom_data.isPeriodic(0) ? 0 : std::max(0, scal_hi[0]-dom_hi[0]);
  const int ntop = geom_data.isPeriodic(1) ? 0 : std::max(0, scal_hi[1]-dom_hi[1]);
  const int nup  = geom_data.isPeriodic(2) ? 0 : std::max(0, scal_hi[2]-dom_hi[2]);

  if (nlft > 0)
  {
    IntVect bx_yz_lo_hi_3D(scal_hi);
    bx_yz_lo_hi_3D[0] = dom_lo[0]-1;
    const Box bx_yz_lo_3D(scal_lo, bx_yz_lo_hi_3D);

    int ilo = dom_lo[0];

    amrex::ParallelFor(bx_yz_lo_3D, nchem_species, [ilo,scal_arr]
      AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
      scal_arr(i,j,k,n) = scal_arr(ilo,j,k,n);

    });
  }

  if (nrgt > 0)
  {
    IntVect bx_yz_hi_lo_3D(scal_lo);
    bx_yz_hi_lo_3D[0] = dom_hi[0]+1;
    const Box bx_yz_hi_3D(bx_yz_hi_lo_3D, scal_hi);

    int ihi = dom_hi[0];

    amrex::ParallelFor(bx_yz_hi_3D, nchem_species, [ihi,scal_arr]
      AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
      scal_arr(i,j,k,n) = scal_arr(ihi,j,k,n);

    });
  }

  if (nbot > 0)
  {
    IntVect bx_xz_lo_hi_3D(scal_hi);
    bx_xz_lo_hi_3D[1] = dom_lo[1]-1;
    const Box bx_xz_lo_3D(scal_lo, bx_xz_lo_hi_3D);

    int jlo = dom_lo[1];

    amrex::ParallelFor(bx_xz_lo_3D, nchem_species, [jlo,scal_arr]
      AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
      scal_arr(i,j,k,n) = scal_arr(i,jlo,k,n);
    });
  }

  if (ntop > 0)
  {
    IntVect bx_xz_hi_lo_3D(scal_lo);
    bx_xz_hi_lo_3D[1] = dom_hi[1]+1;
    const Box bx_xz_hi_3D(bx_xz_hi_lo_3D, scal_hi);

    int jhi = dom_hi[1];

    amrex::ParallelFor(bx_xz_hi_3D, nchem_species, [jhi,scal_arr]
      AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
      scal_arr(i,j,k,n) = scal_arr(i,jhi,k,n);
    });
  }

  if (ndwn > 0)
  {
    IntVect bx_xy_lo_hi_3D(scal_hi);
    bx_xy_lo_hi_3D[2] = dom_lo[2]-1;
    const Box bx_xy_lo_3D(scal_lo, bx_xy_lo_hi_3D);

    int klo = dom_lo[2];

    amrex::ParallelFor(bx_xy_lo_3D, nchem_species, [klo,scal_arr]
      AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
      scal_arr(i,j,k,n) = scal_arr(i,j,klo,n);
    });
  }

  if (nup > 0)
  {
    IntVect bx_xy_hi_lo_3D(scal_lo);
    bx_xy_hi_lo_3D[2] = dom_hi[2]+1;
    const Box bx_xy_hi_3D(bx_xy_hi_lo_3D, scal_hi);

    int khi = dom_hi[2];

    amrex::ParallelFor(bx_xy_hi_3D, nchem_species, [khi,scal_arr]
      AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
      scal_arr(i,j,k,n) = scal_arr(i,j,khi,n);
    });
  }

  Gpu::synchronize();
}

