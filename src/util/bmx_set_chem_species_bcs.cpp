#include <bmx.H>

#include <bmx_fluid_parms.H>
#include <bmx_chem_species_parms.H>

using namespace amrex;

//
// Set the BCs 
//
void
bmx::bmx_set_chem_species_bcs (Real time,
                            Vector< MultiFab* > const& X_gk_in,
                            Vector< MultiFab* > const& D_gk_in)
{
  BL_PROFILE("bmx::bmx_set_chem_species_bcs()");

  for (int lev = 0; lev < nlev; lev++)
  {
     Box domain(geom[lev].Domain());

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
     for (MFIter mfi(*(m_leveldata[lev]->X_gk), false); mfi.isValid(); ++mfi)
     {
        set_chem_species_diffusivities_g_bcs(time, lev, (*D_gk_in[lev])[mfi], domain);
     }

     X_gk_in[lev]->FillBoundary(geom[lev].periodicity());
     D_gk_in[lev]->FillBoundary(geom[lev].periodicity());
  }
}

void 
bmx::set_chem_species_diffusivities_g_bcs (Real time,
                                           const int lev,
                                           FArrayBox& scal_fab,
                                           const Box& domain)
{
  IntVect dom_lo(domain.loVect());
  IntVect dom_hi(domain.hiVect());

  Array4<const int> const& bct_ilo = bc_ilo[lev]->array();
  Array4<const int> const& bct_ihi = bc_ihi[lev]->array();
  Array4<const int> const& bct_jlo = bc_jlo[lev]->array();
  Array4<const int> const& bct_jhi = bc_jhi[lev]->array();
  Array4<const int> const& bct_klo = bc_klo[lev]->array();
  Array4<const int> const& bct_khi = bc_khi[lev]->array();

  Array4<Real> const& scal_arr = scal_fab.array();

  const int nchem_species_g = FLUID::nchem_species;

  Gpu::DeviceVector< Real > D_gk0_d(nchem_species_g);
  Gpu::copyAsync(Gpu::hostToDevice, FLUID::D_gk0.begin(), FLUID::D_gk0.end(), D_gk0_d.begin());

  Real* p_D_gk0 = D_gk0_d.data();

  IntVect scal_lo(scal_fab.loVect());
  IntVect scal_hi(scal_fab.hiVect());

  const int nlft = std::max(0, dom_lo[0]-scal_lo[0]);
  const int nbot = std::max(0, dom_lo[1]-scal_lo[1]);
  const int ndwn = std::max(0, dom_lo[2]-scal_lo[2]);

  const int nrgt = std::max(0, scal_hi[0]-dom_hi[0]);
  const int ntop = std::max(0, scal_hi[1]-dom_hi[1]);
  const int nup  = std::max(0, scal_hi[2]-dom_hi[2]);

  const int minf = bc_list.get_minf();
  const int pout = bc_list.get_pout();

  if (nlft > 0)
  {
    IntVect bx_yz_lo_hi_3D(scal_hi);
    bx_yz_lo_hi_3D[0] = dom_lo[0]-1;
    const Box bx_yz_lo_3D(scal_lo, bx_yz_lo_hi_3D);

    int ilo = dom_lo[0];

    amrex::ParallelFor(bx_yz_lo_3D, nchem_species_g, [bct_ilo,ilo,scal_arr]
      AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
      const int bct = bct_ilo(ilo-1,j,k,0);

      scal_arr(i,j,k,n) = scal_arr(ilo,j,k,n);

    });
  }

  if (nrgt > 0)
  {
    IntVect bx_yz_hi_lo_3D(scal_lo);
    bx_yz_hi_lo_3D[0] = dom_hi[0]+1;
    const Box bx_yz_hi_3D(bx_yz_hi_lo_3D, scal_hi);

    int ihi = dom_hi[0];

    amrex::ParallelFor(bx_yz_hi_3D, nchem_species_g, [bct_ihi,ihi,scal_arr]
      AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
      const int bct = bct_ihi(ihi+1,j,k,0);

      scal_arr(i,j,k,n) = scal_arr(ihi,j,k,n);

    });
  }

  if (nbot > 0)
  {
    IntVect bx_xz_lo_hi_3D(scal_hi);
    bx_xz_lo_hi_3D[1] = dom_lo[1]-1;
    const Box bx_xz_lo_3D(scal_lo, bx_xz_lo_hi_3D);

    int jlo = dom_lo[1];

    amrex::ParallelFor(bx_xz_lo_3D, nchem_species_g, [bct_jlo,jlo,scal_arr]
      AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
      const int bct = bct_jlo(i,jlo-1,k,0);

      scal_arr(i,j,k,n) = scal_arr(i,jlo,k,n);
    });
  }

  if (ntop > 0)
  {
    IntVect bx_xz_hi_lo_3D(scal_lo);
    bx_xz_hi_lo_3D[1] = dom_hi[1]+1;
    const Box bx_xz_hi_3D(bx_xz_hi_lo_3D, scal_hi);

    int jhi = dom_hi[1];

    amrex::ParallelFor(bx_xz_hi_3D, nchem_species_g, [bct_jhi,jhi,scal_arr]
      AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
      const int bct = bct_jhi(i,jhi+1,k,0);

      scal_arr(i,j,k,n) = scal_arr(i,jhi,k,n);
    });
  }

  if (ndwn > 0)
  {
    IntVect bx_xy_lo_hi_3D(scal_hi);
    bx_xy_lo_hi_3D[2] = dom_lo[2]-1;
    const Box bx_xy_lo_3D(scal_lo, bx_xy_lo_hi_3D);

    int klo = dom_lo[2];

    amrex::ParallelFor(bx_xy_lo_3D, nchem_species_g, [bct_klo,klo,scal_arr]
      AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
      const int bct = bct_klo(i,j,klo-1,0);

      scal_arr(i,j,k,n) = scal_arr(i,j,klo,n);
    });
  }

  if (nup > 0)
  {
    IntVect bx_xy_hi_lo_3D(scal_lo);
    bx_xy_hi_lo_3D[2] = dom_hi[2]+1;
    const Box bx_xy_hi_3D(bx_xy_hi_lo_3D, scal_hi);

    int khi = dom_hi[2];

    amrex::ParallelFor(bx_xy_hi_3D, nchem_species_g, [bct_khi,khi,scal_arr]
      AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
      const int bct = bct_khi(i,j,khi+1,0);

      scal_arr(i,j,k,n) = scal_arr(i,j,khi,n);
    });
  }

  Gpu::synchronize();
}

