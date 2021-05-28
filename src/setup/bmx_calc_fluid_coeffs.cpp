#include <bmx_calc_fluid_coeffs.H>
#include <bmx_fluid_parms.H>

void calc_D_k (const Box& bx,
               const Box& domain,
               const Real dx,
               const Real dy,
               const Real dz,
               FArrayBox& D_k_fab)
{
  const int nchem_species = FLUID::nchem_species;
  Gpu::DeviceVector< Real> D_k0_d(nchem_species);
  Gpu::copyAsync(Gpu::hostToDevice, FLUID::D_k0.begin(), FLUID::D_k0.end(), D_k0_d.begin());

  Real* p_D_k0 = D_k0_d.data();

  Array4<Real> const& D_k = D_k_fab.array();

  amrex::Print() << " IN CALC_D_k " << FLUID::D_k0[0] << " " << FLUID::D_k0[1] << std::endl;

  ParmParse pp("geometry");
  amrex::Vector<Real> bmin, bmax;
  pp.getarr("prob_lo",bmin);
  pp.getarr("prob_hi",bmax);

#if 1
  // We set the coeffs at cell centers to D_k everywhere. Diffusion is
  // controlled when we calculate diffusion coefficients on cell faces
  amrex::ParallelFor(bx, nchem_species, [D_k,p_D_k0]
    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
      D_k(i,j,k,n) = p_D_k0[n];
    });

#else
  // We set the coeffs at cell centers to D_k in the lower region and 0 above zhi
  Real zhi = FLUID::surface_location;
  const int ihi = int(zhi/dz+0.00000000001);

  amrex::ParallelFor(bx, nchem_species, [D_k,p_D_k0,ihi]
    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    { if (k < ihi)
         D_k(i,j,k,n) = p_D_k0[n];
      else
         D_k(i,j,k,n) = 0.0; 
    });
#endif

  Gpu::synchronize();
}
