#include <bmx_calc_fluid_coeffs.H>
#include <bmx_fluid_parms.H>

void calc_D_gk (const Box& bx,
                FArrayBox& D_gk_fab)
{
  const int nspecies_g = FLUID::nspecies;
  Gpu::DeviceVector< Real> D_gk0_d(nspecies_g);
  Gpu::copyAsync(Gpu::hostToDevice, FLUID::D_gk0.begin(), FLUID::D_gk0.end(), D_gk0_d.begin());

  Real* p_D_gk0 = D_gk0_d.data();

  Array4<Real> const& D_gk = D_gk_fab.array();

  amrex::Print() << " IN CALC_D_gk " << FLUID::D_gk0[0] << " " << FLUID::D_gk0[1] << std::endl;

  amrex::ParallelFor(bx, nspecies_g, [D_gk,p_D_gk0]
    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    { D_gk(i,j,k,n) = p_D_gk0[n]; });

  Gpu::synchronize();
}
