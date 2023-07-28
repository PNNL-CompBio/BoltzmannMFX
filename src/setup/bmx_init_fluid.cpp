//
//     Copyright (c) 2013 Battelle Memorial Institute
//     Licensed under modified BSD License. A copy of this license can be found
//     in the LICENSE file in the top level directory of this distribution.
//
#include <bmx_init_fluid.H>

#include <bmx_calc_fluid_coeffs.H>
#include <bmx_calc_cell.H>
#include <bmx_fluid_parms.H>
#include <bmx_chem_species_parms.H>

using namespace amrex;

//!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
//!                                                                      !
//!  Purpose: Set fluid chem_species mass fractions initial conditions.       !
//!                                                                      !
//!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
void set_ic_chem_species (const Box& sbx,
                          const Box& /*domain*/,
                          const Real   /*dx*/,
                          const Real   /*dy*/,
                          const Real   dz,
                          const GpuArray<Real, 3>& p_lo,
                          const GpuArray<Real, 3>& /*p_hi*/,
                          const Array4<Real> X_k_arr)
{
  Vector<Real> h_vec(NUM_CHEM_COMPONENTS,0.0);
  Gpu::DeviceVector<Real> d_vec(NUM_CHEM_COMPONENTS);
  int ii;
  // Initialize h_vec (on host)
  for (ii=0; ii<NUM_CHEM_COMPONENTS; ii++) h_vec[ii] = FLUID::init_conc[ii];
  // Copy host to device
  Gpu::copy(Gpu::hostToDevice,h_vec.begin(),h_vec.end(),d_vec.begin());
  auto vec_ptr = d_vec.data();

  const int nchem_species = X_k_arr.nComp();

  // We set the coeffs at cell centers to D_k in the lower region and 0 above zhi
  Real zhi = FLUID::surface_location;

  ParallelFor(sbx, nchem_species, [=]
       AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
      { 
          // Real x   = (i+.5)*dx+p_lo[0];
          // Real y   = (j+.5)*dy+p_lo[1];
          Real z   = (k+.5)*dz+p_lo[2];

#if 0
          Real ra = 0.1;
          Real rb = 0.2;

          Real rsq = (x-.5)*(x-.5) +  (y-.5)*(y-.5) +  (z-.5)*(z-.5); 
          // non-zero only within fixed radius
          if (n == 0)
          {
              if (rsq < ra*ra) 
                  X_k_arr(i,j,k,n) = 1.0;
              else
                  X_k_arr(i,j,k,n) = 0.0;
          } else {
              if (rsq < rb*rb) 
                  X_k_arr(i,j,k,n) = 2.0;
              else
                  X_k_arr(i,j,k,n) = 0.0;
          }
#else
         // two layers
         if (z > zhi)
            X_k_arr(i,j,k,n) = 0.0;
         else
            X_k_arr(i,j,k,n) = vec_ptr[n];
#endif
      }); 
}
