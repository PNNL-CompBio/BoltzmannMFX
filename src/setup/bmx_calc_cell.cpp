//
//     Copyright (c) 2013 Battelle Memorial Institute
//     Licensed under modified BSD License. A copy of this license can be found
//     in the LICENSE file in the top level directory of this distribution.
//
//vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
//                                                                      !
//  Subroutine: calc_cell_ic                                            !
//  Purpose: calculate the i, j or k cell index for IC regions.         !
//                                                                      !
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

#include <AMReX_Math.H>
#include <bmx_calc_cell.H>


using namespace amrex;

void calc_cell_ic(const Real dx,
                  const Real dy,
                  const Real dz,
                  const amrex::Real* lo,
                  const amrex::Real* hi,
                  const amrex::Real* plo,
                  int& i_w,
                  int& i_e,
                  int& j_s,
                  int& j_n,
                  int& k_b,
                  int& k_t)
{
  i_w = amrex::Math::floor((lo[0]-plo[0])/dx + .5);
  i_e = amrex::Math::floor((hi[0]-plo[0])/dx + .5) - 1;

  j_s = amrex::Math::floor((lo[1]-plo[1])/dy + .5);
  j_n = amrex::Math::floor((hi[1]-plo[1])/dy + .5) - 1;

  k_b = amrex::Math::floor((lo[2]-plo[2])/dz + .5);
  k_t = amrex::Math::floor((hi[2]-plo[2])/dz + .5) - 1;
}
