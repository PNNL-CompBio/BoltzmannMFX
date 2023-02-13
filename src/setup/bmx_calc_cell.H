//
//     Copyright (c) 2013 Battelle Memorial Institute
//     Licensed under modified BSD License. A copy of this license can be found
//     in the LICENSE file in the top level directory of this distribution.
//
#ifndef _BMX_CALCCELL_H_
#define _BMX_CALCCELL_H_

#include <AMReX_REAL.H>

void calc_cell_ic(const amrex::Real dx,
                  const amrex::Real dy,
                  const amrex::Real dz,
                  const amrex::Real* lo,
                  const amrex::Real* hi,
                  const amrex::Real* plo,
                  int& i_w,
                  int& i_e,
                  int& j_s,
                  int& j_n,
                  int& k_b,
                  int& k_t);

#endif
