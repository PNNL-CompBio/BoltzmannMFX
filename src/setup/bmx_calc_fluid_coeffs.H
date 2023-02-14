//
//     Copyright (c) 2013 Battelle Memorial Institute
//     Licensed under modified BSD License. A copy of this license can be found
//     in the LICENSE file in the top level directory of this distribution.
//
#ifndef _BMX_CALC_FLUID_COEFFS_HPP_
#define _BMX_CALC_FLUID_COEFFS_HPP_

#include <bmx.H>

#include <AMReX_REAL.H>
#include <AMReX_BLFort.H>
#include <AMReX_SPACE.H>

// Compute fluid chem_species diffusivity
void calc_D_k (const Box& bx, const Box& domain,
               const Real dx,
               const Real dy,
               const Real dz,
               const amrex::Array4<amrex::Real> D_k_arr);
#endif
