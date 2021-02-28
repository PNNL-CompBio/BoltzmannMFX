#ifndef _BMX_CALC_FLUID_COEFFS_HPP_
#define _BMX_CALC_FLUID_COEFFS_HPP_

#include <bmx.H>

#include <AMReX_REAL.H>
#include <AMReX_BLFort.H>
#include <AMReX_SPACE.H>

// Compute fluid species diffusivity
void calc_D_gk (const Box& bx,
               FArrayBox& D_gk_fab);
#endif
