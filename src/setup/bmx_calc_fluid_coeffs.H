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
               FArrayBox& D_k_fab);
#endif
