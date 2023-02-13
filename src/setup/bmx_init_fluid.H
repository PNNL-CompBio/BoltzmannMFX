//
//     Copyright (c) 2013 Battelle Memorial Institute
//     Licensed under modified BSD License. A copy of this license can be found
//     in the LICENSE file in the top level directory of this distribution.
//
#ifndef __BMX_INIT_FLUID_HPP_
#define __BMX_INIT_FLUID_HPP_

#include <bmx.H>

#include <AMReX_REAL.H>
#include <AMReX_BLFort.H>
#include <AMReX_SPACE.H>
#include <AMReX_GpuContainers.H>

void set_ic_chem_species (const amrex::Box& sbx, const amrex::Box& domain,
                          const amrex::Real dx, const amrex::Real dy, const amrex::Real dz,
                          const amrex::GpuArray<amrex::Real, 3>& p_lo,
                          const amrex::GpuArray<amrex::Real, 3>& p_hi,
                          const amrex::Array4<amrex::Real> X_k_arr);
#endif
