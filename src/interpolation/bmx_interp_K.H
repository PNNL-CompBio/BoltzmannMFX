//
//     Copyright (c) 2013 Battelle Memorial Institute
//     Licensed under modified BSD License. A copy of this license can be found
//     in the LICENSE file in the top level directory of this distribution.
//
#ifndef _BMX_INTERP_K_H_
#define _BMX_INTERP_K_H_

#include "bmx_pc.H"

/****************************************************************
 *                                                              *
 *                                                              *
 *                                                              *
 *                                                              *
 ***************************************************************/
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void trilinear_interp (const amrex::RealVect& pos,
                       amrex::Real* interp_vals,
                       amrex::Array4<amrex::Real const> const& interp_array,
                       const amrex::RealVect& plo,
                       const amrex::RealVect& dxi,
                       const int interp_comp)
{

    const amrex::RealVect lx = (pos - plo)*dxi + 0.5;
    const amrex::IntVect ijk = lx.floor();

    int i = ijk[0]; int j = ijk[1]; int k = ijk[2];

    // Weights
    const amrex::RealVect sx_hi = lx - ijk;
    const amrex::RealVect sx_lo = 1 - sx_hi;

    for (int n = 0; n < interp_comp; n++)
       interp_vals[n] = sx_lo[0]*sx_lo[1]*sx_lo[2]*interp_array(i-1, j-1, k-1,n) +
                        sx_lo[0]*sx_lo[1]*sx_hi[2]*interp_array(i-1, j-1, k  ,n) +
                        sx_lo[0]*sx_hi[1]*sx_lo[2]*interp_array(i-1, j  , k-1,n) +
                        sx_lo[0]*sx_hi[1]*sx_hi[2]*interp_array(i-1, j  , k  ,n) +
                        sx_hi[0]*sx_lo[1]*sx_lo[2]*interp_array(i  , j-1, k-1,n) +
                        sx_hi[0]*sx_lo[1]*sx_hi[2]*interp_array(i  , j-1, k  ,n) +
                        sx_hi[0]*sx_hi[1]*sx_lo[2]*interp_array(i  , j  , k-1,n) +
                        sx_hi[0]*sx_hi[1]*sx_hi[2]*interp_array(i  , j  , k  ,n);
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void one_to_one_interp (const amrex::RealVect& pos,
                        amrex::Real* interp_vals,
                        amrex::Array4<amrex::Real const> const& interp_array,
                        const amrex::RealVect& plo,
                        const amrex::RealVect& dxi,
                        const int interp_comp)
{

    const amrex::RealVect lx = (pos - plo)*dxi;
    const amrex::IntVect ijk = lx.floor();

    int i = ijk[0]; int j = ijk[1]; int k = ijk[2];

    for (int n = 0; n < interp_comp; n++)
       interp_vals[n] = interp_array(i,j,k,n);
}
#endif
