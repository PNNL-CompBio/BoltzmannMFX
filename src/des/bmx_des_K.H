//
//     Copyright (c) 2013 Battelle Memorial Institute
//     Licensed under modified BSD License. A copy of this license can be found
//     in the LICENSE file in the top level directory of this distribution.
//
#include "bmx_pc.H"

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void cross_product (const amrex::RealVect& a,
                    const amrex::RealVect& b,
                    amrex::RealVect& c)
{
    c[0] = a[1]*b[2] - a[2]*b[1];
    c[1] = a[2]*b[0] - a[0]*b[2];
    c[2] = a[0]*b[1] - a[1]*b[0];
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
amrex::Real dot_product (const amrex::RealVect& a,
                         const amrex::RealVect& b)
{
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void cfrelvel (const amrex::RealVect& vel1, 
               const amrex::RealVect& vel2,
               const amrex::Real& radius1,
               const amrex::Real& radius2,
               const amrex::RealVect& omega1,
               const amrex::RealVect& omega2,
               amrex::Real& vrn, 
               amrex::RealVect& vslip,
               const amrex::RealVect& normal, 
               const amrex::Real dist_mag)
{
    // translational relative velocity
    amrex::RealVect vreltrans(0.);
    vreltrans[0] =  vel1[0] - vel2[0];
    vreltrans[1] =  vel1[1] - vel2[1];
    vreltrans[2] =  vel1[2] - vel2[2];

    // distance from the contact point to the particle centers
    amrex::Real dist_c1 = (dist_mag*dist_mag 
                          + radius1*radius1 - radius2*radius2) / (2.0*dist_mag);

    amrex::Real dist_c2 = dist_mag - dist_c1;

    amrex::RealVect omega_sum(0.);
    omega_sum[0] = omega1[0] * dist_c1 + omega2[0] * dist_c2;
    omega_sum[1] = omega1[1] * dist_c1 + omega2[1] * dist_c2;
    omega_sum[2] = omega1[2] * dist_c1 + omega2[2] * dist_c2;

    // the rotational relative velocity
    amrex::RealVect v_rot(0.);
    cross_product(omega_sum, normal, v_rot);

    // total relative velocity
    vreltrans[0] += v_rot[0];
    vreltrans[1] += v_rot[1];
    vreltrans[2] += v_rot[2];

    // normal component of relative velocity
    vrn = dot_product(vreltrans, normal);

    vslip[0] = vreltrans[0] - vrn*normal[0];
    vslip[1] = vreltrans[1] - vrn*normal[1];
    vslip[2] = vreltrans[2] - vrn*normal[2];    
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
amrex::Real interp_level_set (const amrex::RealVect& pos,
                              const int n_refine,
                              amrex::Array4<amrex::Real const> const& phi,
                              amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& plo,
                              amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxi)
{
    amrex::RealVect scaled_dxi(0.);
    scaled_dxi[0] = n_refine * dxi[0];
    scaled_dxi[1] = n_refine * dxi[1];
    scaled_dxi[2] = n_refine * dxi[2];
    
    amrex::Real x = (pos[0] - plo[0]) * scaled_dxi[0];
    amrex::Real y = (pos[1] - plo[1]) * scaled_dxi[1];
    amrex::Real z = (pos[2] - plo[2]) * scaled_dxi[2];

    int i = static_cast<int>(amrex::Math::floor(x));
    int j = static_cast<int>(amrex::Math::floor(y));
    int k = static_cast<int>(amrex::Math::floor(z));
    
    amrex::Real wx_hi = x - i;
    amrex::Real wy_hi = y - j;
    amrex::Real wz_hi = z - k;
    
    amrex::Real wx_lo = 1.0 - wx_hi;
    amrex::Real wy_lo = 1.0 - wy_hi;
    amrex::Real wz_lo = 1.0 - wz_hi;

    amrex::Real ls_value = 0;
    ls_value += phi(i,   j  , k  ) * wx_lo * wy_lo * wz_lo;
    ls_value += phi(i+1, j  , k  ) * wx_hi * wy_lo * wz_lo;
    ls_value += phi(i,   j+1, k  ) * wx_lo * wy_hi * wz_lo;
    ls_value += phi(i+1, j+1, k  ) * wx_hi * wy_hi * wz_lo;
    ls_value += phi(i,   j  , k+1) * wx_lo * wy_lo * wz_hi;
    ls_value += phi(i+1, j  , k+1) * wx_hi * wy_lo * wz_hi;
    ls_value += phi(i  , j+1, k+1) * wx_lo * wy_hi * wz_hi;
    ls_value += phi(i+1, j+1, k+1) * wx_hi * wy_hi * wz_hi;
    
    return ls_value;
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void level_set_normal (const amrex::RealVect& pos,
                       const int n_refine,
                       amrex::RealVect& normal,
                       amrex::Array4<amrex::Real const> const& phi,
                       amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& plo,
                       amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxi)
{
    amrex::RealVect scaled_dxi(0.);
    scaled_dxi[0] = n_refine * dxi[0];
    scaled_dxi[1] = n_refine * dxi[1];
    scaled_dxi[2] = n_refine * dxi[2];
    
    amrex::Real x = (pos[0] - plo[0]) * scaled_dxi[0];
    amrex::Real y = (pos[1] - plo[1]) * scaled_dxi[1];
    amrex::Real z = (pos[2] - plo[2]) * scaled_dxi[2];

    int i = static_cast<int>(amrex::Math::floor(x));
    int j = static_cast<int>(amrex::Math::floor(y));
    int k = static_cast<int>(amrex::Math::floor(z));
    
    amrex::Real wx_hi = x - i;
    amrex::Real wy_hi = y - j;
    amrex::Real wz_hi = z - k;
    
    amrex::Real wx_lo = 1.0 - wx_hi;
    amrex::Real wy_lo = 1.0 - wy_hi;
    amrex::Real wz_lo = 1.0 - wz_hi;

    normal[0] = 0.0;
    normal[1] = 0.0;
    normal[2] = 0.0;

    normal[0] -= phi(i,   j  , k  ) * scaled_dxi[0] * wy_lo * wz_lo;
    normal[0] += phi(i+1, j  , k  ) * scaled_dxi[0] * wy_lo * wz_lo;
    normal[0] -= phi(i,   j+1, k  ) * scaled_dxi[0] * wy_hi * wz_lo;
    normal[0] += phi(i+1, j+1, k  ) * scaled_dxi[0] * wy_hi * wz_lo;
    normal[0] -= phi(i,   j  , k+1) * scaled_dxi[0] * wy_lo * wz_hi;
    normal[0] += phi(i+1, j  , k+1) * scaled_dxi[0] * wy_lo * wz_hi;
    normal[0] -= phi(i  , j+1, k+1) * scaled_dxi[0] * wy_hi * wz_hi;
    normal[0] += phi(i+1, j+1, k+1) * scaled_dxi[0] * wy_hi * wz_hi;

    normal[1] -= phi(i,   j  , k  ) * scaled_dxi[1] * wx_lo * wz_lo;
    normal[1] += phi(i  , j+1, k  ) * scaled_dxi[1] * wx_lo * wz_lo;
    normal[1] -= phi(i+1, j  , k  ) * scaled_dxi[1] * wx_hi * wz_lo;
    normal[1] += phi(i+1, j+1, k  ) * scaled_dxi[1] * wx_hi * wz_lo;
    normal[1] -= phi(i,   j  , k+1) * scaled_dxi[1] * wx_lo * wz_hi;
    normal[1] += phi(i  , j+1, k+1) * scaled_dxi[1] * wx_lo * wz_hi;
    normal[1] -= phi(i+1, j  , k+1) * scaled_dxi[1] * wx_hi * wz_hi;
    normal[1] += phi(i+1, j+1, k+1) * scaled_dxi[1] * wx_hi * wz_hi;

    normal[2] -= phi(i  , j  , k  ) * scaled_dxi[2] * wx_lo * wy_lo;
    normal[2] += phi(i  , j  , k+1) * scaled_dxi[2] * wx_lo * wy_lo;
    normal[2] -= phi(i+1, j  , k  ) * scaled_dxi[2] * wx_hi * wy_lo;
    normal[2] += phi(i+1, j  , k+1) * scaled_dxi[2] * wx_hi * wy_lo;
    normal[2] -= phi(i,   j+1, k  ) * scaled_dxi[2] * wx_lo * wy_hi;
    normal[2] += phi(i  , j+1, k+1) * scaled_dxi[2] * wx_lo * wy_hi;
    normal[2] -= phi(i+1, j+1, k  ) * scaled_dxi[2] * wx_hi * wy_hi;
    normal[2] += phi(i+1, j+1, k+1) * scaled_dxi[2] * wx_hi * wy_hi;
    
    amrex::Real inv_norm = 1.0 / sqrt(dot_product(normal, normal));
    normal[0] *= inv_norm;
    normal[1] *= inv_norm;
    normal[2] *= inv_norm;
}
