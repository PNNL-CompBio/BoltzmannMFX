//
//     Copyright (c) 2013 Battelle Memorial Institute
//     Licensed under modified BSD License. A copy of this license can be found
//     in the LICENSE file in the top level directory of this distribution.
//
#include <math.h>
#include <bmx_fluid_parms.H>
#include <bmx_cell_interaction.H>

/**
 * Calculate the maximum separation distance for which there is any
 * interaction
 * @param p1 parameters describing particle 1
 * @param p2 parameters describing particle 2
 * @return maximum interatction distance between particles 1 and 2
 */
AMREX_GPU_HOST_DEVICE AMREX_INLINE
Real 
maxInteractionDistance(const Real *p1, const Real *p2, Real bndry_width)
{
  Real r1 = p1[realIdx::a_size];
  Real r2 = p2[realIdx::a_size];
  return (r1+r2+bndry_width);
}
/**
 * Calculate the force between 2 particles. r12 is computed outside this routine
 * and periodic boundary conditions, if applicable, have already been applied
 * @param r12 position of particle 2 - position of particle 1
 * @param par1 parameters describing particle 1
 * @param par2 parameters describing particle 2
 * @param frc force between particle 1 and particle 2
 */
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void 
evaluateForce(const Real *r12, Real *par1, Real *par2, Real *frc,
              Real bndry_width, Real stiffness )
{
  Real rS = par1[realIdx::a_size] + par2[realIdx::a_size];
  Real rA = rS + bndry_width;
  Real rx, ry, rz, rn;
  rx = r12[0];
  ry = r12[1];
  rz = r12[2];
  rn = sqrt(rx*rx+ry*ry+rz*rz);
  rx /= rn;
  ry /= rn;
  rz /= rn;

  Real F;
  if (rn < rA) {
    F= stiffness*(rn-rA)*(rn-rA)*(rn-rS);
  } else {
    F = 0.0;
  }

  frc[0] = F*rx;
  frc[1] = F*ry;
  frc[2] = F*rz;
}

/**
 * Calculate force between growth surface and particle
 * @param pos position of particle
 * @param par parameters describing particle
 * @param frc force on particle from wall
 */
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void 
evaluateSurfaceForce(const Real *pos, Real *par, Real *frc,
                     Real z_bndry_width, Real z_stiffness, Real z_wall, Real z_gravity)
{
  frc[0] = 0.0;
  frc[1] = 0.0;
  Real z = pos[2];
  Real dz = z - z_wall; 
  Real rS = par[realIdx::a_size];
  Real rA = rS + z_bndry_width;
  if (dz > 0.0 && dz < rA) {
#if 1
    frc[2] = -z_stiffness*(dz-rA)*(dz-rA)*(dz-rS);
#else
    if (dz < rS) {
      frc[2] = -p_z_r_stiffness*(dz-rA)*(dz-rA)*(dz-rS);
    } else if (dz < rA) {
      frc[2] = -p_z_a_stiffness*(dz-rA)*(dz-rA)*(dz-rS);
    } else {
      frc[2] = 0.0;
    }
#endif
  } else if (dz <= 0) {
#if 1
    frc[2] = z_stiffness*rA*rA*rS;
#else
    frc[2] = p_z_r_stiffness*rA*rA*rS;
#endif
  } else {
    frc[2] = 0.0;
  }
#if 1
  frc[2] -= z_gravity*z_stiffness*rA*rA*rS;
#else
  frc[2] -= p_z_gravity*p_z_r_stiffness*rA*rA*rS;
#endif
}
