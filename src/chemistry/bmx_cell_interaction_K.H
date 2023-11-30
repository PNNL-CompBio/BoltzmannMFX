//
//     Copyright (c) 2013 Battelle Memorial Institute
//     Licensed under modified BSD License. A copy of this license can be found
//     in the LICENSE file in the top level directory of this distribution.
//
#ifndef _BMX_CELL_INTERACTION_K
#define _BMX_CELL_INTERACTION_K

#include <math.h>
#include <bmx_fluid_parms.H>
#include <bmx_cell_interaction.H>

/**
 * Calculate the maximum separation distance for which there is any
 * interaction
 * @param p1 parameters describing particle 1
 * @param p2 parameters describing particle 2
 * @param i1 integer parameters describing particle 1
 * @param i2 integer parameters describing particle 2
 * @param fpar force parameters
 * @return maximum interaction distance between particles 1 and 2
 */
AMREX_GPU_HOST_DEVICE AMREX_INLINE
Real 
maxInteractionDistance(const Real *p1, const Real *p2, const int *i1,
    const int *i2, Real *fpar)
{
  if (i1[intIdx::cell_type] == cellType::YEAST &&
      i2[intIdx::cell_type] == cellType::YEAST) {
    Real bndry_width = fpar[2];
    Real r1 = p1[realIdx::radius];
    Real r2 = p2[realIdx::radius];
    return (r1+r2+bndry_width);
  } else if ((i1[intIdx::cell_type] == cellType::FUNGI &&
      i2[intIdx::cell_type] == cellType::YEAST) ||
      (i1[intIdx::cell_type] == cellType::YEAST &&
      i2[intIdx::cell_type] == cellType::FUNGI)) {
    Real bndry_width = fpar[2];
    Real r1 = p1[realIdx::radius];
    Real r2 = p2[realIdx::radius];
    return (r1+r2+bndry_width);
  } else if (i1[intIdx::cell_type] == cellType::FUNGI &&
      i2[intIdx::cell_type] == cellType::FUNGI) {
    Real bndry_width = fpar[0];
    Real r1 = p1[realIdx::c_length];
    Real r2 = p2[realIdx::c_length];
    return (1.5*(r1+r2)+bndry_width);
  }
  return 0.0;
}

/**
 * calculate the minimum distance between two line segments defined by
 * the points (xi1,yi1,zi1) and (xi2,yi2,zi2) on segment i and
 * the points (xj1,yj1,zj1) and (xj2,yj2,zj2) on segment j
 * @param xi1,yi1,zi1,xi2,yi2,zi2 endpoints of segment i
 * @param xj1,yj1,zj1,xj2,yj2,zj2 endpoints of segment j
 * @param tau_i, tau_j location of along segments representing location
 *                     of minimum separation
 */
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void 
evaluateMinimumDistance(Real xi1, Real yi1, Real zi1,
    Real xi2, Real yi2, Real zi2, Real xj1, Real yj1, Real zj1,
    Real xj2, Real yj2, Real zj2, Real *tau_i, Real *tau_j)
{
  // Calculate segment vectors
  Real rho_i[3];
  Real rho_j[3];
  rho_i[0] = xi2-xi1;
  rho_i[1] = yi2-yi1;
  rho_i[2] = zi2-zi1;
  rho_j[0] = xj2-xj1;
  rho_j[1] = yj2-yj1;
  rho_j[2] = zj2-zj1;
  Real rirj = 0.0;

  for (int ii=0; ii<3; ii++) {
    rirj += rho_i[ii]*rho_j[ii];
  } // ii

  // Find tau_i and tau_j corresponding to minimum distance between
  // cylinder segments
  Real Aii = rho_i[0]*rho_i[0]+rho_i[1]*rho_i[1]+rho_i[2]*rho_i[2];
  Real Aij = -rirj;
  Real Aji = Aij;
  Real Ajj = rho_j[0]*rho_j[0]+rho_j[1]*rho_j[1]+rho_j[2]*rho_j[2];
  Real det = Aii*Ajj-Aij*Aji;
  Real Bi = -((xi1-xj1)*rho_i[0]+(yi1-yj1)*rho_i[1]+(zi1-zj1)*rho_i[2]);
  Real Bj = (xi1-xj1)*rho_j[0]+(yi1-yj1)*rho_j[1]+(zi1-zj1)*rho_j[2];
  // First two cases are unlikely
  bool search = false;
  if (Aij == 0.0) {
    // segments are perpendicular to each other
    *tau_i = Bi/Aii; 
    *tau_j = Bj/Ajj; 
  } else if (det == 0.0) {
    // segments are parallel to each other. Try setting tau_i to 0 and 1
    // tau_i = 0.0;
    Real tau_j0 = Bj/Ajj;
    if (tau_j0 >= 0.0 && tau_j0 <= 1.0) {
      *tau_i = 0.0;
      *tau_j = tau_j0;
    } else {
      // tau_i = 1.0;
      Real tau_j1 = (Bj+rirj)/Ajj;
      if (tau_j1 >= 0.0 && tau_j1 <= 1.0) {
        *tau_i = 1.0;
        *tau_j = tau_j1;
      } else {
        // Check all four combinations of endpoints
        search = true;
      }
    }
  } else {
    // Find minimum distance ignoring constraints
    *tau_i = (Ajj*Bi-Aij*Bj)/det;
    *tau_j = (Aii*Bj-Aji*Bi)/det;
    if (*tau_i >= 0.0 && *tau_i <= 1.0) {
      // tau_i satisfies constraints. Check tau_j
      if (*tau_j < 0.0 ||  *tau_j > 1.0) {
        // tau_j is outside limits, Need to try tau_j = 0.0 and
        // tau_j = 1.0 and find minimum separation
        Real tau_i0 = Bi/Aii;
        Real tau_i1 = (Bi+rirj)/Aii;
        if (tau_i0 >= 0.0 && tau_i0 <= 1.0 && tau_i1 >= 0.0 && tau_i1 <= 1.0) {
          // both new values of tau_i are within limits. Compare distances
          Real rix = xi1 + tau_i0*rho_i[0];
          Real riy = yi1 + tau_i0*rho_i[1];
          Real riz = zi1 + tau_i0*rho_i[2];
          Real sep0 = (rix-xj1)*(rix-xj1)+(riy-yj1)*(riy-yj1)
            +(riz-zj1)*(riz-zj1);
          rix = xi1 + tau_i1*rho_i[0];
          riy = yi1 + tau_i1*rho_i[1];
          riz = zi1 + tau_i1*rho_i[2];
          Real sep1 = (rix-xj2)*(rix-xj2)+(riy-yj2)*(riy-yj2)
            +(riz-zj2)*(riz-zj2);
          if (sep0 < sep1) {
            *tau_j = 0.0;
            *tau_i = tau_i0;
          } else {
            *tau_j = 1.0;
            *tau_i = tau_i1;
          }
        } else if (tau_i0 >= 0.0 && tau_i0 <= 1.0) {
          // only tau_i0 is within bounds
          *tau_j = 0.0;
          *tau_i = tau_i0;
        } else if (tau_i1 >= 0.0 && tau_i1 <= 1.0) {
          // only tau_i1 is within bounds
          *tau_j = 1.0;
          *tau_i = tau_i1;
        } else {
          // Neither solution is within bounds
          search = true;
        }
      }
    } else if (*tau_j >= 0.0 && *tau_j <= 1.0) {
      // already know that tao_i is out of bounds
      Real tau_j0 = Bj/Ajj;
      Real tau_j1 = (Bj+rirj)/Ajj;
      if (tau_j0 >= 0.0 && tau_j0 <= 1.0 && tau_j1 >= 0.0 && tau_j1 <= 1.0) {
        // both new values of tau_j are within limits. Compare distances
        Real rjx = xj1 + tau_j0*rho_j[0];
        Real rjy = yj1 + tau_j0*rho_j[1];
        Real rjz = zj1 + tau_j0*rho_j[2];
        Real sep0 = (rjx-xi1)*(rjx-xi1)+(rjy-yi1)*(rjy-yi1)+(rjz-zi1)*(rjz-zi1);
        rjx = xj1 + tau_j1*rho_j[0];
        rjy = yj1 + tau_j1*rho_j[1];
        rjz = zj1 + tau_j1*rho_j[2];
        Real sep1 = (rjx-xi2)*(rjx-xi2)+(rjy-yi2)*(rjy-yi2)+(rjz-zi2)*(rjz-zi2);
        if (sep0 < sep1) {
          *tau_i = 0.0;
          *tau_j = tau_j0;
        } else {
          *tau_i = 1.0;
          *tau_j = tau_j1;
        }
      } else if (tau_j0 >= 0.0 && tau_j0 <= 1.0) {
        // only tau_j0 is within bounds
        *tau_i = 0.0;
        *tau_j = tau_j0;
      } else if (tau_j1 >= 0.0 && tau_j1 <= 1.0) {
        // only tau_j1 is within bounds
        *tau_i = 1.0;
        *tau_j = tau_j1;
      } else {
        // Neither solution is within bounds
        search = true;
      }
    } else {
      // if we reach this point, both tau_i and tau_j are outside
      // the constraints. Try all four end points for minimum distance
      search = true;
    }
  }
  if (search) {
    // try all four endpoints to find minimum distance
    Real min0 = (xi1-xj1)*(xi1-xj1)+(yi1-yj1)*(yi1-yj1)+(zi1-zj1)*(zi1-zj1);
    Real min1 = (xi1-xj2)*(xi1-xj2)+(yi1-yj2)*(yi1-yj2)+(zi1-zj2)*(zi1-zj2);
    Real min2 = (xi2-xj1)*(xi2-xj1)+(yi2-yj1)*(yi2-yj1)+(zi2-zj1)*(zi2-zj1);
    Real min3 = (xi2-xj2)*(xi2-xj2)+(yi2-yj2)*(yi2-yj2)+(zi2-zj2)*(zi2-zj2);
    int imin = 0;
    Real min = min0;
    if (min1 < min) {
      imin = 1;
      min = min1;
    }
    if (min2 < min) {
      imin = 2;
      min = min2;
    }
    if (min3 < min) {
      imin = 3;
      min = min3;
    }
    if (imin == 0) {
      *tau_i = 0.0;
      *tau_j = 0.0;
    } else if (imin == 1) {
      *tau_i = 0.0;
      *tau_j = 1.0;
    } else if (imin == 2) {
      *tau_i = 1.0;
      *tau_j = 0.0;
    } else if (imin == 3) {
      *tau_i = 1.0;
      *tau_j = 1.0;
    }
  }
}

/**
 * calculate the minimum distance between point (px,py,pz)
 * and a line segment defined by the two points (x1,y1,z1)
 * and (x2,y2,z2)
 * @param px, py, pz coordinates of point
 * @param x1,y1,z1,x2,y2,z2 endpoints of segment
 * @param tau location along segment representing location
 *            of minimum separation
 */
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void 
evaluatePointSegmentMinimumDistance(Real px, Real py, Real pz,
    Real x1, Real y1, Real z1, Real x2, Real y2, Real z2,
    Real *tau)
{
  Real rho[3];
  rho[0] = x2-x1;
  rho[1] = y2-y1;
  rho[2] = z2-z1;
  Real rho2 = rho[0]*rho[0]+rho[1]*rho[1]+rho[2]*rho[2];
  Real rp = rho[0]*(px-x1)+rho[1]*(py-y1)+rho[2]*(pz-z1);
  *tau = rp/rho2;
  // If tau falls outside of allowed range [0,1], set minimum
  // distance to endpoint
  if (*tau < 0.0) *tau = 0.0;
  if (*tau > 1.0) *tau = 1.0;
}

/**
 * Calculate the force between 2 particles. r12 is computed outside this routine
 * and periodic boundary conditions, if applicable, have already been applied
 * @param r12 position of particle 2 - position of particle 1
 * @param par1 real parameters describing particle 1
 * @param par2 real parameters describing particle 2
 * @param ipar1 integer parameters describing particle 1
 * @param ipar2 integer parameters describing particle 2
 * @param vel1 velocity of center of mass of segment 1
 * @param vel2 velocity of center of mass of segment 2
 * @param rotv1 rotational velocity of segment 1
 * @param rotv2 rotational velocity of segment 2
 * @param fpar array of params used in force calculation
 */
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void 
evaluateForce(const Real *r12, Real *par1, Real *par2, int *ipar1, int *ipar2,
              Real *vel1, Real *vel2, Real *rotv1, Real *rotv2, Real *fpar, int /*me*/,
              int i_index, int j_index, int step)
{
  if (ipar1[intIdx::cell_type] == cellType::YEAST &&
      ipar2[intIdx::cell_type] == cellType::YEAST) {
    Real bndry_width = fpar[2];
    Real stiffness = fpar[3];
    Real rS = par1[realIdx::radius] + par2[realIdx::radius];
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
      F= -stiffness*(rn-rA)*(rn-rA)*(rn-rS);
    } else {
      F = 0.0;
    }

    vel1[0] = F*rx;
    vel1[1] = F*ry;
    vel1[2] = F*rz;
    vel2[0] = -F*rx;
    vel2[1] = -F*ry;
    vel2[2] = -F*rz;
    rotv1[0] = 0.0;
    rotv1[1] = 0.0;
    rotv1[2] = 0.0;
    rotv2[0] = 0.0;
    rotv2[1] = 0.0;
    rotv2[2] = 0.0;
  } else if ((ipar1[intIdx::cell_type] == cellType::YEAST &&
      ipar2[intIdx::cell_type] == cellType::FUNGI) ||
      (ipar1[intIdx::cell_type] == cellType::FUNGI &&
       ipar2[intIdx::cell_type] == cellType::YEAST)) {
    Real fungi_bndry_width = fpar[0];
    Real fungi_stiffness = fpar[1];
    Real cell_bndry_width = fpar[2];
    Real cell_stiffness = fpar[3];
    Real bndry_width = 0.5*(fungi_bndry_width+cell_bndry_width);
    Real stiffness = std::sqrt(fungi_stiffness*cell_stiffness);
    Real rS = par1[realIdx::radius] + par2[realIdx::radius];
    Real rA = rS + bndry_width;
    Real visc = fpar[12];
    Real f1[3];
    Real f2_s1[3];
    Real f2_s2[3];
    Real px, py, pz;
    Real rx, ry, rz, rn, r;
    Real x1, y1, z1, x2, y2, z2;
    bool part1_seg2;
    rx = r12[0];
    ry = r12[1];
    rz = r12[2];
    rn = sqrt(rx*rx+ry*ry+rz*rz);
    rx /= rn;
    ry /= rn;
    rz /= rn;
    int ii, jj, kk;

    if (ipar1[intIdx::cell_type] == cellType::YEAST &&
              ipar2[intIdx::cell_type] == cellType::FUNGI) {
      part1_seg2 = true;
    } else {
      part1_seg2 = false;
    }

    for (int ii=0; ii<3; ii++) {
      f1[ii] = 0.0;
      f2_s1[ii] = 0.0;
      f2_s2[ii] = 0.0;
    }

    Real c_length, theta, phi;
    if (part1_seg2) {
      c_length = par2[realIdx::c_length];
      theta = par2[realIdx::theta];
      phi = par2[realIdx::phi];
    } else {
      c_length = par1[realIdx::c_length];
      theta = par1[realIdx::theta];
      phi = par1[realIdx::phi];
    }
    Real ct = cos(theta);
    Real st = sin(theta);
    Real cp = cos(phi);
    Real sp = sin(phi);
    Real nx = st*cp;
    Real ny = st*sp;
    Real nz = ct;

    // Find interaction centers on segment 2
    if (part1_seg2) {
      x1 = -0.5*nx*c_length;
      y1 = -0.5*ny*c_length;
      z1 = -0.5*nz*c_length;
      x2 = 0.5*nx*c_length;
      y2 = 0.5*ny*c_length;
      z2 = 0.5*nz*c_length;
      px = r12[0];
      py = r12[1];
      pz = r12[2];
    } else {
      x1 = r12[0] - 0.5*nx*c_length;
      y1 = r12[1] - 0.5*ny*c_length;
      z1 = r12[2] - 0.5*nz*c_length;
      x2 = r12[0] + 0.5*nx*c_length;
      y2 = r12[1] + 0.5*ny*c_length;
      z2 = r12[2] + 0.5*nz*c_length;
      px = 0.0;
      py = 0.0;
      pz = 0.0;
    }

    // Construct rotation matrices for this segment (for later use)
    Real rot[3][3];
    rot[0][0] = cp*st;
    rot[0][1] = -sp;
    rot[0][2] = -cp*ct;
    rot[1][0] = sp*st;
    rot[1][1] = cp;
    rot[1][2] = -sp*ct;
    rot[2][0] = ct;
    rot[2][1] = 0.0;
    rot[2][2] = st;

    Real tau;
    evaluatePointSegmentMinimumDistance(px, py, pz,
        x1, y1, z1, x2, y2, z2, &tau);
    // Calculate segment vectors
    Real rho[3];
    rho[0] = x2-x1;
    rho[1] = y2-y1;
    rho[2] = z2-z1;

    // Evaluate forces and distribute forces to individual sites
    Real rt[3];
    rt[0] = x1 + tau*rho[0];
    rt[1] = y1 + tau*rho[1];
    rt[2] = z1 + tau*rho[2];
    rx = px - rt[0];
    ry = py - rt[1];
    rz = pz - rt[2];
    r = sqrt(rx*rx+ry*ry+rz*rz);

    if (r > 1.0e-8) {
      rx /= r;
      ry /= r;
      rz /= r;
    } else {
      /* Particle intersects segment (this is a corner case).
       * Pick an arbitrary separation vector perpendicular to rho */
      Real Rt = std::sqrt(rho[0]*rho[0] + rho[1]*rho[1]
          + rho[2]*rho[2]);
      for (int ii=0; ii<3; ii++) {
        rho[ii] /= Rt;
      }

      Real qx, qy, qz;
      if (rho[0] > 0.0 || rho[1] > 0.0) {
        qx = 0.0;
        qy = 0.0;
        qz = 1.0;
      } else {
        qx = 1.0;
        qy = 0.0;
        qz = 0.0;
      }

      /* Take cross product of vectors and multiply it by 1.0e-8 */
      Real rx = qy*rho[2]-qz*rho[1];
      Real ry = qz*rho[0]-qx*rho[2];
      Real rz = qx*rho[1]-qy*rho[0];
      rx *= 1.0e-08;
      ry *= 1.0e-08;
      rz *= 1.0e-08;
    }

    Real F;
    if (r < rA) {
      F= -stiffness*(r-rA)*(r-rA)*(r-rS);
    } else {
      F = 0.0;
    }
    if (F != 0.0) {
      if (part1_seg2) {
        f1[0] += F*rx;
        f1[1] += F*ry;
        f1[2] += F*rz;
        f2_s1[0] -= F*rx*(1.0-tau);
        f2_s1[1] -= F*ry*(1.0-tau);
        f2_s1[2] -= F*rz*(1.0-tau);
        f2_s2[0] -= F*rx*tau;
        f2_s2[1] -= F*ry*tau;
        f2_s2[2] -= F*rz*tau;
      } else {
        f1[0] -= F*rx;
        f1[1] -= F*ry;
        f1[2] -= F*rz;
        f2_s1[0] += F*rx*(1.0-tau);
        f2_s1[1] += F*ry*(1.0-tau);
        f2_s1[2] += F*rz*(1.0-tau);
        f2_s2[0] -= F*rx*tau;
        f2_s2[1] += F*ry*tau;
        f2_s2[2] += F*rz*tau;
      }
    }

    if (part1_seg2) {
      vel1[0] = f1[0];
      vel1[1] = f1[1];
      vel1[2] = f1[2];
      rotv1[0] = 0.0;
      rotv1[1] = 0.0;
      rotv1[2] = 0.0;
    } else {
      vel2[0] = f1[0];
      vel2[1] = f1[1];
      vel2[2] = f1[2];
      rotv2[0] = 0.0;
      rotv2[1] = 0.0;
      rotv2[2] = 0.0;
    }

    // Subtract out the part of force that is moving sites together
    // and evaluate velocity of CM and rotational velocity

    rx = x1-x2;
    ry = y1-y2;
    rz = z1-z2;
    r = sqrt(rx*rx+ry*ry+rz*rz);
    nx = rx/r;
    ny = ry/r;
    nz = rz/r;
    Real fx, fy, fz, fn;
    fx = f2_s1[0]-f2_s2[0];
    fy = f2_s1[1]-f2_s2[1];
    fz = f2_s1[2]-f2_s2[2];
    fn = fx*nx + fy*ny + fz*nz;
    f2_s1[0] -= 0.5*fn*nx;
    f2_s1[1] -= 0.5*fn*ny;
    f2_s1[2] -= 0.5*fn*nz;
    f2_s2[0] += 0.5*fn*nx;
    f2_s2[1] += 0.5*fn*ny;
    f2_s2[2] += 0.5*fn*nz;

    Real dii[3][3];
    Real dij[3][3];
    Real dtmp[3][3];
    // Calculate friction tensor for segment. Start by calculating
    // size of spheres with same exposed surface as cylinder
    Real area, arad;
    if (part1_seg2) {
      area = par2[realIdx::area];
    } else {
      area = par1[realIdx::area];
    }
    area /= 2.0;
    arad = sqrt(area/(4.0*M_PI));
    // Assume the cylinder is oriented along the x-axis
    for (int ii=0; ii<3; ii++) {
      for (int jj=0; jj<3; jj++) {
        if (ii==jj) {
          dii[ii][jj] = 1.0/(6.0*M_PI*arad*visc);
        } else {
          dii[ii][jj] = 0.0;
        }
        dij[ii][jj] = 0.0;
      } // jj
    } // ii
    if (c_length > 2.0*arad) {
      dij[0][0] = (c_length*c_length-2.0*arad*arad/3.0)
        /(4.0*M_PI*pow(c_length,3)*visc);
      dij[1][1] = (c_length*c_length+2.0*arad*arad/3.0)
        /(8.0*M_PI*pow(c_length,3)*visc);
      dij[2][2] = dij[1][1];
    } else {
      dij[0][0] = (1.0-3.0*c_length/(16.0*arad))/(6.0*M_PI*arad*visc);
      dij[1][1] = (1.0-9.0*c_length/(32.0*arad))/(6.0*M_PI*arad*visc);
      dij[2][2] = dij[1][1];
    }
    // Reorient dii and dij so that it reflects the correct direction (dii is
    // diagonal but not proportional to the identity, so changes under rotation).
    // Note the transpose of the rotation is equal to its inverse
    for (int ii=0; ii<3; ii++) {
      for (int jj=0; jj<3; jj++) {
        dtmp[ii][jj] = 0.0;
        for (int kk=0; kk<3; kk++) {
          // use inverse rotation here
          dtmp[ii][jj] += dii[ii][kk]*rot[jj][kk];
        }
      }
    }
    for (int ii=0; ii<3; ii++) {
      for (int jj=0; jj<3; jj++) {
        dii[ii][jj] = 0.0;
        for (int kk=0; kk<3; kk++) {
          dii[ii][jj] += rot[ii][kk]*dtmp[kk][jj];
        } // kk
      } // jj
    } // ii

    for (int ii=0; ii<3; ii++) {
      for (int jj=0; jj<3; jj++) {
        dtmp[ii][jj] = 0.0;
        for (int kk=0; kk<3; kk++) {
          // use inverse rotation here
          dtmp[ii][jj] += dij[ii][kk]*rot[jj][kk];
        }
      }
    }
    for (int ii=0; ii<3; ii++) {
      for (int jj=0; jj<3; jj++) {
        dij[ii][jj] = 0.0;
        for (int kk=0; kk<3; kk++) {
          dij[ii][jj] += rot[ii][kk]*dtmp[kk][jj];
        } // kk
      } // jj
    } // ii

    // Calculate velocities on sites 1 and 2
    Real v_1[3], v_2[3];
    for (int ii=0; ii<3; ii++) {
      v_1[ii] = 0.0;
      v_2[ii] = 0.0;
      for (int jj=0; jj<3; jj++) {
        v_1[ii] += dii[ii][jj]*f2_s1[jj]+dij[ii][jj]*f2_s2[jj];
        v_2[ii] += dii[ii][jj]*f2_s2[jj]+dij[ii][jj]*f2_s1[jj];
      } // jj
    } // ii

    // Get translational and rotational velocities of each segment
    if (part1_seg2) {
      vel2[0] = 0.5*(v_1[0]+v_2[0]);
      vel2[1] = 0.5*(v_1[1]+v_2[1]);
      vel2[2] = 0.5*(v_1[2]+v_2[2]);

      rotv2[0] = 0.5*(ry*v_1[2]-rz*v_1[1]);
      rotv2[1] = 0.5*(rz*v_1[0]-rx*v_1[2]);
      rotv2[2] = 0.5*(rx*v_1[1]-ry*v_1[0]);
      rotv2[0] += -0.5*(ry*v_2[2]-rz*v_2[1]);
      rotv2[1] += -0.5*(rz*v_2[0]-rx*v_2[2]);
      rotv2[2] += -0.5*(rx*v_2[1]-ry*v_2[0]);
    } else {
      vel1[0] = 0.5*(v_1[0]+v_2[0]);
      vel1[1] = 0.5*(v_1[1]+v_2[1]);
      vel1[2] = 0.5*(v_1[2]+v_2[2]);

      rotv1[0] = 0.5*(ry*v_1[2]-rz*v_1[1]);
      rotv1[1] = 0.5*(rz*v_1[0]-rx*v_1[2]);
      rotv1[2] = 0.5*(rx*v_1[1]-ry*v_1[0]);
      rotv1[0] += -0.5*(ry*v_2[2]-rz*v_2[1]);
      rotv1[1] += -0.5*(rz*v_2[0]-rx*v_2[2]);
      rotv1[2] += -0.5*(rx*v_2[1]-ry*v_2[0]);
    }

  } else if (ipar1[intIdx::cell_type] == cellType::FUNGI &&
      ipar2[intIdx::cell_type] == cellType::FUNGI) {
    // Non-bonding interactions
    Real bndry_width = fpar[0];
    Real stiffness = fpar[1];
    Real rS = par1[realIdx::radius] + par2[realIdx::radius];
    Real rA = rS + bndry_width;
    // Bonded interactions
    bool bonded = false;
    // Scale bond interaction in the case of a fusion event
    Real scale = par1[realIdx::bond_scale]*par2[realIdx::bond_scale];
    Real bond_strength = fpar[10];
    Real bond_cutoff = fpar[11];
    Real visc = fpar[12];
    Real f1_s1[3];
    Real f1_s2[3];
    Real f2_s1[3];
    Real f2_s2[3];
    int i_id[4];
    int i_cpu[4];
    int j_id[4];
    int j_cpu[4];
    int i_site[4];
    int j_site[4];
    int i_bnds = ipar1[intIdx::n_bnds];
    i_id[0] = ipar1[intIdx::seg1_id1];
    i_id[1] = ipar1[intIdx::seg2_id1];
    i_id[2] = ipar1[intIdx::seg3_id1];
    i_id[3] = ipar1[intIdx::seg4_id1];

    i_cpu[0] = ipar1[intIdx::seg1_id2];
    i_cpu[1] = ipar1[intIdx::seg2_id2];
    i_cpu[2] = ipar1[intIdx::seg3_id2];
    i_cpu[3] = ipar1[intIdx::seg4_id2];

    j_id[0] = ipar2[intIdx::seg1_id1];
    j_id[1] = ipar2[intIdx::seg2_id1];
    j_id[2] = ipar2[intIdx::seg3_id1];
    j_id[3] = ipar2[intIdx::seg4_id1];

    j_cpu[0] = ipar2[intIdx::seg1_id2];
    j_cpu[1] = ipar2[intIdx::seg2_id2];
    j_cpu[2] = ipar2[intIdx::seg3_id2];
    j_cpu[3] = ipar2[intIdx::seg4_id2];

    i_site[0] = ipar1[intIdx::site1];
    i_site[1] = ipar1[intIdx::site2];
    i_site[2] = ipar1[intIdx::site3];
    i_site[3] = ipar1[intIdx::site4];

    j_site[0] = ipar2[intIdx::site1];
    j_site[1] = ipar2[intIdx::site2];
    j_site[2] = ipar2[intIdx::site3];
    j_site[3] = ipar2[intIdx::site4];
    int j_bnds = ipar2[intIdx::n_bnds];
    int iid = ipar1[intIdx::id];
    int icpu = ipar1[intIdx::cpu];
    int jid = ipar2[intIdx::id];
    int jcpu = ipar2[intIdx::cpu];
    // Fix up status of interior location if attached TIP has converted to an
    // interior segment
    if (ipar1[intIdx::position] == siteLocation::SECOND_1 ||
        ipar1[intIdx::position] == siteLocation::SECOND_2) {
      int site;
      if (ipar1[intIdx::position] == siteLocation::SECOND_1) {
        // Growth tip is attached to site 1
        site = 1;
      } else {
        // Growth tip is attached to site 2
        site = 2;
      }
      // find out if segment j is attached to appropriate site
      int ib;
      for (ib=0; ib<i_bnds; ib++) {
        if (i_id[ib] == jid && i_cpu[ib] == jcpu && i_site[ib] == site) {
          break;
        }
      }
      // segment j is attached to segment i at site specified by
      // siteLocation value. Is segment j a growth tip? If not,
      // change position to INTERIOR
      if (ib<i_bnds && ipar2[intIdx::position] != siteLocation::TIP) {
#ifndef AMREX_USE_GPU
        printf("Converting site position on segment id: %d cpu: %d\n",iid,icpu);
#endif
        ipar1[intIdx::position] = siteLocation::INTERIOR;
      }
    }
    // clean up old flags
    ipar1[intIdx::dltd_site1] = -1;
    ipar1[intIdx::dltd_site2] = -1;
    ipar1[intIdx::dltd_id1] = -1;
    ipar1[intIdx::dltd_id2] = -1;
    ipar1[intIdx::dltd_cpu1] = -1;
    ipar1[intIdx::dltd_cpu2] = -1;
    // Is i bonded to j and is j bonded to i? Take action if only one particle
    // recognized the bond
    bool i_bonded_to_j = false;
    bool j_bonded_to_i = false;
    int iidx = -1;
    int jidx = -1;
    for (int ib = 0; ib<i_bnds; ib++) {
      if (i_id[ib] == jid && i_cpu[ib] == jcpu) {
        iidx = ib;
        i_bonded_to_j = true;
        break;
      }
    }
    for (int jb = 0; jb<j_bnds; jb++) {
      if (j_id[jb] == iid && j_cpu[jb] == icpu) {
        jidx = jb;
        j_bonded_to_i = true;
        break;
      }
    }
    if (!i_bonded_to_j && j_bonded_to_i) {
      // j has registered i as a bonding partner but i
      // has not registered j as a bonding partner. Fix up i
#ifndef AMREX_USE_GPU
      printf("I id: %d cpu: %d nbnds: %d id1: %d cpu1: %d id2: %d cpu: %d"
          " id3: %d cpu3: %d id4: %d cpu4: %d fix_site: %d\n",ipar1[intIdx::id],
          ipar1[intIdx::cpu],ipar1[intIdx::n_bnds],
          ipar1[intIdx::seg1_id1],ipar1[intIdx::seg1_id2],
          ipar1[intIdx::seg2_id1],ipar1[intIdx::seg2_id2],
          ipar1[intIdx::seg3_id1],ipar1[intIdx::seg3_id2],
          ipar1[intIdx::seg4_id1],ipar1[intIdx::seg4_id2],
          ipar1[intIdx::fix_site]);
      printf("J id: %d cpu: %d nbnds: %d id1: %d cpu1: %d id2: %d cpu: %d"
          " id3: %d cpu3: %d id4: %d cpu4: %d fix_site: %d\n",ipar2[intIdx::id],
          ipar2[intIdx::cpu],ipar2[intIdx::n_bnds],
          ipar2[intIdx::seg1_id1],ipar2[intIdx::seg1_id2],
          ipar2[intIdx::seg2_id1],ipar2[intIdx::seg2_id2],
          ipar2[intIdx::seg3_id1],ipar2[intIdx::seg3_id2],
          ipar2[intIdx::seg4_id1],ipar2[intIdx::seg4_id2],
          ipar2[intIdx::fix_site]);
#endif
      int ib0 = i_bnds;
      // Figure out where which site to add new bond
      int isite;
      if (ipar1[intIdx::fix_site]==0) {
        if (ipar1[intIdx::site1] == 1) {
          isite = 2;
        } else {
          isite = 1;
        }
      } else {
        isite = ipar1[intIdx::fix_site];
      }

#ifndef AMREX_USE_GPU
      printf("FIXING up i id: %d cpu: %d attached to id: %d cpu: %d ib: %d i_bnds: %d isite: %d\n",
          iid,icpu,jid,jcpu,ib0,i_bnds,isite);
      printf("pre FIXING up i id: %d cpu: %d nbnds: %d id1: %d cpu1: %d site1: %d"
          " id2: %d cpu2: %d site2: %d id3: %d cpu3: %d site3: %d id4: %d cpu4: %d site4: %d\n",
          ipar1[intIdx::id],ipar1[intIdx::cpu],ipar1[intIdx::n_bnds],
          ipar1[intIdx::seg1_id1],ipar1[intIdx::seg1_id2],ipar1[intIdx::site1],
          ipar1[intIdx::seg2_id1],ipar1[intIdx::seg2_id2],ipar1[intIdx::site2],
          ipar1[intIdx::seg3_id1],ipar1[intIdx::seg3_id2],ipar1[intIdx::site3],
          ipar1[intIdx::seg4_id1],ipar1[intIdx::seg4_id2],ipar1[intIdx::site4]);
#endif
      if (i_bnds == 0) {
        ipar1[intIdx::seg1_id1] = jid;
        ipar1[intIdx::seg1_id2] = jcpu;
        ipar1[intIdx::site1] = isite;
        i_site[i_bnds] = isite;
        i_id[i_bnds] = jid;
        i_cpu[i_bnds] = jcpu;
        i_bnds++;
        ipar1[intIdx::n_bnds] = i_bnds;
      } else if (i_bnds == 1) {
        ipar1[intIdx::seg2_id1] = jid;
        ipar1[intIdx::seg2_id2] = jcpu;
        ipar1[intIdx::site2] = isite;
        i_site[i_bnds] = isite;
        i_id[i_bnds] = jid;
        i_cpu[i_bnds] = jcpu;
        i_bnds++;
        ipar1[intIdx::n_bnds] = i_bnds;
      } else if (i_bnds == 2) {
        ipar1[intIdx::seg3_id1] = jid;
        ipar1[intIdx::seg3_id2] = jcpu;
        ipar1[intIdx::site3] = isite;
        i_site[i_bnds] = isite;
        i_id[i_bnds] = jid;
        i_cpu[i_bnds] = jcpu;
        i_bnds++;
        ipar1[intIdx::n_bnds] = i_bnds;
      } else if (i_bnds == 3) {
        ipar1[intIdx::seg4_id1] = jid;
        ipar1[intIdx::seg4_id2] = jcpu;
        ipar1[intIdx::site4] = isite;
        i_site[i_bnds] = isite;
        i_id[i_bnds] = jid;
        i_cpu[i_bnds] = jcpu;
        i_bnds++;
        ipar1[intIdx::n_bnds] = i_bnds;
      }
#ifndef AMREX_USE_GPU
      printf("post FIXING up i id: %d cpu: %d nbnds: %d id1: %d cpu1: %d site1: %d"
          " id2: %d cpu2: %d site2: %d id3: %d cpu3: %d site3: %d id4: %d cpu4: %d site4: %d\n",
          ipar1[intIdx::id],ipar1[intIdx::cpu],ipar1[intIdx::n_bnds],
          ipar1[intIdx::seg1_id1],ipar1[intIdx::seg1_id2],ipar1[intIdx::site1],
          ipar1[intIdx::seg2_id1],ipar1[intIdx::seg2_id2],ipar1[intIdx::site2],
          ipar1[intIdx::seg3_id1],ipar1[intIdx::seg3_id2],ipar1[intIdx::site3],
          ipar1[intIdx::seg4_id1],ipar1[intIdx::seg4_id2],ipar1[intIdx::site4]);
#endif
    } else if (i_bonded_to_j && !j_bonded_to_i) {
      // i has registered j as a bonding partner but j
      // has not registered i as a bonding partner. Fix up j
#ifndef AMREX_USE_GPU
      printf("I id: %d cpu: %d nbnds: %d id1: %d cpu1: %d id2: %d cpu: %d"
          " id3: %d cpu3: %d id4: %d cpu4: %d fix_site: %d\n",ipar1[intIdx::id],
          ipar1[intIdx::cpu],ipar1[intIdx::n_bnds],
          ipar1[intIdx::seg1_id1],ipar1[intIdx::seg1_id2],
          ipar1[intIdx::seg2_id1],ipar1[intIdx::seg2_id2],
          ipar1[intIdx::seg3_id1],ipar1[intIdx::seg3_id2],
          ipar1[intIdx::seg4_id1],ipar1[intIdx::seg4_id2],
          ipar1[intIdx::fix_site]);
      printf("J id: %d cpu: %d nbnds: %d id1: %d cpu1: %d id2: %d cpu: %d"
          " id3: %d cpu3: %d id4: %d cpu4: %d fix_site: %d\n",ipar2[intIdx::id],
          ipar2[intIdx::cpu],ipar2[intIdx::n_bnds],
          ipar2[intIdx::seg1_id1],ipar2[intIdx::seg1_id2],
          ipar2[intIdx::seg2_id1],ipar2[intIdx::seg2_id2],
          ipar2[intIdx::seg3_id1],ipar2[intIdx::seg3_id2],
          ipar2[intIdx::seg4_id1],ipar2[intIdx::seg4_id2],
          ipar2[intIdx::fix_site]);
#endif
      int jb0 = j_bnds;
      // Figure out where which site to add new bond
      int jsite;
      if (ipar2[intIdx::fix_site]==0) {
        if (ipar2[intIdx::site1] == 1) {
          jsite = 2;
        } else {
          jsite = 1;
        }
      } else {
        jsite = ipar2[intIdx::fix_site];
      }
#ifndef AMREX_USE_GPU
      printf("FIXING up j id: %d cpu: %d attached to id: %d cpu: %d jb: %d j_bnds: %d jsite: %d\n",
          jid,jcpu,iid,icpu,jb0,j_bnds,jsite);
      printf("pre FIXING up j id: %d cpu: %d nbnds: %d id1: %d cpu1: %d site1: %d"
          " id2: %d cpu2: %d site2: %d id3: %d cpu3: %d site3: %d id4: %d cpu4: %d site4: %d\n",
          ipar2[intIdx::id],ipar2[intIdx::cpu],ipar2[intIdx::n_bnds],
          ipar2[intIdx::seg1_id1],ipar2[intIdx::seg1_id2],ipar2[intIdx::site1],
          ipar2[intIdx::seg2_id1],ipar2[intIdx::seg2_id2],ipar2[intIdx::site2],
          ipar2[intIdx::seg3_id1],ipar2[intIdx::seg3_id2],ipar2[intIdx::site3],
          ipar2[intIdx::seg4_id1],ipar2[intIdx::seg4_id2],ipar2[intIdx::site4]);
#endif
      if (j_bnds == 0) {
        ipar2[intIdx::seg1_id1] = iid;
        ipar2[intIdx::seg1_id2] = icpu;
        ipar2[intIdx::site1] = jsite;
        j_site[j_bnds] = jsite;
        j_id[j_bnds] = iid;
        j_cpu[j_bnds] = icpu;
        j_bnds++;
        ipar2[intIdx::n_bnds] = j_bnds;
      } else if (j_bnds == 1) {
        ipar2[intIdx::seg2_id1] = iid;
        ipar2[intIdx::seg2_id2] = icpu;
        ipar2[intIdx::site2] = jsite;
        j_site[j_bnds] = jsite;
        j_id[j_bnds] = iid;
        j_cpu[j_bnds] = icpu;
        j_bnds++;
        ipar2[intIdx::n_bnds] = j_bnds;
      } else if (j_bnds == 2) {
        ipar2[intIdx::seg3_id1] = iid;
        ipar2[intIdx::seg3_id2] = icpu;
        ipar2[intIdx::site3] = jsite;
        j_site[j_bnds] = jsite;
        j_id[j_bnds] = iid;
        j_cpu[j_bnds] = icpu;
        j_bnds++;
        ipar2[intIdx::n_bnds] = j_bnds;
      } else if (j_bnds == 3) {
        ipar2[intIdx::seg4_id1] = iid;
        ipar2[intIdx::seg4_id2] = icpu;
        ipar2[intIdx::site4] = jsite;
        j_site[j_bnds] = jsite;
        j_id[j_bnds] = iid;
        j_cpu[j_bnds] = icpu;
        j_bnds++;
        ipar2[intIdx::n_bnds] = j_bnds;
      }
#ifndef AMREX_USE_GPU
      printf("post FIXING up j id: %d cpu: %d nbnds: %d id1: %d cpu1: %d site1: %d"
          " id2: %d cpu2: %d site2: %d id3: %d cpu3: %d site3: %d id4: %d cpu4: %d site4: %d\n",
          ipar2[intIdx::id],ipar2[intIdx::cpu],ipar2[intIdx::n_bnds],
          ipar2[intIdx::seg1_id1],ipar2[intIdx::seg1_id2],ipar2[intIdx::site1],
          ipar2[intIdx::seg2_id1],ipar2[intIdx::seg2_id2],ipar2[intIdx::site2],
          ipar2[intIdx::seg3_id1],ipar2[intIdx::seg3_id2],ipar2[intIdx::site3],
          ipar2[intIdx::seg4_id1],ipar2[intIdx::seg4_id2],ipar2[intIdx::site4]);
#endif
    }
    for (int ii=0; ii<3; ii++) {
      f1_s1[ii] = 0.0;
      f1_s2[ii] = 0.0;
      f2_s1[ii] = 0.0;
      f2_s2[ii] = 0.0;
    }
    // find locations of segment end points
    Real xi1, yi1, zi1, xi2, yi2, zi2;
    Real xj1, yj1, zj1, xj2, yj2, zj2;
    Real c_length1 = par1[realIdx::c_length];
    Real c_length2 = par2[realIdx::c_length];
    Real theta = par1[realIdx::theta];
    Real phi = par1[realIdx::phi];
    Real ct = cos(theta);
    Real st = sin(theta);
    Real cp = cos(phi);
    Real sp = sin(phi);
    Real nx = st*cp;
    Real ny = st*sp;
    Real nz = ct;
    // Find interaction centers on segment i
    xi1 = r12[0] - 0.5*nx*c_length1;
    yi1 = r12[1] - 0.5*ny*c_length1;
    zi1 = r12[2] - 0.5*nz*c_length1;
    xi2 = r12[0] + 0.5*nx*c_length1;
    yi2 = r12[1] + 0.5*ny*c_length1;
    zi2 = r12[2] + 0.5*nz*c_length1;
    // Construct rotation matrices for this segment (for later use)
    Real rot_i[3][3];
    rot_i[0][0] = cp*st;
    rot_i[0][1] = -sp;
    rot_i[0][2] = -cp*ct;
    rot_i[1][0] = sp*st;
    rot_i[1][1] = cp;
    rot_i[1][2] = -sp*ct;
    rot_i[2][0] = ct;
    rot_i[2][1] = 0.0;
    rot_i[2][2] = st;

    theta = par2[realIdx::theta];
    phi = par2[realIdx::phi];
    ct = cos(theta);
    st = sin(theta);
    cp = cos(phi);
    sp = sin(phi);
    nx = st*cp;
    ny = st*sp;
    nz = ct;
    // Construct rotation matrices for this segment
    Real rot_j[3][3];
    rot_j[0][0] = cp*st;
    rot_j[0][1] = -sp;
    rot_j[0][2] = -cp*ct;
    rot_j[1][0] = sp*st;
    rot_j[1][1] = cp;
    rot_j[1][2] = -sp*ct;
    rot_j[2][0] = ct;
    rot_j[2][1] = 0.0;
    rot_j[2][2] = st;

    // find interaction centers on segment j
    xj1 = -0.5*nx*c_length2;
    yj1 = -0.5*ny*c_length2;
    zj1 = -0.5*nz*c_length2;
    xj2 = 0.5*nx*c_length2;
    yj2 = 0.5*ny*c_length2;
    zj2 = 0.5*nz*c_length2;

    // Look for bonds connecting these two segments
    int ib, jb;
    for (ib=0; ib<i_bnds; ib++) {
      if (i_id[ib] == jid && i_cpu[ib] == jcpu) break;
    }
    for (jb=0; jb<j_bnds; jb++) {
      if (j_id[jb] == iid && j_cpu[jb] == icpu) break;
    }
    // Found a connection
    Real rx, ry, rz, r;
    if (ib <i_bnds && jb < j_bnds) {
      int itip = ipar1[intIdx::fuse_tip];
      int jtip = ipar2[intIdx::fuse_tip];
      // Four possibilities for bonded interactions
      if ((itip > 0 && itip == i_site[ib]) ||
          (jtip > 0 && jtip == j_site[jb])) {
        bond_strength = scale*bond_strength;
#if 0
        printf("SCALING check i_bnd: %d is1: %d is2: %d is3: %d is4: %d\n",
            i_bnds,i_site[0],i_site[1],i_site[2],i_site[3]);
        printf("SCALING check j_bnd: %d js1: %d js2: %d js3: %d js4: %d\n",
            j_bnds,j_site[0],j_site[1],j_site[2],j_site[3]);
        printf("SCALING interaction between %d,%d and %d,%d scale: %f"
            " ib: %d jb: %d isite: %d jsite: %d itip: %d jtip: %d\n",
            iid,icpu,jid,jcpu,scale,ib,jb,i_site[ib],j_site[jb],itip,jtip);
#endif
      }
      if ((i_site[ib] == 1 && j_site[jb] == 1)) {
        // 1-1 interaction
        rx = xi1 - xj1;
        ry = yi1 - yj1;
        rz = zi1 - zj1;
        r = sqrt(rx*rx+ry*ry+rz*rz);
        f1_s1[0] -= bond_strength*rx;
        f1_s1[1] -= bond_strength*ry;
        f1_s1[2] -= bond_strength*rz;
        f2_s2[0] += bond_strength*rx;
        f2_s2[1] += bond_strength*ry;
        f2_s2[2] += bond_strength*rz;
        bonded = true;
      } // if ((i_site[ib] == 1 && j_site[j] == 1)) ...

      if ((i_site[ib] == 1 && j_site[jb] == 2 && !bonded)) {
        // 1-2 interaction
        rx = xi1 - xj2;
        ry = yi1 - yj2;
        rz = zi1 - zj2;
        r = sqrt(rx*rx+ry*ry+rz*rz);
        f1_s1[0] -= bond_strength*rx;
        f1_s1[1] -= bond_strength*ry;
        f1_s1[2] -= bond_strength*rz;
        f2_s2[0] += bond_strength*rx;
        f2_s2[1] += bond_strength*ry;
        f2_s2[2] += bond_strength*rz;
        bonded = true;
      } // if ((i_site[i] == 1 && j_site[j] == 2 && !bonded)) ...
      if ((i_site[ib] == 2 && j_site[jb] == 1 && !bonded)) {
        // 2-1 interaction
        rx = xi2 - xj1;
        ry = yi2 - yj1;
        rz = zi2 - zj1;
        r = sqrt(rx*rx+ry*ry+rz*rz);
        f1_s2[0] -= bond_strength*rx;
        f1_s2[1] -= bond_strength*ry;
        f1_s2[2] -= bond_strength*rz;
        f2_s1[0] += bond_strength*rx;
        f2_s1[1] += bond_strength*ry;
        f2_s1[2] += bond_strength*rz;
        bonded = true;
      }
      if ((i_site[ib] == 2 && j_site[jb] == 2 && !bonded)) {
        // 2-2 interaction
        rx = xi2 - xj2;
        ry = yi2 - yj2;
        rz = zi2 - zj2;
        r = sqrt(rx*rx+ry*ry+rz*rz);
        f1_s2[0] -= bond_strength*rx;
        f1_s2[1] -= bond_strength*ry;
        f1_s2[2] -= bond_strength*rz;
        f2_s1[0] += bond_strength*rx;
        f2_s1[1] += bond_strength*ry;
        f2_s1[2] += bond_strength*rz;
        bonded = true;
      } // if ((i_site[i] == 2 && j_site[j] == 2 && !bonded)) ...
    } // if (i<i_bnds && j_bnds) ...
    // Evaluate interactions between cylinders. Only calculate these
    // interactions if cylinders are not bonded.
    if (!bonded) {
      Real tau_i, tau_j;
      evaluateMinimumDistance(xi1, yi1, zi1, xi2, yi2, zi2,
          xj1, yj1, zj1, xj2, yj2, zj2, &tau_i, &tau_j);
      // Calculate segment vectors
      Real rho_i[3];
      Real rho_j[3];
      rho_i[0] = xi2-xi1;
      rho_i[1] = yi2-yi1;
      rho_i[2] = zi2-zi1;
      rho_j[0] = xj2-xj1;
      rho_j[1] = yj2-yj1;
      rho_j[2] = zj2-zj1;
      // Evaluate forces and distribute forces to individual sites
      Real ri[3], rj[3];
      ri[0] = xi1 + tau_i*rho_i[0];
      ri[1] = yi1 + tau_i*rho_i[1];
      ri[2] = zi1 + tau_i*rho_i[2];
      rj[0] = xj1 + tau_j*rho_j[0];
      rj[1] = yj1 + tau_j*rho_j[1];
      rj[2] = zj1 + tau_j*rho_j[2];
      rx = ri[0] - rj[0];
      ry = ri[1] - rj[1];
      rz = ri[2] - rj[2];
      r = sqrt(rx*rx+ry*ry+rz*rz);
      if (r > 1.0e-8) {
        rx /= r;
        ry /= r;
        rz /= r;
      } else {
        /* Segments intersect (this is a corner case). Find an
         * appropriate separation vector */
        Real Ri = std::sqrt(rho_i[0]*rho_i[0] + rho_i[1]*rho_i[1]
            + rho_i[2]*rho_i[2]);
        Real Rj = std::sqrt(rho_j[0]*rho_j[0] + rho_j[1]*rho_j[1]
            + rho_j[2]*rho_j[2]);

        for (int ii=0; ii<3; ii++) {
          rho_i[ii] /= Ri;
          rho_j[ii] /= Rj;
        }

        /* Take cross product of vectors and multiply it by 1.0e-8 */
        Real rx = rho_i[1]*rho_j[2]-rho_i[2]*rho_j[1];
        Real ry = rho_i[2]*rho_j[0]-rho_i[0]*rho_j[2];
        Real rz = rho_i[0]*rho_j[1]-rho_i[1]*rho_j[0];
        rx *= 1.0e-08;
        ry *= 1.0e-08;
        rz *= 1.0e-08;
      }

      Real F;
      if (r < rA) {
        F= -stiffness*(r-rA)*(r-rA)*(r-rS);
      } else {
        F = 0.0;
      }
      if (F != 0.0) {
        f1_s1[0] += F*rx*(1.0-tau_i);
        f1_s1[1] += F*ry*(1.0-tau_i);
        f1_s1[2] += F*rz*(1.0-tau_i);
        f1_s2[0] += F*rx*tau_i;
        f1_s2[1] += F*ry*tau_i;
        f1_s2[2] += F*rz*tau_i;
        f2_s1[0] -= F*rx*(1.0-tau_j);
        f2_s1[1] -= F*ry*(1.0-tau_j);
        f2_s1[2] -= F*rz*(1.0-tau_j);
        f2_s2[0] -= F*rx*tau_j;
        f2_s2[1] -= F*ry*tau_j;
        f2_s2[2] -= F*rz*tau_j;
      }
    }

    // Subtract out the part of force that is moving sites together
    // and evaluate velocity of CM and rotational velocity
    Real rx_i = xi1-xi2;
    Real ry_i = yi1-yi2;
    Real rz_i = zi1-zi2;
    r = sqrt(rx_i*rx_i+ry_i*ry_i+rz_i*rz_i);
    nx = rx_i/r;
    ny = ry_i/r;
    nz = rz_i/r;
    Real fx = f1_s1[0]-f1_s2[0];
    Real fy = f1_s1[1]-f1_s2[1];
    Real fz = f1_s1[2]-f1_s2[2];
    Real fn = fx*nx + fy*ny + fz*nz;
    f1_s1[0] -= 0.5*fn*nx;
    f1_s1[1] -= 0.5*fn*ny;
    f1_s1[2] -= 0.5*fn*nz;
    f1_s2[0] += 0.5*fn*nx;
    f1_s2[1] += 0.5*fn*ny;
    f1_s2[2] += 0.5*fn*nz;

    Real rx_j = xj1-xj2;
    Real ry_j = yj1-yj2;
    Real rz_j = zj1-zj2;
    r = sqrt(rx_j*rx_j+ry_j*ry_j+rz_j*rz_j);
    nx = rx_j/r;
    ny = ry_j/r;
    nz = rz_j/r;
    fx = f2_s1[0]-f2_s2[0];
    fy = f2_s1[1]-f2_s2[1];
    fz = f2_s1[2]-f2_s2[2];
    fn = fx*nx + fy*ny + fz*nz;
    f2_s1[0] -= 0.5*fn*nx;
    f2_s1[1] -= 0.5*fn*ny;
    f2_s1[2] -= 0.5*fn*nz;
    f2_s2[0] += 0.5*fn*nx;
    f2_s2[1] += 0.5*fn*ny;
    f2_s2[2] += 0.5*fn*nz;

    // Calculate friction tensor for segment i. Start by calculating
    // size of spheres with same exposed surface as cylinder
    Real area = par1[realIdx::area];
    area /= 2.0;
    Real arad = sqrt(area/(4.0*M_PI));
    // Assume the cylinder is oriented along the x-axis
    Real dii[3][3];
    Real dij[3][3];
    for (int ii=0; ii<3; ii++) {
      for (int jj=0; jj<3; jj++) {
        if (ii==jj) {
          dii[ii][jj] = 1.0/(6.0*M_PI*arad*visc);
        } else {
          dii[ii][jj] = 0.0;
        }
        dij[ii][jj] = 0.0;
      }
    }
    if (c_length1 > 2.0*arad) {
      dij[0][0] = (c_length1*c_length1-2.0*arad*arad/3.0)
        /(4.0*M_PI*pow(c_length1,3)*visc);
      dij[1][1] = (c_length1*c_length1+2.0*arad*arad/3.0)
        /(8.0*M_PI*pow(c_length1,3)*visc);
      dij[2][2] = dij[1][1];
    } else {
      dij[0][0] = (1.0-3.0*c_length1/(16.0*arad))/(6.0*M_PI*arad*visc);
      dij[1][1] = (1.0-9.0*c_length1/(32.0*arad))/(6.0*M_PI*arad*visc);
      dij[2][2] = dij[1][1];
    }

    // Reorient dii and dij so that it reflects the correct direction (dii is
    // diagonal but not proportional to the identity, so changes under rotation).
    // Note the transpose of the rotation is equal to its inverse
    Real dtmp[3][3];
    for (int ii=0; ii<3; ii++) {
      for (int jj=0; jj<3; jj++) {
        dtmp[ii][jj] = 0.0;
        for (int kk=0; kk<3; kk++) {
          // use inverse rotation here
          dtmp[ii][jj] += dii[ii][kk]*rot_i[jj][kk];
        } // kk
      } // jj 
    } // ii

    for (int ii=0; ii<3; ii++) {
      for (int jj=0; jj<3; jj++) {
        dii[ii][jj] = 0.0;
        for (int kk=0; kk<3; kk++) {
          dii[ii][jj] += rot_i[ii][kk]*dtmp[kk][jj];
        } // kk
      } // jj
    } // ii

    for (int ii=0; ii<3; ii++) {
      for (int jj=0; jj<3; jj++) {
        dtmp[ii][jj] = 0.0;
        for (int kk=0; kk<3; kk++) {
          // use inverse rotation here
          dtmp[ii][jj] += dij[ii][kk]*rot_i[jj][kk];
        } // kk
      } // jj 
    } // ii

    for (int ii=0; ii<3; ii++) {
      for (int jj=0; jj<3; jj++) {
        dij[ii][jj] = 0.0;
        for (int kk=0; kk<3; kk++) {
          dij[ii][jj] += rot_i[ii][kk]*dtmp[kk][jj];
        } // kk
      } // jj
    } // ii

    // calculate velocities on sites 1 and 2
    Real vi_1[3], vi_2[3];
    for (int ii=0; ii<3; ii++) {
      vi_1[ii] = 0.0;
      vi_2[ii] = 0.0;
      for (int jj=0; jj<3; jj++) {
        vi_1[ii] += dii[ii][jj]*f1_s1[jj]+dij[ii][jj]*f1_s2[jj];
        vi_2[ii] += dii[ii][jj]*f1_s2[jj]+dij[ii][jj]*f1_s1[jj];
      } // jj 
    } // ii

    // Calculate friction tensor for segment j. Start by calculating
    // size of spheres with same exposed surface as cylinder
    area = par2[realIdx::area];
    area /= 2.0;
    arad = sqrt(area/(4.0*M_PI));
    // Assume the cylinder is oriented along the x-axis
    for (int ii=0; ii<3; ii++) {
      for (int jj=0; jj<3; jj++) {
        if (ii==jj) {
          dii[ii][jj] = 1.0/(6.0*M_PI*arad*visc);
        } else {
          dii[ii][jj] = 0.0;
        }
        dij[ii][jj] = 0.0;
      } // jj
    } // ii
    if (c_length2 > 2.0*arad) {
      dij[0][0] = (c_length2*c_length2-2.0*arad*arad/3.0)
        /(4.0*M_PI*pow(c_length2,3)*visc);
      dij[1][1] = (c_length2*c_length2+2.0*arad*arad/3.0)
        /(8.0*M_PI*pow(c_length2,3)*visc);
      dij[2][2] = dij[1][1];
    } else {
      dij[0][0] = (1.0-3.0*c_length2/(16.0*arad))/(6.0*M_PI*arad*visc);
      dij[1][1] = (1.0-9.0*c_length2/(32.0*arad))/(6.0*M_PI*arad*visc);
      dij[2][2] = dij[1][1];
    }
    // Reorient dii and dij so that it reflects the correct direction (dii is
    // diagonal but not proportional to the identity, so changes under rotation).
    // Note the transpose of the rotation is equal to its inverse
    for (int ii=0; ii<3; ii++) {
      for (int jj=0; jj<3; jj++) {
        dtmp[ii][jj] = 0.0;
        for (int kk=0; kk<3; kk++) {
          // use inverse rotation here
          dtmp[ii][jj] += dii[ii][kk]*rot_j[jj][kk];
        }
      }
    }
    for (int ii=0; ii<3; ii++) {
      for (int jj=0; jj<3; jj++) {
        dii[ii][jj] = 0.0;
        for (int kk=0; kk<3; kk++) {
          dii[ii][jj] += rot_j[ii][kk]*dtmp[kk][jj];
        } // kk
      } // jj
    } // ii

    for (int ii=0; ii<3; ii++) {
      for (int jj=0; jj<3; jj++) {
        dtmp[ii][jj] = 0.0;
        for (int kk=0; kk<3; kk++) {
          // use inverse rotation here
          dtmp[ii][jj] += dij[ii][kk]*rot_j[jj][kk];
        }
      }
    }
    for (int ii=0; ii<3; ii++) {
      for (int jj=0; jj<3; jj++) {
        dij[ii][jj] = 0.0;
        for (int kk=0; kk<3; kk++) {
          dij[ii][jj] += rot_j[ii][kk]*dtmp[kk][jj];
        } // kk
      } // jj
    } // ii

    // Calculate velocities on sites 1 and 2
    Real vj_1[3], vj_2[3];
    for (int ii=0; ii<3; ii++) {
      vj_1[ii] = 0.0;
      vj_2[ii] = 0.0;
      for (int jj=0; jj<3; jj++) {
        vj_1[ii] += dii[ii][jj]*f2_s1[jj]+dij[ii][jj]*f2_s2[jj];
        vj_2[ii] += dii[ii][jj]*f2_s2[jj]+dij[ii][jj]*f2_s1[jj];
      } // jj
    } // ii

    // Get translational and rotational velocities of each segment
    vel1[0] = 0.5*(vi_1[0]+vi_2[0]);
    vel1[1] = 0.5*(vi_1[1]+vi_2[1]);
    vel1[2] = 0.5*(vi_1[2]+vi_2[2]);
    vel2[0] = 0.5*(vj_1[0]+vj_2[0]);
    vel2[1] = 0.5*(vj_1[1]+vj_2[1]);
    vel2[2] = 0.5*(vj_1[2]+vj_2[2]);

    rotv1[0] = 0.5*(ry_i*vi_1[2]-rz_i*vi_1[1]);
    rotv1[1] = 0.5*(rz_i*vi_1[0]-rx_i*vi_1[2]);
    rotv1[2] = 0.5*(rx_i*vi_1[1]-ry_i*vi_1[0]);
    rotv1[0] += -0.5*(ry_i*vi_2[2]-rz_i*vi_2[1]);
    rotv1[1] += -0.5*(rz_i*vi_2[0]-rx_i*vi_2[2]);
    rotv1[2] += -0.5*(rx_i*vi_2[1]-ry_i*vi_2[0]);

    rotv2[0] = 0.5*(ry_j*vj_1[2]-rz_j*vj_1[1]);
    rotv2[1] = 0.5*(rz_j*vj_1[0]-rx_j*vj_1[2]);
    rotv2[2] = 0.5*(rx_j*vj_1[1]-ry_j*vj_1[0]);
    rotv2[0] += -0.5*(ry_j*vj_2[2]-rz_j*vj_2[1]);
    rotv2[1] += -0.5*(rz_j*vj_2[0]-rx_j*vj_2[2]);
    rotv2[2] += -0.5*(rx_j*vj_2[1]-ry_j*vj_2[0]);
  }
}

/**
 * Calculate force between growth surface and particle
 * @param pos position of particle
 * @param par Real parameters describing particle
 * @param ipar integer parameters describing particle
 * @param vel velocity of center of mass of segment
 * @param rotv rotation velocity of segment
 * @param fpar force model parameters
 * @param engine random number engine
 */
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void 
evaluateSurfaceForce(const Real *pos, Real *par, int *ipar, Real *vel, Real *rotv, Real *fpar, amrex::RandomEngine const& engine)
{
  if (ipar[intIdx::cell_type] == cellType::YEAST) {
    Real z_wall = fpar[8];
    Real z_bndry_width = fpar[6];
    Real z_stiffness = fpar[7];
    Real z_gravity = fpar[9];
    Real scale = fpar[13];
    Real mix = fpar[14];
    Real fx, fy;
    vel[0] = 0.0;
    vel[1] = 0.0;
    rotv[0] = 0.0;
    rotv[1] = 0.0;
    rotv[2] = 0.0;
    Real z = pos[2];
    Real dz = z - z_wall; 
    Real rS = par[realIdx::radius];
    Real rA = rS + z_bndry_width;
    if (dz > 0.0 && dz < rA) {
      vel[2] = -z_stiffness*(dz-rA)*(dz-rA)*(dz-rS);
    } else if (dz <= 0) {
      vel[2] = z_stiffness*rA*rA*rS;
    } else {
      vel[2] = 0.0;
    }
    vel[2] -= z_gravity*z_stiffness*rA*rA*rS;
    // Only add random fluctuations in xy plane
    Real vx, vy;
    vx = amrex::RandomNormal(0.0,scale,engine);
    vy = amrex::RandomNormal(0.0,scale,engine);
    fx = (1.0-mix)*par[25]+mix*vx;
    fy = (1.0-mix)*par[26]+mix*vy;
    vel[0] = fx;
    vel[1] = fy;
    par[25] = fx;
    par[26] = fy;
    par[27] = 0.0;
  } else if (ipar[intIdx::cell_type] == cellType::FUNGI) {
    Real z_wall = fpar[8];
    Real z_bndry_width = fpar[4];
    Real z_stiffness = fpar[5];
    Real z_gravity = fpar[9];
    Real visc = fpar[12];
    Real f_s1[3];
    Real f_s2[3];
#if 0
    vel[0] = 0.0;
    vel[1] = 0.0;
    vel[2] = 0.0;
    rotv[0] = 0.0;
    rotv[1] = 0.0;
    rotv[2] = 0.0;
    return;
#endif
#if 0
    // increment bond scaling parameter
    if (par[realIdx::bond_scale] < 1.0) {
      printf("SCALING is less than 1.0 for id: %d cpu: %d scale: %f inc: %f\n",
          ipar[intIdx::id],ipar[intIdx::cpu],par[realIdx::bond_scale],fpar[9]);
      par[realIdx::bond_scale] += fpar[9];
      if (par[realIdx::bond_scale] > 1.0)
        par[realIdx::bond_scale] = 1.0;
    }
    if (par[realIdx::bond_scale] >= 1.0)
      ipar[intIdx::fuse_tip] = 0;
#endif
    ipar[intIdx::fix_site] = 0;

    int i, j, k;
    for (i=0; i<3; i++) {
      f_s1[i] = 0.0;
      f_s2[i] = 0.0;
    }
    // find locations of segment end points
    Real x1, y1, z1, x2, y2, z2;
    Real c_length = par[realIdx::c_length];
    Real theta = par[realIdx::theta];
    Real phi = par[realIdx::phi];
    Real ct = cos(theta);
    Real st = sin(theta);
    Real cp = cos(phi);
    Real sp = sin(phi);
    Real nx = st*cp;
    Real ny = st*sp;
    Real nz = ct;
    // Find interaction centers on segment
    x1 = pos[0] - 0.5*nx*c_length;
    y1 = pos[1] - 0.5*ny*c_length;
    z1 = pos[2] - 0.5*nz*c_length;
    x2 = pos[0] + 0.5*nx*c_length;
    y2 = pos[1] + 0.5*ny*c_length;
    z2 = pos[2] + 0.5*nz*c_length;
    //DBG printf("LOCATION z1: %f z2: %f\n",z1,z2);
    // Construct rotation matrices for this segment (for later use)
    Real rot[3][3];
    rot[0][0] = cp*st;
    rot[0][1] = -sp;
    rot[0][2] = -cp*ct;
    rot[1][0] = sp*st;
    rot[1][1] = cp;
    rot[1][2] = -sp*ct;
    rot[2][0] = ct;
    rot[2][1] = 0.0;
    rot[2][2] = st;
    // Calculate force on sites 1 and 2
    Real z = z1;
    Real dz = z - z_wall; 
    Real rS = par[realIdx::radius];
    Real rA = rS + z_bndry_width;
#if 1
    if (dz > 0.0 && dz < rA) {
      f_s1[2] = -z_stiffness*(dz-rA)*(dz-rA)*(dz-rS);
    } else if (dz <= 0.0) {
      f_s1[2] = z_stiffness*rA*rA*rS;
    } else {
      f_s1[2] = 0.0;
    }
    f_s1[2] -= z_gravity*z_stiffness*rA*rA*rS;
    z = z2;
    dz = z - z_wall; 
    if (dz > 0.0 && dz < rA) {
      f_s2[2] = -z_stiffness*(dz-rA)*(dz-rA)*(dz-rS);
    } else if (dz <= 0.0) {
      f_s2[2] = z_stiffness*rA*rA*rS;
    } else {
      f_s2[2] = 0.0;
    }
    f_s2[2] -= z_gravity*z_stiffness*rA*rA*rS;
#else
    if (dz > 0.0) {
      f_s1[2] = -z_stiffness*rA*rA*rS;
    } else {
      f_s1[2] = 0.99*z_stiffness*rA*rA*rS;
    }
    z = z2;
    dz = z - z_wall; 
    if (dz > 0.0) {
      f_s2[2] = -z_stiffness*rA*rA*rS;
    } else {
      f_s2[2] = 0.99*z_stiffness*rA*rA*rS;
    }
#endif
    //DBG printf("DELTA rS: %f rA: %f stffns: %f dz1: %f dz2: %f\n",rS,rA,z_stiffness,z1-z_wall,z2-z_wall);

    // Subtract out the part of force that is moving sites together
    // and evaluate velocity of CM and rotational velocity
    Real rx = x1-x2;
    Real ry = y1-y2;
    Real rz = z1-z2;
    Real r = sqrt(rx*rx+ry*ry+rz*rz);
    //DBG printf("SEPARATION: %e\n",r);
    nx = rx/r;
    ny = ry/r;
    nz = rz/r;
    Real fx = f_s1[0]-f_s2[0];
    Real fy = f_s1[1]-f_s2[1];
    Real fz = f_s1[2]-f_s2[2];
    Real fn = fx*nx + fy*ny + fz*nz;
    //DBG printf("FORCE f1x: %e f1y: %e f1z: %e f2x: %e f2y: %e f2z: %e\n",
    //DBG  f_s1[0],f_s1[1],f_s1[2],f_s2[0],f_s2[1],f_s2[2]);
    f_s1[0] -= 0.5*fn*nx;
    f_s1[1] -= 0.5*fn*ny;
    f_s1[2] -= 0.5*fn*nz;
    f_s2[0] += 0.5*fn*nx;
    f_s2[1] += 0.5*fn*ny;
    f_s2[2] += 0.5*fn*nz;
    //DBG printf("N FORCE f1x: %e f1y: %e f1z: %e f2x: %e f2y: %e f2z: %e\n",
    //DBG  f_s1[0],f_s1[1],f_s1[2],f_s2[0],f_s2[1],f_s2[2]);
#if 0
    {
      Real df = 0.0;
      Real nh[3];
      nh[0] = nx;
      nh[1] = ny;
      nh[2] = nz;
      for (i=0; i<3; i++) {
	df += (f_s1[i]-f_s2[i])*nh[i];
      }
      printf("DIFFERENCE: %e\n",df);
    }
#endif

#if 1
    // Calculate friction tensor for segment. Start by calculating
    // size of spheres with same exposed surface as cylinder
    Real area = par[realIdx::area];
    area /= 2.0;
    Real arad = sqrt(area/(4.0*M_PI));
    // Start by assuming the cylinder is oriented along the x-axis
    Real dii[3][3];
    Real dij[3][3];
    for (i=0; i<3; i++) {
      for (j=0; j<3; j++) {
	if (i==j) {
	  dii[i][j] = 1.0/(6.0*M_PI*arad*visc);
	} else {
	  dii[i][j] = 0.0;
	}
	dij[i][j] = 0.0;
      }
    }
    if (c_length > 2.0*arad) {
      dij[0][0] = (c_length*c_length-2.0*arad*arad/3.0)/(4.0*M_PI*pow(c_length,3)*visc);
      dij[1][1] = (c_length*c_length+2.0*arad*arad/3.0)/(8.0*M_PI*pow(c_length,3)*visc);
      dij[2][2] = dij[1][1];
    } else {
      dij[0][0] = (1.0-3.0*c_length/(16.0*arad))/(6.0*M_PI*arad*visc);
      dij[1][1] = (1.0-9.0*c_length/(32.0*arad))/(6.0*M_PI*arad*visc);
      dij[2][2] = dij[1][1];
    }
    // Reorient dij so that it reflects the correct direction (dii is diagonal
    // so it doesn't change under reorientation). Note the transpose of the
    // rotation is equal to its inverse
    Real dtmp[3][3];
    for (i=0; i<3; i++) {
      for (j=0; j<3; j++) {
	dtmp[i][j] = 0.0;
	for (k=0; k<3; k++) {
	  // use inverse rotation here
	  dtmp[i][j] += dij[i][k]*rot[j][k];
	}
      }
    }
    for (i=0; i<3; i++) {
      for (j=0; j<3; j++) {
	dij[i][j] = 0.0;
	for (k=0; k<3; k++) {
	  dij[i][j] += rot[i][k]*dtmp[k][j];
	}
      }
    }
#if 0
    printf("Dii11: %e Dii12: %e Dii13: %e\n",dii[0][0],dii[0][1],dii[0][2]);
    printf("Dii21: %e Dii22: %e Dii23: %e\n",dii[1][0],dii[1][1],dii[1][2]);
    printf("Dii31: %e Dii32: %e Dii33: %e\n",dii[2][0],dii[2][1],dii[2][2]);
    printf("\n");
    printf("Dij11: %e Dij12: %e Dij13: %e\n",dij[0][0],dij[0][1],dij[0][2]);
    printf("Dij21: %e Dij22: %e Dij23: %e\n",dij[1][0],dij[1][1],dij[1][2]);
    printf("Dij31: %e Dij32: %e Dij33: %e\n",dij[2][0],dij[2][1],dij[2][2]);
#endif
    // calculate velocities on sites 1 and 2
    Real v_1[3], v_2[3];
    for (i=0; i<3; i++) {
      v_1[i] = 0.0;
      v_2[i] = 0.0;
      for (j=0; j<3; j++) {
	v_1[i] += dii[i][j]*f_s1[j]+dij[i][j]*f_s2[j];
	v_2[i] += dii[i][j]*f_s2[j]+dij[i][j]*f_s1[j];
      }
    } 
    //DBG printf("V1X: %e V1Y: %e V1Z: %e V2X: %e V2Y: %e V2Z: %e\n",
    //DBG v_1[0],v_1[1],v_1[2],v_2[0],v_2[1],v_2[2]);
    // Get translational and rotational velocities for segment
    vel[0] = 0.5*(v_1[0]+v_2[0]);
    vel[1] = 0.5*(v_1[1]+v_2[1]);
    vel[2] = 0.5*(v_1[2]+v_2[2]);
    //DBG printf("VELOCITY VX: %e VY: %e VZ: %e\n",vel[0],vel[1],vel[2]);

    // Remove force on center of mass
    for (i=0; i<3; i++) {
      v_1[i] -= vel[i];
      v_2[i] -= vel[i];
    }
    // Calculate angular momentum
    rotv[0] = 0.5*(ry*v_1[2]-rz*v_1[1]);
    rotv[1] = 0.5*(rz*v_1[0]-rx*v_1[2]);
    rotv[2] = 0.5*(rx*v_1[1]-ry*v_1[0]);
    rotv[0] += -0.5*(ry*v_2[2]-rz*v_2[1]);
    rotv[1] += -0.5*(rz*v_2[0]-rx*v_2[2]);
    rotv[2] += -0.5*(rx*v_2[1]-ry*v_2[0]);
    //DBG printf("ANG MOM LX: %e LY: %e LZ: %e\n",rotv[0],rotv[1],rotv[2]);
#endif
    vel[0] = 0.5*(f_s1[0]+f_s2[0]);
    vel[1] = 0.5*(f_s1[1]+f_s2[1]);
    vel[2] = 0.5*(f_s1[2]+f_s2[2]);

    // Remove force on center of mass
    for (i=0; i<3; i++) {
      f_s1[i] -= 0.5*vel[i];
      f_s2[i] -= 0.5*vel[i];
    }

    rotv[0] = 0.5*(ry*f_s1[2]-rz*f_s1[1]);
    rotv[1] = 0.5*(rz*f_s1[0]-rx*f_s1[2]);
    rotv[2] = 0.5*(rx*f_s1[1]-ry*f_s1[0]);
    rotv[0] += -0.5*(ry*f_s2[2]-rz*f_s2[1]);
    rotv[1] += -0.5*(rz*f_s2[0]-rx*f_s2[2]);
    rotv[2] += -0.5*(rx*f_s2[1]-ry*f_s2[0]);
  }
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void 
fixBonds(Real *r12, Real *ripar, int *ipar, Real *rjpar, int *jpar, Real rlim)
{
  if (ipar[intIdx::cell_type] == cellType::FUNGI &&
      jpar[intIdx::cell_type] == cellType::FUNGI) {
    bool bonded = false;
    int *i_id;
    int *i_cpu;
    int *j_id;
    int *j_cpu;
    int *i_site;
    int *j_site;
    int i_bnds = ipar[intIdx::n_bnds];
    int j_bnds = jpar[intIdx::n_bnds];

    i_id = &ipar[intIdx::seg1_id1];
    i_cpu = &ipar[intIdx::seg1_id2];
    i_site = &ipar[intIdx::site1];

    j_id = &jpar[intIdx::seg1_id1];
    j_cpu = &jpar[intIdx::seg1_id2];
    j_site = &jpar[intIdx::site1];

    int iid = ipar[intIdx::id];
    int icpu = ipar[intIdx::cpu];
    int jid = jpar[intIdx::id];
    int jcpu = jpar[intIdx::cpu];

    // find locations of segment end points
    Real xi1, yi1, zi1, xi2, yi2, zi2;
    Real xj1, yj1, zj1, xj2, yj2, zj2;
    Real c_length1 = ripar[realIdx::c_length];
    Real c_length2 = rjpar[realIdx::c_length];

    // calculate orientation of segment i
    Real theta = ripar[realIdx::theta];
    Real phi = ripar[realIdx::phi];
    Real ct = cos(theta);
    Real st = sin(theta);
    Real cp = cos(phi);
    Real sp = sin(phi);
    Real nx = st*cp;
    Real ny = st*sp;
    Real nz = ct;

    // Find interaction centers on segment i
    xi1 = r12[0] - 0.5*nx*c_length1;
    yi1 = r12[1] - 0.5*ny*c_length1;
    zi1 = r12[2] - 0.5*nz*c_length1;
    xi2 = r12[0] + 0.5*nx*c_length1;
    yi2 = r12[1] + 0.5*ny*c_length1;
    zi2 = r12[2] + 0.5*nz*c_length1;

    // calculate orientation of segment j
    theta = rjpar[realIdx::theta];
    phi = rjpar[realIdx::phi];
    ct = cos(theta);
    st = sin(theta);
    cp = cos(phi);
    sp = sin(phi);
    nx = st*cp;
    ny = st*sp;
    nz = ct;

    // find interaction centers on segment j
    xj1 = -0.5*nx*c_length2;
    yj1 = -0.5*ny*c_length2;
    zj1 = -0.5*nz*c_length2;
    xj2 = 0.5*nx*c_length2;
    yj2 = 0.5*ny*c_length2;
    zj2 = 0.5*nz*c_length2;

    // check four possible combinations to see if anyone is bonded
    Real rx, ry, rz, r;
    // site 1 - site 1
    rx = xi1 - xj1;
    ry = yi1 - yj1;
    rz = zi1 - zj1;
    r = sqrt(rx*rx+ry*ry+rz*rz);
    printf("R11: %e\n",r);
    if (r < rlim) {
      // Add bond to i and j
      i_id[i_bnds] = jid;
      i_cpu[i_bnds] = jcpu;
      i_site[i_bnds] = 1;
      ipar[intIdx::n_bnds] = i_bnds+1;
      j_id[j_bnds] = iid;
      j_cpu[j_bnds] = icpu;
      j_site[j_bnds] = 1;
      jpar[intIdx::n_bnds] = j_bnds+1;
      bonded = true;
#ifndef AMREX_USE_GPU
      printf("Adding bond between I: %d,%d and J: %d,%d\n",iid,icpu,jid,jcpu);
#endif
    }
    if (!bonded) {
      // site 2 - site 1
      rx = xi2 - xj1;
      ry = yi2 - yj1;
      rz = zi2 - zj1;
      r = sqrt(rx*rx+ry*ry+rz*rz);
      if (r < rlim) {
        // Add bond to i and j
        i_id[i_bnds] = jid;
        i_cpu[i_bnds] = jcpu;
        i_site[i_bnds] = 2;
        ipar[intIdx::n_bnds] = i_bnds+1;
        j_id[j_bnds] = iid;
        j_cpu[j_bnds] = icpu;
        j_site[j_bnds] = 1;
        jpar[intIdx::n_bnds] = j_bnds+1;
        bonded = true;
#ifndef AMREX_USE_GPU
      printf("Adding bond between I: %d,%d and J: %d,%d\n",iid,icpu,jid,jcpu);
#endif
      }
    }
    if (!bonded) {
      // site 1 - site 2
      rx = xi1 - xj2;
      ry = yi1 - yj2;
      rz = zi1 - zj2;
      r = sqrt(rx*rx+ry*ry+rz*rz);
      if (r < rlim) {
        // Add bond to i and j
        i_id[i_bnds] = jid;
        i_cpu[i_bnds] = jcpu;
        i_site[i_bnds] = 1;
        ipar[intIdx::n_bnds] = i_bnds+1;
        j_id[j_bnds] = iid;
        j_cpu[j_bnds] = icpu;
        j_site[j_bnds] = 2;
        jpar[intIdx::n_bnds] = j_bnds+1;
        bonded = true;
#ifndef AMREX_USE_GPU
      printf("Adding bond between I: %d,%d and J: %d,%d\n",iid,icpu,jid,jcpu);
#endif
      }
    }
    if (!bonded) {
      // site 1 - site 2
      rx = xi2 - xj2;
      ry = yi2 - yj2;
      rz = zi2 - zj2;
      r = sqrt(rx*rx+ry*ry+rz*rz);
      if (r < rlim) {
        // Add bond to i and j
        i_id[i_bnds] = jid;
        i_cpu[i_bnds] = jcpu;
        i_site[i_bnds] = 2;
        ipar[intIdx::n_bnds] = i_bnds+1;
        j_id[j_bnds] = iid;
        j_cpu[j_bnds] = icpu;
        j_site[j_bnds] = 2;
        jpar[intIdx::n_bnds] = j_bnds+1;
        bonded = true;
#ifndef AMREX_USE_GPU
      printf("Adding bond between I: %d,%d and J: %d,%d\n",iid,icpu,jid,jcpu);
#endif
      }
    }
  }
}
#endif
