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
 * @return maximum interatction distance between particles 1 and 2
 */
AMREX_GPU_HOST_DEVICE AMREX_INLINE
Real 
maxInteractionDistance(const Real *p1, const Real *p2, const int *i1,
    const int *i2, Real bndry_width)
{
  if (i1[intIdx::cell_type] == cellType::YEAST &&
      i2[intIdx::cell_type] == cellType::YEAST) {
    Real r1 = p1[realIdx::radius];
    Real r2 = p2[realIdx::radius];
    return (r1+r2+bndry_width);
  } else if (i1[intIdx::cell_type] == cellType::FUNGI &&
      i2[intIdx::cell_type] == cellType::FUNGI) {
    Real r1 = p1[realIdx::c_length];
    Real r2 = p2[realIdx::c_length];
    return (r1+r2+bndry_width);
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
    Real bndry_width = fpar[0];
    Real stiffness = fpar[1];
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
      F= stiffness*(rn-rA)*(rn-rA)*(rn-rS);
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
  } else if (ipar1[intIdx::cell_type] == cellType::FUNGI &&
      ipar2[intIdx::cell_type] == cellType::FUNGI) {
    // Non-bonding interactions
    Real bndry_width = fpar[0];
    Real stiffness = fpar[1];
    Real rS = par1[realIdx::radius] + par2[realIdx::radius];
    Real rA = rS + bndry_width;
    // Bonded interactions
    bool bonded = false;
    Real bond_strength = fpar[6];
    Real bond_cutoff = fpar[7];
    Real visc = fpar[8];
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
    //printf("PAIR iid: %d icpu: %d jid: %d jcpu: %d\n",iid,icpu,jid,jcpu);
#if 0
    printf("BONDS on segment id: %d cpu: %d\n",iid,icpu);
    for (i=0; i<i_bnds; i++) {
      printf("Bonded to id: %d cpu: %d at site %d\n",
          i_id[i],i_cpu[i],i_site[i]);
    }
    printf("BONDS on segment id: %d cpu: %d\n",jid,jcpu);
    for (i=0; i<j_bnds; i++) {
      printf("Bonded to id: %d cpu: %d at site %d\n",
          j_id[i],j_cpu[i],j_site[i]);
    }
#endif
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
        printf("Converting site position on segment id: %d cpu: %d\n",iid,icpu);
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
#if 1
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
      printf("I id: %d cpu: %d nbnds: %d id1: %d cpu1: %d id2: %d cpu: %d"
          " id3: %d cpu3: %d id4: %d cpu4: %d\n",ipar1[intIdx::id],
          ipar1[intIdx::cpu],ipar1[intIdx::n_bnds],
          ipar1[intIdx::seg1_id1],ipar1[intIdx::seg1_id2],
          ipar1[intIdx::seg2_id1],ipar1[intIdx::seg2_id2],
          ipar1[intIdx::seg3_id1],ipar1[intIdx::seg3_id2],
          ipar1[intIdx::seg4_id1],ipar1[intIdx::seg4_id2]);
      printf("J id: %d cpu: %d nbnds: %d id1: %d cpu1: %d id2: %d cpu: %d"
          " id3: %d cpu3: %d id4: %d cpu4: %d\n",ipar2[intIdx::id],
          ipar2[intIdx::cpu],ipar2[intIdx::n_bnds],
          ipar2[intIdx::seg1_id1],ipar2[intIdx::seg1_id2],
          ipar2[intIdx::seg2_id1],ipar2[intIdx::seg2_id2],
          ipar2[intIdx::seg3_id1],ipar2[intIdx::seg3_id2],
          ipar2[intIdx::seg4_id1],ipar2[intIdx::seg4_id2]);
      if (i_bnds == 1) {
        ipar1[intIdx::seg2_id1] = jid;
        ipar1[intIdx::seg2_id2] = jcpu;
        if (j_site[jidx] == 1) {
          ipar1[intIdx::site2] = 2;
          i_site[i_bnds] = 2;
        } else {
          ipar1[intIdx::site2] = 1;
          i_site[i_bnds] = 1;
        }
        i_id[i_bnds] = jid;
        i_cpu[i_bnds] = jcpu;
        i_bnds++;
        ipar1[intIdx::n_bnds] = i_bnds;
      } else if (i_bnds == 2) {
        ipar1[intIdx::seg3_id1] = jid;
        ipar1[intIdx::seg3_id2] = jcpu;
        if (j_site[jidx] == 1) {
          ipar1[intIdx::site3] = 2;
          i_site[i_bnds] = 2;
        } else {
          ipar1[intIdx::site3] = 1;
          i_site[i_bnds] = 1;
        }
        i_id[i_bnds] = jid;
        i_cpu[i_bnds] = jcpu;
        i_bnds++;
        ipar1[intIdx::n_bnds] = i_bnds;
      } else if (i_bnds == 3) {
        ipar1[intIdx::seg4_id1] = jid;
        ipar1[intIdx::seg4_id2] = jcpu;
        if (j_site[jidx] == 1) {
          ipar1[intIdx::site4] = 2;
          i_site[i_bnds] = 2;
        } else {
          ipar1[intIdx::site4] = 1;
          i_site[i_bnds] = 1;
        }
        i_id[i_bnds] = jid;
        i_cpu[i_bnds] = jcpu;
        i_bnds++;
        ipar1[intIdx::n_bnds] = i_bnds;
      }
      printf("FIXING up i id: %d cpu: %d attached to id: %d cpu: %d ipar: %p\n",
          iid,icpu,jid,jcpu,ipar1);
    } else if (i_bonded_to_j && !j_bonded_to_i) {
      // i has registered j as a bonding partner but j
      // has not registered i as a bonding partner. Fix up j
      printf("I id: %d cpu: %d nbnds: %d id1: %d cpu1: %d id2: %d cpu: %d"
          " id3: %d cpu3: %d id4: %d cpu4: %d\n",ipar1[intIdx::id],
          ipar1[intIdx::cpu],ipar1[intIdx::n_bnds],
          ipar1[intIdx::seg1_id1],ipar1[intIdx::seg1_id2],
          ipar1[intIdx::seg2_id1],ipar1[intIdx::seg2_id2],
          ipar1[intIdx::seg3_id1],ipar1[intIdx::seg3_id2],
          ipar1[intIdx::seg4_id1],ipar1[intIdx::seg4_id2]);
      printf("J id: %d cpu: %d nbnds: %d id1: %d cpu1: %d id2: %d cpu: %d"
          " id3: %d cpu3: %d id4: %d cpu4: %d\n",ipar2[intIdx::id],
          ipar2[intIdx::cpu],ipar2[intIdx::n_bnds],
          ipar2[intIdx::seg1_id1],ipar2[intIdx::seg1_id2],
          ipar2[intIdx::seg2_id1],ipar2[intIdx::seg2_id2],
          ipar2[intIdx::seg3_id1],ipar2[intIdx::seg3_id2],
          ipar2[intIdx::seg4_id1],ipar2[intIdx::seg4_id2]);
      if (j_bnds == 1) {
        ipar2[intIdx::seg2_id1] = iid;
        ipar2[intIdx::seg2_id2] = icpu;
        if (j_site[iidx] == 1) {
          ipar2[intIdx::site2] = 2;
          j_site[j_bnds] = 2;
        } else {
          ipar2[intIdx::site2] = 1;
          j_site[j_bnds] = 1;
        }
        j_id[j_bnds] = iid;
        j_cpu[j_bnds] = icpu;
        j_bnds++;
        ipar2[intIdx::n_bnds] = j_bnds;
      } else if (j_bnds == 2) {
        ipar2[intIdx::seg3_id1] = iid;
        ipar2[intIdx::seg3_id2] = icpu;
        if (j_site[iidx] == 1) {
          ipar2[intIdx::site3] = 2;
          j_site[j_bnds] = 2;
        } else {
          ipar2[intIdx::site3] = 1;
          j_site[j_bnds] = 1;
        }
        j_id[j_bnds] = iid;
        j_cpu[j_bnds] = icpu;
        j_bnds++;
        ipar2[intIdx::n_bnds] = j_bnds;
      } else if (j_bnds == 3) {
        ipar2[intIdx::seg4_id1] = iid;
        ipar2[intIdx::seg4_id2] = icpu;
        if (j_site[iidx] == 1) {
          ipar2[intIdx::site4] = 2;
          j_site[j_bnds] = 2;
        } else {
          ipar2[intIdx::site4] = 1;
          j_site[j_bnds] = 1;
        }
        j_id[j_bnds] = iid;
        j_cpu[j_bnds] = icpu;
        j_bnds++;
        ipar2[intIdx::n_bnds] = j_bnds;
      }
      printf("FIXING up j id: %d cpu: %d attached to id: %d cpu: %d jpar: %p\n",
          jid,jcpu,iid,icpu,ipar2);
    }
#endif
//    printf("p[%d] ID: %d CPU: %d I_ID: %d %d %d I_CPU: %d %d %d\n",
//        me,iid,icpu,i_id[0],i_id[1],i_id[2],i_cpu[0],
//        i_cpu[1],i_cpu[2]);
//    printf("p[%d] nbnd: %d i_site[0]: %d i_site[1]: %d i_site[2]: %d\n",me,i_bnds,i_site[0],i_site[1],i_site[2]);
//    printf("p[%d] nbnd: %d j_site[0]: %d j_site[1]: %d j_site[2]: %d\n",me,j_bnds,j_site[0],j_site[1],j_site[2]);
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
//    printf("p[%d] id: %d cpu: %d xi1: %e yi1: %e zi1: %e xi2: %e yi2: %e zi2: %e\n",
//      me,ipar1[intIdx::id],ipar1[intIdx::cpu],xi1,yi1,zi1,xi2,yi2,zi2);
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
//    printf("RI11: %e RI12: %e RI13: %e\n",rot_i[0][0],rot_i[0][1],rot_i[0][2]);
//    printf("RI21: %e RI22: %e RI23: %e\n",rot_i[1][0],rot_i[1][1],rot_i[1][2]);
//    printf("RI31: %e RI32: %e RI33: %e\n",rot_i[2][0],rot_i[2][1],rot_i[2][2]);

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
//    printf("RJ11: %e RJ12: %e RJ13: %e\n",rot_j[0][0],rot_j[0][1],rot_j[0][2]);
//    printf("RJ21: %e RJ22: %e RJ23: %e\n",rot_j[1][0],rot_j[1][1],rot_j[1][2]);
//    printf("RJ31: %e RJ32: %e RJ33: %e\n",rot_j[2][0],rot_j[2][1],rot_j[2][2]);

    // find interaction centers on segment j
    xj1 = -0.5*nx*c_length2;
    yj1 = -0.5*ny*c_length2;
    zj1 = -0.5*nz*c_length2;
    xj2 = 0.5*nx*c_length2;
    yj2 = 0.5*ny*c_length2;
    zj2 = 0.5*nz*c_length2;
//     printf("p[%d] id: %d cpu: %d xj1: %e yj1: %e zj1: %e xj2: %e yj2: %e zj2: %e\n",
//       me,ipar2[intIdx::id],ipar2[intIdx::cpu],xj1,yj1,zj1,xj2,yj2,zj2);

    // printf("p[%d] r12x: %e r12y: %e r12z: %e r12: %e clen1: %e clen2: %e\n",
    //   me,r12[0],r12[1],r12[2],sqrt(r12[0]*r12[0]+r12[1]*r12[1]+r12[2]*r12[2]),
    //   c_length1,c_length2);
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
      // Four possibilities for bonded interactions
      if ((i_site[ib] == 1 && j_site[jb] == 1)) {
        // 1-1 interaction
        rx = xi1 - xj1;
        ry = yi1 - yj1;
        rz = zi1 - zj1;
        r = sqrt(rx*rx+ry*ry+rz*rz);
        if (r<bond_cutoff) {
          // Might be bonded. Check further
//          printf("Found bond iid: %d icpu: %d jid: %d jcpu: %d\n",iid,icpu,jid,jcpu);
          f1_s1[0] -= bond_strength*rx;
          f1_s1[1] -= bond_strength*ry;
          f1_s1[2] -= bond_strength*rz;
          f2_s2[0] += bond_strength*rx;
          f2_s2[1] += bond_strength*ry;
          f2_s2[2] += bond_strength*rz;
          //          printf("p[%d] f1_s1: %e x1: %e f2_s2: %e x2: %e rx: %e r: %e\n",me,f1_s1[0],xi1,f2_s2[0],xj2,rx,r);
          bonded = true;
        } // r < bond_cutoff
      } // if ((i_site[ib] == 1 && j_site[j] == 1)) ...

      if ((i_site[ib] == 1 && j_site[jb] == 2 && !bonded)) {
        // 1-2 interaction
        rx = xi1 - xj2;
        ry = yi1 - yj2;
        rz = zi1 - zj2;
        r = sqrt(rx*rx+ry*ry+rz*rz);
        //    printf("p[%d] 12 Bond_cutoff: %e r: %e\n",me,bond_cutoff,r);
        if (r<bond_cutoff) {
          // Might be bonded. Check further
          //      printf("p[%d] 12 n_bonds: %d i_site[0]: %d i_id[0]: %d i_cpu[0]: %d j_id: %d j_cpu: %d\n",
          //         me,i_bnds,i_site[0],iid,icpu,jid,jcpu);
          //      printf("p[%d] DIFF (12) rx: %e ry: %e rz: %e r: %e i: %d j: %d\n",me,rx,ry,rz,r,i,j);
//          printf("Found bond iid: %d icpu: %d jid: %d jcpu: %d\n",iid,icpu,jid,jcpu);
          f1_s1[0] -= bond_strength*rx;
          f1_s1[1] -= bond_strength*ry;
          f1_s1[2] -= bond_strength*rz;
          f2_s2[0] += bond_strength*rx;
          f2_s2[1] += bond_strength*ry;
          f2_s2[2] += bond_strength*rz;
          //          printf("p[%d] f1_s1: %e x1: %e f2_s2: %e x2: %e rx: %e r: %e\n",me,f1_s1[0],xi1,f2_s2[0],xj2,rx,r);
          bonded = true;
        } // r < bond_cutoff
      } // if ((i_site[i] == 1 && j_site[j] == 2 && !bonded)) ...
      if ((i_site[ib] == 2 && j_site[jb] == 1 && !bonded)) {
        // 2-1 interaction
        rx = xi2 - xj1;
        ry = yi2 - yj1;
        rz = zi2 - zj1;
        r = sqrt(rx*rx+ry*ry+rz*rz);
        //    printf("p[%d] 21 Bond_cutoff: %e r: %e\n",me,bond_cutoff,r);
        if (r<bond_cutoff) {
          // Might be bonded. Check further
          //      printf("p[%d] 21 n_bonds: %d i_site[0]: %d i_id[0]: %d i_cpu[0]: %d j_id: %d j_cpu: %d\n",
          //         me,i_bnds,i_site[0],iid,icpu,jid,jcpu);
//          printf("p[%d] DIFF (21) rx: %e ry: %e rz: %e r: %e i: %d j: %d ibnd: %d jbnd: %d\n",
//              me,rx,ry,rz,r,i,j,i_bnds,j_bnds);
//          printf("Found bond iid: %d icpu: %d jid: %d jcpu: %d\n",iid,icpu,jid,jcpu);
          f1_s2[0] -= bond_strength*rx;
          f1_s2[1] -= bond_strength*ry;
          f1_s2[2] -= bond_strength*rz;
          f2_s1[0] += bond_strength*rx;
          f2_s1[1] += bond_strength*ry;
          f2_s1[2] += bond_strength*rz;
          //          printf("p[%d] f1_s2: %e x1: %e f2_s1: %e x2: %e rx: %e r: %e\n",me,f1_s2[0],xi2,f2_s1[0],xj1,rx,r);
          bonded = true;
        }
      }
      if ((i_site[ib] == 2 && j_site[jb] == 2 && !bonded)) {
        // 2-2 interaction
        rx = xi2 - xj2;
        ry = yi2 - yj2;
        rz = zi2 - zj2;
        r = sqrt(rx*rx+ry*ry+rz*rz);
//        printf("p[%d] 22 Bond_cutoff: %e r: %e\n",me,bond_cutoff,r);
        if (r<bond_cutoff) {
          // Might be bonded. Check further
          //      printf("p[%d] 22 n_bonds: %d i_site[0]: %d i_id[0]: %d i_cpu[0]: %d j_id: %d j_cpu: %d\n",
          //         me,i_bnds,i_site[0],iid,icpu,jid,jcpu);
//          printf("p[%d] DIFF (22) rx: %e ry: %e rz: %e r: %e i: %d j: %d ibnd: %d jbnd: %d\n",
//              me,rx,ry,rz,r,i,j,i_bnds,j_bnds);
//          printf("Found bond iid: %d icpu: %d jid: %d jcpu: %d\n",iid,icpu,jid,jcpu);
          f1_s2[0] -= bond_strength*rx;
          f1_s2[1] -= bond_strength*ry;
          f1_s2[2] -= bond_strength*rz;
          f2_s1[0] += bond_strength*rx;
          f2_s1[1] += bond_strength*ry;
          f2_s1[2] += bond_strength*rz;
          bonded = true;
        } // if (r<bond_cutoff)  ...
      } // if ((i_site[i] == 2 && j_site[j] == 2 && !bonded)) ...
    } // if (i<i_bnds && j_bnds) ...
#if 1
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
#if 0
      Real rirj = 0.0;
//    Real rtest = sqrt(r12[0]*r12[0]+r12[1]*r12[1]+r12[2]*r12[0]);
//    if (rtest > 0.003) {
//      printf("p[%d] INTERACTION PAIR IID: %d ICPU: %d JID: %d JCPU: %d\n",
//          me,iid,icpu,jid,jcpu);
//      printf("p[%d] SEP RX: %e RY: %e RZ: %e R: %e RA: %e\n",me,r12[0],r12[1],r12[2],
//          rtest,rA);
//    }

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
        tau_i = Bi/Aii; 
        tau_j = Bj/Ajj; 
      } else if (det == 0.0) {
        // segments are parallel to each other. Try setting tau_i to 0 and 1
        // tau_i = 0.0;
        Real tau_j0 = Bj/Ajj;
        if (tau_j0 >= 0.0 && tau_j0 <= 1.0) {
          tau_i = 0.0;
          tau_j = tau_j0;
        } else {
          // tau_i = 1.0;
          Real tau_j1 = (Bj+rirj)/Ajj;
          if (tau_j1 >= 0.0 && tau_j1 <= 1.0) {
            tau_i = 1.0;
            tau_j = tau_j1;
          } else {
            // Check all four combinations of endpoints
            search = true;
          }
        }
      } else {
        // Find minimum distance ignoring constraints
        tau_i = (Ajj*Bi-Aij*Bj)/det;
        tau_j = (Aii*Bj-Aji*Bi)/det;
        if (tau_i >= 0.0 && tau_i <= 1.0) {
          // tau_i satisfies constraints. Check tau_j
          if (tau_j < 0.0 ||  tau_j > 1.0) {
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
                tau_j = 0.0;
                tau_i = tau_i0;
              } else {
                tau_j = 1.0;
                tau_i = tau_i1;
              }
            } else if (tau_i0 >= 0.0 && tau_i0 <= 1.0) {
              // only tau_i0 is within bounds
              tau_j = 0.0;
              tau_i = tau_i0;
            } else if (tau_i1 >= 0.0 && tau_i1 <= 1.0) {
              // only tau_i1 is within bounds
              tau_j = 1.0;
              tau_i = tau_i1;
            } else {
              // Neither solution is within bounds
              search = true;
            }
          }
        } else if (tau_j >= 0.0 && tau_j <= 1.0) {
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
              tau_i = 0.0;
              tau_j = tau_j0;
            } else {
              tau_i = 1.0;
              tau_j = tau_j1;
            }
          } else if (tau_j0 >= 0.0 && tau_j0 <= 1.0) {
            // only tau_j0 is within bounds
            tau_i = 0.0;
            tau_j = tau_j0;
          } else if (tau_j1 >= 0.0 && tau_j1 <= 1.0) {
            // only tau_j1 is within bounds
            tau_i = 1.0;
            tau_j = tau_j1;
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
          tau_i = 0.0;
          tau_j = 0.0;
        } else if (imin == 1) {
          tau_i = 0.0;
          tau_j = 1.0;
        } else if (imin == 2) {
          tau_i = 1.0;
          tau_j = 0.0;
        } else if (imin == 3) {
          tau_i = 1.0;
          tau_j = 1.0;
        }
      }
#endif
      // Evaluate minimum distance and calculate forces. Then
      // distribute forces to individual sites
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
//      if (rtest > 0.003) {
//      printf("p[%d] rx: %e ry: %e rz: %e r: %e tau_i: %f tau_j: %f\n",
//          me,rx,ry,rz,r,tau_i,tau_j);
//      }
      if (r > 1.0e-8) {
        rx /= r;
        ry /= r;
        rz /= r;
      } else {
        /* Segments intersect. Find an appropriate separation vector */
        Real Ri = std::sqrt(rho_i[0]*rho_i[0] + rho_i[1]*rho_i[1] + rho_i[2]*rho_i[2]);
        Real Rj = std::sqrt(rho_j[0]*rho_j[0] + rho_j[1]*rho_j[1] + rho_j[2]*rho_j[2]);

        for (int ii=0; ii<3; ii++) {
          rho_i[ii] /= Ri;
          rho_j[ii] /= Rj;
        }

        Real rdot = rho_i[0]*rho_j[0] + rho_i[1]*rho_j[1] + rho_i[2]*rho_j[2];
        if (rdot < 0.0) {
          for (int ii=0; ii<3; ii++) {
            rho_j[ii] = -rho_j[ii];
          } // ii
        } // if rdot
        rx = rho_i[0]-rho_j[0];
        ry = rho_i[1]-rho_j[1];
        rz = rho_i[2]-rho_j[2];
        Real inv_rn = 1. / sqrt(rx*rx+ry*ry+rz*rz);
        rx *= inv_rn;
        ry *= inv_rn;
        rz *= inv_rn;
      }

      Real F;
      if (r < rA) {
        F= -stiffness*(r-rA)*(r-rA)*(r-rS);
      } else {
        F = 0.0;
      }
//      if (F!= 0.0) {
//      printf("p[%d] F: %e r: %e tau_i: %f tau_j: %f rA: %f rS: %f stiff: %f\n",
//          me,F,r,tau_i,tau_j,rA,rS,stiffness);
//      }
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
//      printf("fi1x: %e fi1y: %e fi1z: %e fi2x: %e fi2y: %e fi2z: %e\n",
//          f1_s1[0],f1_s1[1],f1_s1[2],f1_s2[0],f1_s2[1],f1_s1[2]);
//      printf("fj1x: %e fj1y: %e fj1z: %e fj2x: %e fj2y: %e fj2z: %e\n",
//          f2_s1[0],f2_s1[1],f2_s1[2],f2_s2[0],f2_s2[1],f2_s1[2]);
    }
    /*
    else {
      printf("BONDED PAIR IID: %d ICPU: %d JID: %d JCPU: %d\n",
          iid,icpu,jid,jcpu);
    }
    */
#endif

    // Subtract out the part of force that is moving sites together
    // and evaluate velocity of CM and rotational velocity
    Real rx_i = xi1-xi2;
    Real ry_i = yi1-yi2;
    Real rz_i = zi1-zi2;
    r = sqrt(rx_i*rx_i+ry_i*ry_i+rz_i*rz_i);
//    printf("rx_i: %f ry_i: %f rz_i: %f r: %f\n",rx_i,ry_i,rz_i,r);
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
//    printf("SEGI bond_strength: %f f1x: %e f1y: %e f1z: %e f2x: %e f2y: %e f2z: %e\n",
//      bond_strength,f1_s1[0],f1_s1[1],f1_s1[2],f1_s2[0],f1_s2[1],f1_s2[2]);

    Real rx_j = xj1-xj2;
    Real ry_j = yj1-yj2;
    Real rz_j = zj1-zj2;
    r = sqrt(rx_j*rx_j+ry_j*ry_j+rz_j*rz_j);
//    printf("rx_j: %f ry_j: %f rz_j: %f r: %f\n",rx_j,ry_j,rz_j,r);
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
//    printf("SEGJ bond_strength: %f f1x: %e f1y: %e f1z: %e f2x: %e f2y: %e f2z: %e\n",
//	bond_strength,f2_s1[0],f2_s1[1],f2_s1[2],f2_s2[0],f2_s2[1],f2_s2[2]);

#if 1
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
      dij[0][0] = (c_length1*c_length1-2.0*arad*arad/3.0)/(4.0*M_PI*pow(c_length1,3)*visc);
      dij[1][1] = (c_length1*c_length1+2.0*arad*arad/3.0)/(8.0*M_PI*pow(c_length1,3)*visc);
      dij[2][2] = dij[1][1];
    } else {
      dij[0][0] = (1.0-3.0*c_length1/(16.0*arad))/(6.0*M_PI*arad*visc);
      dij[1][1] = (1.0-9.0*c_length1/(32.0*arad))/(6.0*M_PI*arad*visc);
      dij[2][2] = dij[1][1];
    }
#if 0
    printf("(1) Dii11: %e Dii12: %e Dii13: %e\n",dii[0][0],dii[0][1],dii[0][2]);
    printf("(1) Dii21: %e Dii22: %e Dii23: %e\n",dii[1][0],dii[1][1],dii[1][2]);
    printf("(1) Dii31: %e Dii32: %e Dii33: %e\n",dii[2][0],dii[2][1],dii[2][2]);
    printf("\n");
    printf("(1) Dij11: %e Dij12: %e Dij13: %e\n",dij[0][0],dij[0][1],dij[0][2]);
    printf("(1) Dij21: %e Dij22: %e Dij23: %e\n",dij[1][0],dij[1][1],dij[1][2]);
    printf("(1) Dij31: %e Dij32: %e Dij33: %e\n",dij[2][0],dij[2][1],dij[2][2]);
#endif

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
      dij[0][0] = (c_length2*c_length2-2.0*arad*arad/3.0)/(4.0*M_PI*pow(c_length2,3)*visc);
      dij[1][1] = (c_length2*c_length2+2.0*arad*arad/3.0)/(8.0*M_PI*pow(c_length2,3)*visc);
      dij[2][2] = dij[1][1];
    } else {
      dij[0][0] = (1.0-3.0*c_length2/(16.0*arad))/(6.0*M_PI*arad*visc);
      dij[1][1] = (1.0-9.0*c_length2/(32.0*arad))/(6.0*M_PI*arad*visc);
      dij[2][2] = dij[1][1];
    }
#if 0
    printf("\n");
    printf("(2) Dii11: %e Dii12: %e Dii13: %e\n",dii[0][0],dii[0][1],dii[0][2]);
    printf("(2) Dii21: %e Dii22: %e Dii23: %e\n",dii[1][0],dii[1][1],dii[1][2]);
    printf("(2) Dii31: %e Dii32: %e Dii33: %e\n",dii[2][0],dii[2][1],dii[2][2]);
    printf("\n");
    printf("(2) Dij11: %e Dij12: %e Dij13: %e\n",dij[0][0],dij[0][1],dij[0][2]);
    printf("(2) Dij21: %e Dij22: %e Dij23: %e\n",dij[1][0],dij[1][1],dij[1][2]);
    printf("(2) Dij31: %e Dij32: %e Dij33: %e\n",dij[2][0],dij[2][1],dij[2][2]);
#endif
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
#else
    vel1[0] = 0.5*(f1_s1[0]+f1_s2[0]);
    vel1[1] = 0.5*(f1_s1[1]+f1_s2[1]);
    vel1[2] = 0.5*(f1_s1[2]+f1_s2[2]);
    vel2[0] = 0.5*(f2_s1[0]+f2_s2[0]);
    vel2[1] = 0.5*(f2_s1[1]+f2_s2[1]);
    vel2[2] = 0.5*(f2_s1[2]+f2_s2[2]);

    // remove force on center of mass
    for (int ii=0; ii<3; ii++) {
      f1_s1[ii] -= 0.5*vel1[ii];
      f1_s2[ii] -= 0.5*vel1[ii];
      f2_s1[ii] -= 0.5*vel2[ii];
      f2_s2[ii] -= 0.5*vel2[ii];
    }
    rotv1[0] = 0.5*(ry_i*f1_s1[2]-rz_i*f1_s1[1]);
    rotv1[1] = 0.5*(rz_i*f1_s1[0]-rx_i*f1_s1[2]);
    rotv1[2] = 0.5*(rx_i*f1_s1[1]-ry_i*f1_s1[0]);
    rotv1[0] += -0.5*(ry_i*f1_s2[2]-rz_i*f1_s2[1]);
    rotv1[1] += -0.5*(rz_i*f1_s2[0]-rx_i*f1_s2[2]);
    rotv1[2] += -0.5*(rx_i*f1_s2[1]-ry_i*f1_s2[0]);

    rotv2[0] = 0.5*(ry_j*f2_s1[2]-rz_j*f2_s1[1]);
    rotv2[1] = 0.5*(rz_j*f2_s1[0]-rx_j*f2_s1[2]);
    rotv2[2] = 0.5*(rx_j*f2_s1[1]-ry_j*f2_s1[0]);
    rotv2[0] += -0.5*(ry_j*f2_s2[2]-rz_j*f2_s2[1]);
    rotv2[1] += -0.5*(rz_j*f2_s2[0]-rx_j*f2_s2[2]);
    rotv2[2] += -0.5*(rx_j*f2_s2[1]-ry_j*f2_s2[0]);
#endif
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
 */
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void 
evaluateSurfaceForce(const Real *pos, Real *par, int *ipar, Real *vel, Real *rotv, Real *fpar)
{
  if (ipar[intIdx::cell_type] == cellType::YEAST) {
    Real z_wall = fpar[2];
    Real z_bndry_width = fpar[3];
    Real z_stiffness = fpar[4];
    Real z_gravity = fpar[5];
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
  } else if (ipar[intIdx::cell_type] == cellType::FUNGI) {
    Real z_wall = fpar[2];
    Real z_bndry_width = fpar[3];
    Real z_stiffness = fpar[4];
    Real z_gravity = fpar[5];
    Real visc = fpar[8];
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
fixBonds(int *ipar, int *jpar)
{
  int i_id[4];
  int i_cpu[4];
  int j_id[4];
  int j_cpu[4];
  int i_site[4];
  int j_site[4];
  int i_bnds = ipar[intIdx::n_bnds];
  i_id[0] = ipar[intIdx::seg1_id1];
  i_id[1] = ipar[intIdx::seg2_id1];
  i_id[2] = ipar[intIdx::seg3_id1];
  i_id[3] = ipar[intIdx::seg4_id1];

  i_cpu[0] = ipar[intIdx::seg1_id2];
  i_cpu[1] = ipar[intIdx::seg2_id2];
  i_cpu[2] = ipar[intIdx::seg3_id2];
  i_cpu[3] = ipar[intIdx::seg4_id2];

  i_site[0] = ipar[intIdx::site1];
  i_site[1] = ipar[intIdx::site2];
  i_site[2] = ipar[intIdx::site3];
  i_site[3] = ipar[intIdx::site4];


  int j_bnds = jpar[intIdx::n_bnds];
  j_id[0] = jpar[intIdx::seg1_id1];
  j_id[1] = jpar[intIdx::seg2_id1];
  j_id[2] = jpar[intIdx::seg3_id1];
  j_id[3] = jpar[intIdx::seg4_id1];

  j_cpu[0] = jpar[intIdx::seg1_id2];
  j_cpu[1] = jpar[intIdx::seg2_id2];
  j_cpu[2] = jpar[intIdx::seg3_id2];
  j_cpu[3] = jpar[intIdx::seg4_id2];

  j_site[0] = jpar[intIdx::site1];
  j_site[1] = jpar[intIdx::site2];
  j_site[2] = jpar[intIdx::site3];
  j_site[3] = jpar[intIdx::site4];

  int iid = ipar[intIdx::id];
  int icpu = ipar[intIdx::cpu];
  int jid = jpar[intIdx::id];
  int jcpu = jpar[intIdx::cpu];

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
    if (i_bnds == 1) {
      ipar[intIdx::seg2_id1] = jid;
      ipar[intIdx::seg2_id2] = jcpu;
      if (j_site[jidx] == 1) {
        ipar[intIdx::site2] = 2;
      } else {
        ipar[intIdx::site2] = 1;
      }
      i_bnds++;
      ipar[intIdx::n_bnds] = i_bnds;
    } else if (i_bnds == 2) {
      ipar[intIdx::seg3_id1] = jid;
      ipar[intIdx::seg3_id2] = jcpu;
      if (j_site[jidx] == 1) {
        ipar[intIdx::site3] = 2;
      } else {
        ipar[intIdx::site3] = 1;
      }
      i_bnds++;
      ipar[intIdx::n_bnds] = i_bnds;
    } else if (i_bnds == 3) {
      ipar[intIdx::seg4_id1] = jid;
      ipar[intIdx::seg4_id2] = jcpu;
      if (j_site[jidx] == 1) {
        ipar[intIdx::site4] = 2;
      } else {
        ipar[intIdx::site4] = 1;
      }
      i_bnds++;
      ipar[intIdx::n_bnds] = i_bnds;
    }
    printf("FIXING up i id: %d cpu: %d attached to id: %d cpu: %d\n",
        iid,icpu,jid,jcpu);
  } else if (i_bonded_to_j && !j_bonded_to_i) {
    // i has registered j as a bonding partner but j
    // has not registered i as a bonding partner. Fix up j
    if (j_bnds == 1) {
      jpar[intIdx::seg2_id1] = iid;
      jpar[intIdx::seg2_id2] = icpu;
      if (j_site[iidx] == 1) {
        jpar[intIdx::site2] = 2;
      } else {
        jpar[intIdx::site2] = 1;
      }
      j_bnds++;
      jpar[intIdx::n_bnds] = j_bnds;
    } else if (j_bnds == 2) {
      jpar[intIdx::seg3_id1] = iid;
      jpar[intIdx::seg3_id2] = icpu;
      if (j_site[iidx] == 1) {
        jpar[intIdx::site3] = 2;
      } else {
        jpar[intIdx::site3] = 1;
      }
      j_bnds++;
      jpar[intIdx::n_bnds] = j_bnds;
    } else if (j_bnds == 3) {
      jpar[intIdx::seg4_id1] = iid;
      jpar[intIdx::seg4_id2] = icpu;
      if (j_site[iidx] == 1) {
        jpar[intIdx::site4] = 2;
      } else {
        jpar[intIdx::site4] = 1;
      }
      j_bnds++;
      jpar[intIdx::n_bnds] = j_bnds;
    }
    printf("FIXING up j id: %d cpu: %d attached to id: %d cpu: %d\n",
        jid,jcpu,iid,icpu);
  }
}
#endif
