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
  Real r1 = p1[realIdx::radius];
  Real r2 = p2[realIdx::radius];
  return (r1+r2+bndry_width);
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
              Real *vel1, Real *vel2, Real *rotv1, Real *rotv2, Real *fpar )
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
    Real bond_strength = fpar[6];
    Real bond_cutoff = fpar[7];
    Real f1_s1[3];
    Real f1_s2[3];
    Real f2_s1[3];
    Real f2_s2[3];
    int i_id[3];
    int i_cpu[3];
    int i_site[3];
    int n_bnds = ipar1[intIdx::n_bnds];
    i_id[0] = ipar1[intIdx::seg1_id1];
    i_id[1] = ipar1[intIdx::seg2_id1];
    i_id[2] = ipar1[intIdx::seg3_id1];
    i_cpu[0] = ipar1[intIdx::seg1_id2];
    i_cpu[1] = ipar1[intIdx::seg2_id2];
    i_cpu[2] = ipar1[intIdx::seg3_id2];
    i_site[0] = ipar1[intIdx::site1];
    i_site[1] = ipar1[intIdx::site2];
    i_site[2] = ipar1[intIdx::site3];
    int j_id = ipar2[intIdx::id];
    int j_cpu = ipar2[intIdx::cpu];
    int i, j, k;
    for (i=0; i<3; i++) {
      f1_s1[i] = 0.0;
      f1_s2[i] = 0.0;
      f2_s1[i] = 0.0;
      f2_s2[i] = 0.0;
    }
    // find locations of segment end points
    Real xi1, yi1, zi1, xi2, yi2, zi2;
    Real xj1, yj1, zj1, xj2, yj2, zj2;
    Real c_length = par1[realIdx::c_length];
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
    xi1 = 0.5*nx*c_length;
    yi1 = 0.5*ny*c_length;
    zi1 = 0.5*nz*c_length;
    xi2 = -0.5*nx*c_length;
    yi2 = -0.5*ny*c_length;
    zi2 = -0.5*nz*c_length;
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

    c_length = par2[realIdx::c_length];
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
    xj1 = r12[0] + 0.5*nx*c_length;
    yj1 = r12[1] + 0.5*ny*c_length;
    zj1 = r12[2] + 0.5*nz*c_length;
    xj2 = r12[0] - 0.5*nx*c_length;
    yj2 = r12[1] - 0.5*ny*c_length;
    zj2 = r12[2] - 0.5*nz*c_length;

    // Only two possibilities for bonded interactions
    // since only sites 1 and 2 form bonds. Ignore possible 1-1 and 2-2
    // interactions.  Start with 1-2 interaction
    Real rx = xi1 - xj2;
    Real ry = yi1 - yj2;
    Real rz = zi1 - zj2;
    Real r = sqrt(rx*rx+ry*ry+rz*rz);
    if (r<bond_cutoff) {
      // Might be bonded. Check further
      for (i=0; i<n_bnds; i++) {
	if (i_site[i] == 1 && i_id[i] == j_id && i_cpu[i] == j_cpu) {
          f1_s1[0] += -bond_strength*rx;
          f1_s1[1] += -bond_strength*ry;
          f1_s1[2] += -bond_strength*rz;
          f2_s2[0] += bond_strength*rx;
          f2_s2[1] += bond_strength*ry;
          f2_s2[2] += bond_strength*rz;
	  break;
	}
      }
    }
    // 2-1 interaction
    rx = xi2 - xj1;
    ry = yi2 - yj1;
    rz = zi2 - zj1;
    r = sqrt(rx*rx+ry*ry+rz*rz);
    if (r<bond_cutoff) {
      // Might be bonded. Check further
      for (i=0; i<n_bnds; i++) {
	if (i_site[i] == 2 && i_id[i] == j_id && i_cpu[i] == j_cpu) {
          f1_s2[0] += -bond_strength*rx;
          f1_s2[1] += -bond_strength*ry;
          f1_s2[2] += -bond_strength*rz;
          f2_s1[0] += bond_strength*rx;
          f2_s1[1] += bond_strength*ry;
          f2_s1[2] += bond_strength*rz;
	  break;
	}
      }
    }
    // Evaluate interactions between cylinders
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
    for (i=0; i<3; i++) rirj += rho_i[i]*rho_j[i];
    Real tau_i, tau_j;
    // Find tau_i and tau_j corresponding to minimum distance between
    // cylinder segments
    if (rirj != 0.0) {
      // Find values of tau_i and tau_j that minimize
      // separation between lines defining cylinder
      // segments
      r = -(xi1-xj1)*rho_j[0];
      r -= (yi1-yj1)*rho_j[1];
      r -= (zi1-zj1)*rho_j[2];
      tau_i = -r/rirj;
      if (tau_i < 0.0) tau_i = 0.0;
      if (tau_i > 1.0) tau_i = 1.0;
      r = (xi1-xj1)*rho_i[0];
      r += (yi1-yj1)*rho_i[1];
      r += (zi1-zj1)*rho_i[2];
      tau_j = -r/rirj;
      if (tau_j < 0.0) tau_j = 0.0;
      if (tau_j > 1.0) tau_j = 1.0;
    } else {
      // This is an unlikely corner case
      Real rmin;
      rx = xi1-xj1;
      ry = yi1-yj1;
      rz = zi1-zj1;
      rmin = sqrt(rx*rx+ry*ry+rz*rz);
      tau_i = 0.0;
      tau_j = 0.0;

      rx = xi2-xj1;
      ry = yi2-yj1;
      rz = zi2-zj1;
      r = sqrt(rx*rx+ry*ry+rz*rz);
      if (r < rmin) {
	rmin = r;
	tau_i = 1.0;
	tau_j = 0.0;
      }

      rx = xi1-xj2;
      ry = yi1-yj2;
      rz = zi1-zj2;
      r = sqrt(rx*rx+ry*ry+rz*rz);
      if (r < rmin) {
	rmin = r;
	tau_i = 0.0;
	tau_j = 1.0;
      }

      rx = xi2-xj2;
      ry = yi2-yj2;
      rz = zi2-zj2;
      r = sqrt(rx*rx+ry*ry+rz*rz);
      if (r < rmin) {
	rmin = r;
	tau_i = 1.0;
	tau_j = 1.0;
      }
    }
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
    rx /= r;
    ry /= r;
    rz /= r;

    Real F;
    if (r < rA) {
      F= stiffness*(r-rA)*(r-rA)*(r-rS);
    } else {
      F = 0.0;
    }
    f1_s1[0] = F*rx*(1.0-tau_i);
    f1_s1[1] = F*ry*(1.0-tau_i);
    f1_s1[2] = F*rz*(1.0-tau_i);
    f1_s2[0] = F*rx*tau_i;
    f1_s2[1] = F*ry*tau_i;
    f1_s2[2] = F*rz*tau_i;
    f2_s1[0] = -F*rx*(1.0-tau_j);
    f2_s1[1] = -F*ry*(1.0-tau_j);
    f2_s1[2] = -F*rz*(1.0-tau_j);
    f2_s2[0] = -F*rx*tau_j;
    f2_s2[1] = -F*ry*tau_j;
    f2_s2[2] = -F*rz*tau_j;

    // Subtract out the part of force that is moving sites together
    // and evaluate velocity of CM and rotational velocity
    Real rx_i = xi1-xi2;
    Real ry_i = yi1-yi2;
    Real rz_i = zi1-zi2;
    r = sqrt(rx_i*rx_i+ry_i*ry_i+rz_i*rz_i);
    if (r < 1.0e-10) {
      printf("rx_i: %f ry_i: %f rz_i: %f r: %f\n",rx_i,ry_i,rz_i,r);
      printf("xi1: %f yi1: %f zi1: %f xi2: %f yi2: %f zi2: %f\n",xi1,yi1,zi1,xi2,yi2,zi2);
    }
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
    if (r < 1.0e-10) 
      printf("rx_j: %f ry_j: %f rz_j: %f r: %f\n",rx_j,ry_j,rz_j,r);
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
    for (i=0; i<3; i++) {
      for (j=0; j<3; j++) {
	if (i==j) {
	  dii[i][j] = 1.0/(6.0*M_PI*arad);
	} else {
	  dii[i][j] = 0.0;
	}
	dij[i][j] = 0.0;
      }
    }
    if (c_length > 2.0*arad) {
      dij[0][0] = (c_length*c_length+2.0*arad*arad/3.0)/(8.0*M_PI*pow(c_length,3));
      dij[1][1] = dij[0][0];
      dij[2][2] = (c_length*c_length-2.0*arad*arad/3.0)/(4.0*M_PI*pow(c_length,3));
    } else {
      dij[0][0] = (1.0-9.0*c_length/(32.0*arad))/(6.0*arad);
      dij[1][1] = dij[0][0];
      dij[2][2] = (1.0-3.0*c_length/(16.0*arad))/(6.0*arad);
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
	  dtmp[i][j] += dij[i][k]*rot_i[j][k];
	}
      }
    }
    for (i=0; i<3; i++) {
      for (j=0; j<3; j++) {
	dij[i][j] = 0.0;
	for (k=0; k<3; k++) {
	  dij[i][j] += rot_i[i][k]*dtmp[k][j];
	}
      }
    }
    // calculate velocities on sites 1 and 2
    Real vi_1[3], vi_2[3];
    for (i=0; i<3; i++) {
      vi_1[i] = 0.0;
      vi_2[i] = 0.0;
      for (j=0; j<3; j++) {
	vi_1[i] += dii[i][j]*f1_s1[j]+dij[i][j]*f1_s2[j];
	vi_2[i] += dii[i][j]*f1_s2[j]+dij[i][j]*f1_s1[j];
      }
    } 

    // Calculate friction tensor for segment j. Start by calculating
    // size of spheres with same exposed surface as cylinder
    area = par2[realIdx::area];
    area /= 2.0;
    arad = sqrt(area/(4.0*M_PI));
    // Assume the cylinder is oriented along the x-axis
    for (i=0; i<3; i++) {
      for (j=0; j<3; j++) {
	if (i==j) {
	  dii[i][j] = 1.0/(6.0*M_PI*arad);
	} else {
	  dii[i][j] = 0.0;
	}
	dij[i][j] = 0.0;
      }
    }
    if (c_length > 2.0*arad) {
      dij[0][0] = (c_length*c_length+2.0*arad*arad/3.0)/(8.0*M_PI*pow(c_length,3));
      dij[1][1] = dij[0][0];
      dij[2][2] = (c_length*c_length-2.0*arad*arad/3.0)/(4.0*M_PI*pow(c_length,3));
    } else {
      dij[0][0] = (1.0-9.0*c_length/(32.0*arad))/(6.0*arad);
      dij[1][1] = dij[0][0];
      dij[2][2] = (1.0-3.0*c_length/(16.0*arad))/(6.0*arad);
    }
    // Reorient dij so that it reflects the correct direction (dii is diagonal
    // so it doesn't change under reorientation). Note the transpose of the
    // rotation is equal to its inverse
    for (i=0; i<3; i++) {
      for (j=0; j<3; j++) {
	dtmp[i][j] = 0.0;
	for (k=0; k<3; k++) {
	  // use inverse rotation here
	  dtmp[i][j] += dij[i][k]*rot_j[j][k];
	}
      }
    }
    for (i=0; i<3; i++) {
      for (j=0; j<3; j++) {
	dij[i][j] = 0.0;
	for (k=0; k<3; k++) {
	  dij[i][j] += rot_j[i][k]*dtmp[k][j];
	}
      }
    }
    // calculate velocities on sites 1 and 2
    Real vj_1[3], vj_2[3];
    for (i=0; i<3; i++) {
      vj_1[i] = 0.0;
      vj_2[i] = 0.0;
      for (j=0; j<3; j++) {
	vj_1[i] += dii[i][j]*f2_s1[j]+dij[i][j]*f2_s2[j];
	vj_2[i] += dii[i][j]*f2_s2[j]+dij[i][j]*f2_s1[j];
      }
    } 
    // Get translational and rotational velocities of each
    // segment
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
    Real f_s1[3];
    Real f_s2[3];
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
    if (pos[2] < 0.00478) {
      printf("POS_Z: %e NX: %f NY: %f NZ: %f THETA: %f PHI: %f\n",pos[2],nx,ny,nz,theta,phi);
    }
    // Find interaction centers on segment
    x1 = pos[0] + 0.5*nx*c_length;
    y1 = pos[1] + 0.5*ny*c_length;
    z1 = pos[2] + 0.5*nz*c_length;
    x2 = pos[0] - 0.5*nx*c_length;
    y2 = pos[1] - 0.5*ny*c_length;
    z2 = pos[2] - 0.5*nz*c_length;
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
    if (dz > 0.0 && dz < rA) {
      f_s1[2] = -z_stiffness*(dz-rA)*(dz-rA)*(dz-rS);
    } else if (dz <= 0.0) {
      f_s1[2] = z_stiffness*rA*rA*rS;
    } else {
      f_s1[2] = 0.0;
    }
    f_s1[2] -= z_gravity*z_stiffness*rA*rA*rS;
    if (pos[2] < 0.00478) {
      printf("POS_Z1: %e DZ: %f fz1: %e\n",z1,dz,f_s1[2]);
    }
    z = z2;
    dz = z - z_wall; 
    rS = par[realIdx::radius];
    rA = rS + z_bndry_width;
    if (dz > 0.0 && dz < rA) {
      f_s2[2] = -z_stiffness*(dz-rA)*(dz-rA)*(dz-rS);
    } else if (dz <= 0.0) {
      f_s2[2] = z_stiffness*rA*rA*rS;
    } else {
      f_s2[2] = 0.0;
    }
    f_s2[2] -= z_gravity*z_stiffness*rA*rA*rS;
    if (pos[2] < 0.00478) {
      printf("POS_Z2: %e DZ: %f fz2: %e\n",z2,dz,f_s2[2]);
    }

    // Subtract out the part of force that is moving sites together
    // and evaluate velocity of CM and rotational velocity
    Real rx = x1-x2;
    Real ry = y1-y2;
    Real rz = z1-z2;
    Real r = sqrt(rx*rx+ry*ry+rz*rz);
    if (r == 0) {
      printf("Separation is zero\n");
    }
    nx = rx/r;
    ny = ry/r;
    nz = rz/r;
    Real fx = f_s1[0]-f_s2[0];
    Real fy = f_s1[1]-f_s2[1];
    Real fz = f_s1[2]-f_s2[2];
    Real fn = fx*nx + fy*ny + fz*nz;
    f_s1[0] -= 0.5*fn*nx;
    f_s1[1] -= 0.5*fn*ny;
    f_s1[2] -= 0.5*fn*nz;
    f_s2[0] += 0.5*fn*nx;
    f_s2[1] += 0.5*fn*ny;
    f_s2[2] += 0.5*fn*nz;
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
	  dii[i][j] = 1.0/(6.0*M_PI*arad);
	} else {
	  dii[i][j] = 0.0;
	}
	dij[i][j] = 0.0;
      }
    }
    if (c_length > 2.0*arad) {
      dij[0][0] = (c_length*c_length-2.0*arad*arad/3.0)/(4.0*M_PI*pow(c_length,3));
      dij[1][1] = (c_length*c_length+2.0*arad*arad/3.0)/(8.0*M_PI*pow(c_length,3));
      dij[2][2] = dij[1][1];
    } else {
      dij[0][0] = (1.0-3.0*c_length/(16.0*arad))/(6.0*M_PI*arad);
      dij[1][1] = (1.0-9.0*c_length/(32.0*arad))/(6.0*M_PI*arad);
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
    // Get translational and rotational velocities for segment
    vel[0] = 0.5*(v_1[0]+v_2[0]);
    vel[1] = 0.5*(v_1[1]+v_2[1]);
    vel[2] = 0.5*(v_1[2]+v_2[2]);
    if (pos[2] < 0.00478) {
      printf("VELOCITY_Z: %e Z: %f fz1: %e fz2: %e\n",vel[2],pos[2],f_s1[2],f_s2[2]);
    }

    rotv[0] = 0.5*(ry*v_1[2]-rz*v_1[1]);
    rotv[1] = 0.5*(rz*v_1[0]-rx*v_1[2]);
    rotv[2] = 0.5*(rx*v_1[1]-ry*v_1[0]);
    rotv[0] += -0.5*(ry*v_2[2]-rz*v_2[1]);
    rotv[1] += -0.5*(rz*v_2[0]-rx*v_2[2]);
    rotv[2] += -0.5*(rx*v_2[1]-ry*v_2[0]);
  }
}
