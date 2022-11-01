/**
 * Check to see if particle meets criteria for splitting
 * @param p_par real particle parameters
 * @param p_ipar integer particle parameters
 * @param max_vol maximum volume for YEAST splitting
 * @param max_len maximum segment length for splitting of FUNGI
 * @param max_rad maximum radius of FUNGI segments
 * @param brnch_prob probability per step that a side branch will be created
 * @param engine object for implementing random number generator
 * @return true if particle should split
 */
AMREX_GPU_HOST_DEVICE AMREX_INLINE
bool 
checkSplit(Real *p_par, int *p_ipar, Real max_vol, Real max_len,
    Real max_rad, Real brnch_prob, amrex::RandomEngine const& engine)
{
  bool ret = false;
  int cell_type = p_ipar[intIdx::cell_type];
  if (cell_type == cellType::YEAST) {
    if (p_par[realIdx::vol] > max_vol) ret = true;
  } else if (cell_type == cellType::FUNGI) {
    //TODO: Check fungi values
    int position = p_ipar[intIdx::position];
    if (position == siteLocation::TIP) {
      // Only consider splitting if segment cannot grow radially
      if (p_par[realIdx::radius] >= max_rad) {
        // If segment is a growth tip then split if segment is greater than
        // or equal to maximum length
        if (p_par[realIdx::c_length] >= max_len) {
          printf("LENGTH: %16.8f\n",p_par[realIdx::c_length]);
          ret = true;
        }
      }
    } else if (position == siteLocation::SECOND_1 ||
               position == siteLocation::SECOND_2) {
      // Segment is immediately behind the growth tip. Generate side
      // branch with some probability.
      Real x = amrex::Random(engine);
      if (x < brnch_prob) ret = true;
      printf("CHECKING FOR SIDE BRANCH brnch_prob: %f x: %f\n",brnch_prob,x);
      if (ret) printf("SIDE BRANCH WILL BE CREATED\n");
    }
  }
  return ret;
}

/**
 * Duplicate data from original particle to child when splitting
 * @param p_real_orig pointer to real values from original particle
 * @param p_int_orig pointer to integer values from original particle
 * @param p_real_child pointer to real values on child  particle
 * @param p_int_child pointer to integer values on child particle
 */
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void 
copyChildParameters(Real *p_real_orig, int *p_int_orig,
                   Real *p_real_child, int *p_int_child,
                   int num_reals, int num_ints) 
{
  int i;
  int nreals = num_reals + realIdx::count-1;
  int nints  = num_ints  +  intIdx::count-1;
  // Copy all values to child, fix up details below
  for (i=0; i<nreals; i++) p_real_child[i] = p_real_orig[i];
  for (i=0; i<nints; i++) p_int_child[i] = p_int_orig[i];
}

/**
 * Create positions and parameter values for new particle pair. This consist
 * of the original particle and a new particle
 * @param pos_orig, pos_new positions of original and new particles
 * @param par_orig, par_new real parameter values of original and new particles
 * @param ipar_orig, ipar_new integer parameter values of original and new particles
 * @param overlap degree of overlap of new and old cell
 * @param num_reals, num_ints number of real and integer parameters
 * @param engine object for implementing random number generator
 */
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void 
setNewCell(Real *pos_orig, Real *pos_new, Real *par_orig,
           Real *par_new, int *ipar_orig, int *ipar_new,
           Real overlap, int num_reals, int num_ints,
           amrex::RandomEngine const& engine)
{
  // Copy values from original particle to child
  copyChildParameters(par_orig, ipar_orig, par_new, ipar_new,
      num_reals, num_ints);

  if (ipar_orig[intIdx::cell_type] == cellType::YEAST) {
    Real x,y,z;
    x = pos_orig[0];
    y = pos_orig[1];
    z = pos_orig[2];

    // Find new locations for split particles
    Real radius = par_new[realIdx::radius];
    Real theta = 2.0 * amrex::Random(engine) - 1.0;
    theta = acos(theta);
    Real phi = 2.0 * M_PI * amrex::Random(engine);
    Real nx = sin(theta)*cos(phi);
    Real ny = sin(theta)*sin(phi);
    Real nz = cos(theta);
    Real scale = 1.0 - overlap;
    pos_new[0] = x + scale*0.5*nx*radius;
    pos_new[1] = y + scale*0.5*ny*radius;
    pos_new[2] = z + scale*0.5*nz*radius;
    pos_orig[0] = x - scale*0.5*nx*radius;
    pos_orig[1] = y - scale*0.5*ny*radius;
    pos_orig[2] = z - scale*0.5*nz*radius;

    // Fix up some parameters to reflect split
    Real volume = par_orig[realIdx::vol]/2.0;
    radius = pow((3.0*volume/(4.0*M_PI)),1.0/3.0);
    Real area = 4.0*M_PI*radius*radius;
    Real dvdt = par_orig[realIdx::dvdt];
    Real dadt = 2.0*dvdt/radius;
    par_orig[realIdx::vol] = volume;
    par_new[realIdx::vol] = volume;
    par_orig[realIdx::area] = area;
    par_new[realIdx::area] = area;
    par_orig[realIdx::dadt] = dadt;
    par_new[realIdx::dadt] = dadt;
    par_orig[realIdx::radius] = radius;
    par_new[realIdx::radius] = radius;
  } else if (ipar_orig[intIdx::cell_type] == cellType::FUNGI) {
    // Fungi split handled in setNewSegment
  }
}

/**
 * Create positions and parameter values for new segment pair. This consist
 * of the original segment and a new segment
 * @param pos_orig, pos_new positions of original and new segments
 * @param par_orig, par_new real parameter values of original and new segments
 * @param ipar_orig, ipar_new integer parameter values of original and new segments
 * @param num_reals, num_ints number of real and integer parameters
 * @param split_len length of parent segment after split
 * @param max_rad  maximum radius of parent segment
 * @param id, cpu ID and cpu of new segment
 * @param engine object for implementing random number generator
 */
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void 
setNewSegment(Real *pos_orig, Real *pos_new, Real *par_orig,
              Real *par_new, int *ipar_orig, int *ipar_new, int num_reals,
              int num_ints, Real split_len, Real /*max_rad*/, int id, int cpu,
              amrex::RandomEngine const& engine)
{
  // Copy values from original particle to child
  copyChildParameters(par_orig, ipar_orig, par_new, ipar_new,
      num_reals, num_ints);

  // Get chemical gradient
  Real gx = par_orig[realIdx::gx];
  Real gy = par_orig[realIdx::gy];
  Real gz = par_orig[realIdx::gz];
  if (true) {
    gx = 0.0;
    gy = 0.0;
    gz = 0.0;
  }

  //Check out if this is an interior or terminal segment
  int nbnd = ipar_orig[intIdx::n_bnds];
  int position = ipar_orig[intIdx::position];
  if (nbnd == 0) {
    // Isolated segment
    Real x_orig = pos_orig[0];
    Real y_orig = pos_orig[1];
    Real z_orig = pos_orig[2];
    Real c_length = par_orig[realIdx::c_length];
    // Find orientation of segment
    Real theta = par_orig[realIdx::theta];
    Real phi = par_orig[realIdx::phi];
    Real ct = cos(theta);
    Real st = sin(theta);
    Real cp = cos(phi);
    Real sp = sin(phi);
    Real nx = st*cp;
    Real ny = st*sp;
    Real nz = ct;
    // Calculate location of new segment centers
    pos_orig[0] = x_orig + 0.25*c_length*nx;
    pos_orig[1] = y_orig + 0.25*c_length*ny;
    pos_orig[2] = z_orig + 0.25*c_length*nz;
    pos_new[0] = x_orig - 0.25*c_length*nx;
    pos_new[1] = y_orig - 0.25*c_length*ny;
    pos_new[2] = z_orig - 0.25*c_length*nz;
    par_orig[realIdx::c_length] = 0.5*c_length;
    par_new[realIdx::c_length] = 0.5*c_length;

    /* Fix up topology of new and old segment
              1-----------2 parent
                  /    \
                 /      \
      (child) 1----2  1----2 (parent) */

    ipar_orig[intIdx::n_bnds] = 1;
    ipar_orig[intIdx::seg1_id1] = id;
    ipar_orig[intIdx::seg1_id2] = cpu;
    ipar_orig[intIdx::site1] = 1;
    ipar_orig[intIdx::position] = siteLocation::TIP;

    ipar_new[intIdx::n_bnds] = 1;
    ipar_new[intIdx::seg1_id1] = ipar_orig[intIdx::id];
    ipar_new[intIdx::seg1_id2] = ipar_orig[intIdx::cpu];
    ipar_new[intIdx::site1] = 2;
    ipar_new[intIdx::position] = siteLocation::TIP;
    // Growth is towards site 1
    ipar_new[intIdx::id] = id;
    ipar_new[intIdx::cpu] = cpu;

  } else if (position == siteLocation::TIP) {
    // Terminal segment (growth tip)
    // Find direction of chemical gradient
    Real gn = sqrt(gx*gx+gy*gy+gz*gz);

    // Find out which end is bonded to a segment
    int bsite = ipar_orig[intIdx::site1];
    if (nbnd == 2) {
      // Check to make sure second bonding site is equal to first bonding
      // site (in case where TIP is generated by side branching)
      if (ipar_orig[intIdx::site2] != bsite) {
        amrex::Abort("ERROR: bonds at both ends of tip segment");
      }
    }
    //printf("NBNDS: %d SITE1: %d SITE2: %d SITE3: %d\n",
    //    ipar_orig[intIdx::n_bnds],ipar_orig[intIdx::site1],ipar_orig[intIdx::site2],
    //    ipar_orig[intIdx::site3]);

    // Find orientation of existing segment
    Real theta = par_orig[realIdx::theta];
    Real phi = par_orig[realIdx::phi];
    Real ct = cos(theta);
    Real st = sin(theta);
    Real cp = cos(phi);
    Real sp = sin(phi);
    Real nx = st*cp;
    Real ny = st*sp;
    Real nz = ct;
    if (gn == 0.0) {
      gx = nx;
      gy = ny;
      gz = nz;
      //      printf("THETA: %f PHI: %f GX: %e GY: %e GZ: %e GN: %e\n",theta,phi,gx,gy,gz,gn);
      gn = 1.0;
    }

    // create unit vector in direction of chemical gradient
    gx = gx/gn;
    gy = gy/gn;
    gz = gz/gn;

    //printf("THETA: %f PHI: %f GX: %e GY: %e GZ: %e GN: %e\n",theta,phi,gx,gy,gz,gn);
    Real gtheta = acos(gz);
    Real gst = sin(gtheta);
    Real gct = gz;
    Real gphi = 0.0;
    if (gst != 0.0) {
      gphi = acos(gx/gst);
      if (gy < 0.0) gphi = 2.0 * M_PI - gphi;
    }

    Real gcp = cos(gphi);
    Real gsp = sin(gphi);
    //printf("GTHETA: %f GPHI: %f GX: %e GY: %e GZ: %e GN: %e\n",theta,phi,gx,gy,gz,gn);
    // Find matrix that rotates x-axis to gradient direction
    Real rotx[3][3];
    rotx[0][0] = gcp*gst;
    rotx[0][1] = -gsp;
    rotx[0][2] = -gcp*gct;
    rotx[1][0] = gsp*gst;
    rotx[1][1] = gcp;
    rotx[1][2] = -gsp*gct;
    rotx[2][0] = gct;
    rotx[2][1] = 0.0;
    rotx[2][2] = gst;
    // Choose a direction for the new segment from a distribution of directions
    // around the x-axis
    Real rtheta = 0.1 * amrex::Random(engine);
    Real rphi = 2.0 * M_PI * amrex::Random(engine);
    Real rct = cos(rtheta);
    Real rst = sin(rtheta);
    Real rcp = cos(rphi);
    Real rsp = sin(rphi);
    Real rx = rct;
    Real ry = rst*rsp;
    Real rz = rst*rcp;
    // Rotate random vector so that it is approximately oriented in the direction
    // of chemical gradient
    Real gnx = rotx[0][0]*rx+rotx[0][1]*ry+rotx[0][2]*rz;
    Real gny = rotx[1][0]*rx+rotx[1][1]*ry+rotx[1][2]*rz;
    Real gnz = rotx[2][0]*rx+rotx[2][1]*ry+rotx[2][2]*rz;
    //  printf("RTHETA: %f RPHI: %f GNX: %f GNY: %f GNZ: %f NX: %f NY: %f NZ: %f\n",
    //      rtheta,rphi,gnx,gny,gnz,nx,ny,nz);
    // Now have orientation of old segment (nx,ny,nz) and new segment
    // (gnx,gny,gnz). Finish constructing new configuration.
    Real x = pos_orig[0];
    Real y = pos_orig[1];
    Real z = pos_orig[2];
    Real c_length = par_orig[realIdx::c_length];
    Real x1, y1, z1, x2, y2, z2;
    if (bsite == 1) {
    //  printf("Add to site 2\n");
      // new segment splits off site 2
      x1 = x - 0.5*c_length*nx;
      y1 = y - 0.5*c_length*ny;
      z1 = z - 0.5*c_length*nz;
      x2 = x1 + split_len*nx;
      y2 = y1 + split_len*ny;
      z2 = z1 + split_len*nz;
      pos_orig[0] = 0.5*(x1+x2);
      pos_orig[1] = 0.5*(y1+y2);
      pos_orig[2] = 0.5*(z1+z2);
      par_orig[realIdx::c_length] = split_len;
      if (nbnd == 1) {
        ipar_orig[intIdx::site2] = 2;
      } else {
        ipar_orig[intIdx::site3] = 2;
      }
      ipar_orig[intIdx::position] = siteLocation::SECOND_2;
      Real seg_len = c_length-split_len;
      // find location of Site 2 on new segment (site 1 on new segment equals
      // site 2 of old segment)
      Real ax = x2 + gnx*seg_len;
      Real ay = y2 + gny*seg_len;
      Real az = z2 + gnz*seg_len;
      //printf("B1 X1: %e Y1: %e Z1: %e X2: %e Y2: %e Z2: %e X3: %e Y3: %e Z3: %e\n",
      //    x1,y1,z1,x2,y2,z2,ax,ay,az);
      pos_new[0] = 0.5*(x2 + ax);
      pos_new[1] = 0.5*(y2 + ay);
      pos_new[2] = 0.5*(z2 + az);
      par_new[realIdx::c_length] = seg_len;
      gtheta = acos(gnz);
      gst = sin(gtheta);
      gphi = 0.0;
      if (gst != 0.0) {
        gphi = acos(gnx/gst);
        if (gny < 0.0) gphi = 2.0 * M_PI - gphi;
      }
      //printf("SITE: %d GTHETA: %f GPHI: %f\n",bsite,gtheta,gphi);
      par_new[realIdx::theta] = gtheta;
      par_new[realIdx::phi] = gphi;
      ipar_new[intIdx::site1] = 1;
      ipar_new[intIdx::position] = siteLocation::TIP;
    } else {
      //printf("Add to site 1\n");
      x2 = x + 0.5*c_length*nx;
      y2 = y + 0.5*c_length*ny;
      z2 = z + 0.5*c_length*nz;
      x1 = x2 - split_len*nx;
      y1 = y2 - split_len*ny;
      z1 = z2 - split_len*nz;
      pos_orig[0] = 0.5*(x1+x2);
      pos_orig[1] = 0.5*(y1+y2);
      pos_orig[2] = 0.5*(z1+z2);
      par_orig[realIdx::c_length] = split_len;
      if (nbnd == 1) {
        ipar_orig[intIdx::site2] = 1;
      } else {
        ipar_orig[intIdx::site3] = 1;
      }
      ipar_orig[intIdx::position] = siteLocation::SECOND_1;
      Real seg_len = c_length-split_len;
      // find location of Site 1 on new segment (site 2 on new segment equals
      // site 1 of old segment)
      Real ax = x1 - gnx*seg_len;
      Real ay = y1 - gny*seg_len;
      Real az = z1 - gnz*seg_len;
      //printf("B2 X1: %e Y1: %e Z1: %e X2: %e Y2: %e Z2: %e X3: %e Y3: %e Z3: %e\n",
      //    x1,y1,z1,x2,y2,z2,ax,ay,az);
      pos_new[0] = 0.5*(x1 + ax);
      pos_new[1] = 0.5*(y1 + ay);
      pos_new[2] = 0.5*(z1 + az);
      par_new[realIdx::c_length] = seg_len;
      gtheta = acos(gnz);
      gst = sin(gtheta);
      gphi = 0.0;
      if (gst != 0.0) {
        gphi = acos(gnx/gst);
        if (gny < 0.0) gphi = 2.0 * M_PI - gphi;
      }
      //printf("SITE: %d GTHETA: %f GPHI: %f\n",bsite,gtheta,gphi);
      par_new[realIdx::theta] = gtheta;
      par_new[realIdx::phi] = gphi;
      ipar_new[intIdx::site1] = 2;
      ipar_new[intIdx::position] = siteLocation::TIP;
    }
//    printf("PARENT LENGTH: %e\n",par_orig[realIdx::c_length]);
//    printf("CHILD LENGTH: %e\n",par_new[realIdx::c_length]);
//    printf("Adding new segment ID: %d CPU: %d\n",id,cpu);
    if (nbnd == 1) {
      ipar_orig[intIdx::n_bnds] = 2; 
      ipar_orig[intIdx::seg2_id1] = id;
      ipar_orig[intIdx::seg2_id2] = cpu;
    } else {
      ipar_orig[intIdx::n_bnds] = 3; 
      ipar_orig[intIdx::seg3_id1] = id;
      ipar_orig[intIdx::seg3_id2] = cpu;
    }
    ipar_new[intIdx::n_bnds] = 1; 
    ipar_new[intIdx::seg1_id1] = ipar_orig[intIdx::id];
    ipar_new[intIdx::seg1_id2] = ipar_orig[intIdx::cpu];
  } else if (ipar_orig[intIdx::position] == siteLocation::SECOND_1
        || ipar_orig[intIdx::position] == siteLocation::SECOND_2) {
    // Find direction of chemical gradient
    Real gn = sqrt(gx*gx+gy*gy+gz*gz);
    // Find orientation of existing segment
    Real theta = par_orig[realIdx::theta];
    Real phi = par_orig[realIdx::phi];
    Real ct = cos(theta);
    Real st = sin(theta);
    Real cp = cos(phi);
    Real sp = sin(phi);
    Real nx = st*cp;
    Real ny = st*sp;
    Real nz = ct;
    if (gn == 0.0) {
      // if gradient vanishes, then pick a direction parallel to the xy plane and
      // at a 45 degree angle to the growth direction
      Real da = amrex::Random(engine);
      if (da > 0.5) {
        da = 0.25*M_PI;
      } else {
        da = -0.25*M_PI;
      }
      Real tphi;
#if 0
      if (ipar_orig[intIdx::position] == siteLocation::SECOND_1) {
        printf("Attaching at site 1\n");
        // Growing tip is attached to site 1 so branch is oriented in same
        // general direction as segment
        tphi = phi+da;
      } else {
        printf("Attaching at site 2\n");
        // Growing tip is attached to site 2 so branch is oriented in
        // approximately opposite direction as segment
        tphi = phi+M_PI+da;
      }
#endif
      tphi = phi+da;
      gx = cos(tphi);
      gy = sin(tphi);
      gz = 0;
      gn = 1.0;
//      printf("GRADIENT nx: %f ny: %f nz: %f gx: %f gy: %f gz: %f\n",nx,ny,nz,gx,gy,gz);
    }
    // create unit vector in direction of chemical gradient
    gx = gx/gn;
    gy = gy/gn;
    gz = gz/gn;
    Real gtheta = acos(gz);
    Real gst = sin(gtheta);
    Real gct = gz;
    Real gphi = 0.0;
    if (gst != 0.0) {
      gphi = acos(gx/gst);
      if (gy < 0.0) gphi = 2.0 * M_PI - gphi;
    }
    Real gcp = cos(gphi);
    Real gsp = sin(gphi);
    // Find matrix that rotates x-axis to gradient direction
    Real rotx[3][3];
    rotx[0][0] = gcp*gst;
    rotx[0][1] = -gsp;
    rotx[0][2] = -gcp*gct;
    rotx[1][0] = gsp*gst;
    rotx[1][1] = gcp;
    rotx[1][2] = -gsp*gct;
    rotx[2][0] = gct;
    rotx[2][1] = 0.0;
    rotx[2][2] = gst;
    // Choose a direction for the new segment from a distribution of directions
    // around the x-axis
    Real rtheta = 0.1 * amrex::Random(engine);
    Real rphi = 2.0 * M_PI * amrex::Random(engine);
    Real rct = cos(rtheta);
    Real rst = sin(rtheta);
    Real rcp = cos(rphi);
    Real rsp = sin(rphi);
    Real rx = rct;
    Real ry = rst*rsp;
    Real rz = rst*rcp;
    // Rotate random vector so that it is approximately oriented in the direction
    // of chemical gradient
    Real gnx = rotx[0][0]*rx+rotx[0][1]*ry+rotx[0][2]*rz;
    Real gny = rotx[1][0]*rx+rotx[1][1]*ry+rotx[1][2]*rz;
    Real gnz = rotx[2][0]*rx+rotx[2][1]*ry+rotx[2][2]*rz;
//      printf("RTHETA: %f RPHI: %f GNX: %f GNY: %f GNZ: %f NX: %f NY: %f NZ: %f\n",
//          rtheta,rphi,gnx,gny,gnz,nx,ny,nz);
    // Find theta and phi for new segment
    Real gntheta = acos(gnz);
    Real gnst = sin(gntheta);
    // Real gnct = gz;
    Real gnphi = 0.0;
    if (gnst != 0.0) {
      gnphi = acos(gnx/gnst);
      if (gny < 0.0) gnphi = 2.0 * M_PI - gnphi;
    }
    // Now have orientation of old segment (nx,ny,nz) and new segment
    // (gnx,gny,gnz). Finish constructing new configuration.
    Real x = pos_orig[0];
    Real y = pos_orig[1];
    Real z = pos_orig[2];
    Real c_length = par_orig[realIdx::c_length];
    Real x1, y1, z1, x2, y2, z2;
    // Create a branch
    x1 = x - 0.5*c_length*nx;
    y1 = y - 0.5*c_length*ny;
    z1 = z - 0.5*c_length*nz;
    x2 = x + 0.5*c_length*nx;
    y2 = y + 0.5*c_length*ny;
    z2 = z + 0.5*c_length*nz;
    Real seg_len = 0.0005;
    if (ipar_orig[intIdx::position] == siteLocation::SECOND_1) {
      // Growing tip is attached to site 1 so side branching will be
      // from site 1
      if (nbnd == 2) {
        ipar_orig[intIdx::site3] = 1;
      } else if (nbnd == 3) {
        ipar_orig[intIdx::site4] = 1;
      }
      // find location of Site 1 on new segment (site 2 on new segment equals
      // site 1 of old segment)
      Real ax = x1 - gnx*seg_len;
      Real ay = y1 - gny*seg_len;
      Real az = z1 - gnz*seg_len;
      pos_new[0] = 0.5*(x1 + ax);
      pos_new[1] = 0.5*(y1 + ay);
      pos_new[2] = 0.5*(z1 + az);
      ipar_new[intIdx::site1] = 2;
      ipar_new[intIdx::site2] = 2;
    } else {
      // Growing tip is attached to site 2 so side branching will be
      // from site 2
      if (nbnd == 2) {
        ipar_orig[intIdx::site3] = 2;
      } else if (nbnd == 3) {
        ipar_orig[intIdx::site4] = 2;
      }
      // find location of Site 2 on new segment (site 1 on new segment equals
      // site 2 of old segment)
      Real ax = x2 + gnx*seg_len;
      Real ay = y2 + gny*seg_len;
      Real az = z2 + gnz*seg_len;
      pos_new[0] = 0.5*(x2 + ax);
      pos_new[1] = 0.5*(y2 + ay);
      pos_new[2] = 0.5*(z2 + az);
      ipar_new[intIdx::site1] = 1;
      ipar_new[intIdx::site2] = 1;
    }
    // Volume for new segment comes from shrinking radius of old segment
    Real radius = par_orig[realIdx::radius];
//    printf("OLD RADIUS: %e\n",radius);
    Real old_vol =  c_length*radius*radius*M_PI;
    par_new[realIdx::theta] = gntheta;
    par_new[realIdx::phi] = gnphi;
    // New radius is calculated so volume of old segment and new segment equals
    // old volume
    radius = old_vol/(M_PI*(c_length+par_new[realIdx::c_length]));
    radius = sqrt(radius);
//    printf("NEW RADIUS: %e\n",radius);
    par_orig[realIdx::radius] = radius;
    par_new[realIdx::radius] = radius;
    par_new[realIdx::c_length] = seg_len;
    if (nbnd == 2) {
      ipar_orig[intIdx::n_bnds] = 3; 
      ipar_orig[intIdx::seg3_id1] = id;
      ipar_orig[intIdx::seg3_id2] = cpu;
    } else {
      ipar_orig[intIdx::n_bnds] = 4; 
      ipar_orig[intIdx::seg4_id1] = id;
      ipar_orig[intIdx::seg4_id2] = cpu;
    }
    ipar_new[intIdx::n_bnds] = 2; 
    ipar_new[intIdx::seg1_id1] = ipar_orig[intIdx::id];
    ipar_new[intIdx::seg1_id2] = ipar_orig[intIdx::cpu];
    if (nbnd == 2) {
      ipar_new[intIdx::seg2_id1] = ipar_orig[intIdx::seg2_id1];
      ipar_new[intIdx::seg2_id2] = ipar_orig[intIdx::seg2_id2];
    } else {
      ipar_new[intIdx::seg2_id1] = ipar_orig[intIdx::seg3_id1];
      ipar_new[intIdx::seg2_id2] = ipar_orig[intIdx::seg3_id2];
    }
    ipar_new[intIdx::position] = siteLocation::TIP;
    ipar_orig[intIdx::position] = siteLocation::INTERIOR;
  }
  // Recalculate area and volume
  Real clen = par_orig[realIdx::c_length];
  Real rad = par_orig[realIdx::radius];
  par_orig[realIdx::vol] = M_PI*rad*rad*clen;
  par_orig[realIdx::area] = 2.0*M_PI*rad*clen+2.0*M_PI*rad*rad;
  clen = par_new[realIdx::c_length];
  rad = par_new[realIdx::radius];
  par_new[realIdx::vol] = M_PI*rad*rad*clen;
  par_new[realIdx::area] = 2.0*M_PI*rad*clen+2.0*M_PI*rad*rad;

}

/**
 * Transfer mesh values to internal concentrations and evaluate chemistry
 * @param grid_vol fluid volume in grid cell that contains biological cell
 * @param npart number of particles in grid cell
 * @param cell_par pointer to cell parameter values
 * @param mesh_vals values of concentrations on mesh
 * @param p_vals values of concentrations in particles
 * @param dt time step interval
 * @param nloop number of subintervals to use for chemistry
 * @param chempar vector containing chemistry parameters
 * @param cell_ipar pointer to integer cell parameter values
 */
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void xferMeshToParticleAndUpdateChem(Real grid_vol, Real npart,
    Real *cell_par, Real *mesh_vals, Real *p_vals, Real dt, int nloop,
    Real *chempar, int *cell_ipar)
{
  Real k1, k2, k3, kr1, kr2, kr3, kg, kv;
  Real fluid_vol = grid_vol/npart;
  int cell_type =  cell_ipar[intIdx::cell_type];
  bool is_tip = (cell_ipar[intIdx::position] == siteLocation::TIP);
  Real radius_max, length_max;

  k1 = chempar[0];
  k2 = chempar[1];
  k3 = chempar[2];
  kr1 = chempar[3];
  kr2 = chempar[4];
  kr3 = chempar[5];
  kg = chempar[6];
  kv = chempar[7];
  radius_max = chempar[8];
  length_max = chempar[9];
  // check to see if particle can grow. Set growth rates to zero
  // if it cannot
    printf("KG: %e KV: %e\n",kg,kv);
  if (cell_type == cellType::FUNGI) {
    if (cell_ipar[intIdx::position] != siteLocation::TIP) {
      if ((cell_par[realIdx::radius] >= radius_max &&
            cell_par[realIdx::c_length] >= length_max)
          || cell_ipar[intIdx::position] == siteLocation::INTERIOR) {
        kg = 0.0;
        kv = 0.0;
      }
    }
  }
    printf("AT 2 KG: %e KV: %e\n",kg,kv);
#if 0
  // Stop growth on interior segments that are maximum sized
  if (cell_par[realIdx::c_length] == length_max &&
      cell_par[realIdx::radius] == radius_max &&
      cell_ipar[intIdx::position] != siteLocation::TIP &&
      cell_type == cellType::FUNGI) {
    k2 = 0.0;
    kr2 = 0.0;
  }
#endif
#if 0
  if (cell_type == cellType::FUNGI &&
      cell_ipar[intIdx::position] == siteLocation::TIP) {
    printf("Rate constants: %e %e id: %d cpu: %d\n",
        k2,kr2, cell_par[realIdx::c_length],
        cell_ipar[intIdx::id], cell_ipar[intIdx::cpu]);
  }
#endif

  // cell parameters
  Real orig_cell_vol = cell_par[realIdx::vol];
//  printf("CELL ID: %d CELL CPU: %d\n",cell_ipar[intIdx::id],cell_ipar[intIdx::cpu]);

  // ********************************************************************
  // Define original concentrations
  // ********************************************************************
  // Fluid concentrations
  Real fA, fB, fC;
  fA = mesh_vals[0];
  fB = mesh_vals[1];
  fC = mesh_vals[2];
#if 0
  if (cell_type == cellType::FUNGI &&
      cell_ipar[intIdx::position] == siteLocation::TIP) {
    printf("mesh concentrations: %e %e %e id: %d cpu: %d\n",
        fA,fB,fC, cell_par[realIdx::c_length],
        cell_ipar[intIdx::id], cell_ipar[intIdx::cpu]);
  }
#endif


#if 0
  printf("   Loop size                     : %16d\n",nloop);
  printf("   Initial radius                : %16.8e\n",cell_par[realIdx::radius]);
  printf("   Initial fluid concentration A : %16.8e\n",fA);
  printf("   Initial fluid concentration B : %16.8e\n",fB);
  printf("   Initial fluid concentration C : %16.8e\n",fC);
#endif

#ifndef AMREX_USE_GPU
    if (fA < 0.0) printf("   internal NEGATIVE fluid concentration A : %16.8e\n",fA);
    if (fB < 0.0) printf("   internal NEGATIVE fluid concentration B : %16.8e\n",fB);
    if (fC < 0.0) printf("   internal NEGATIVE fluid concentration C : %16.8e\n",fC);
#endif

   if (fA < 0. || fB < 0. || fC < 0.) amrex::Abort("Negative fluid concentration into particle integrator");

  // Original cell concentrations
  Real cA, cB, cC;
  cA = p_vals[0];
  cB = p_vals[1];
  cC = p_vals[2];
  int i;
  Real dtp = dt/static_cast<Real>(nloop);

  Real fA_tmp = fA;
  Real fB_tmp = fB;
  Real fC_tmp = fC;

  p_vals[3] = cA;
  p_vals[4] = cB;
  p_vals[5] = cC;

  p_vals[6] = 0.0;
  p_vals[7] = 0.0;
  p_vals[8] = 0.0;

  Real new_cell_vol = cell_par[realIdx::vol];
  Real new_cell_area = cell_par[realIdx::area];
#if 0
  if (cell_type == cellType::FUNGI &&
      cell_ipar[intIdx::position] == siteLocation::TIP) {
    printf("cell concentrations: %e %e %e vol: %e area: %e id: %d cpu: %d\n",
        cA,cB,cC,new_cell_vol,new_cell_area,cell_par[realIdx::c_length],
        cell_ipar[intIdx::id], cell_ipar[intIdx::cpu]);
  }
#endif
  bool is_negative = false;
  // ********************************************************************
  // Take multiple shorter steps to evaluate chemistry
  // ********************************************************************
  for (i=0; i<nloop; i++) {
    // ********************************************************************
    // Compute exchange with fluid
    // 1) Update cA, cB, cC
    // 2) Update fA, fB, fC (approximate)
    // ********************************************************************
    Real dA1, dB1, dC1;
    dA1 = 0.5*dtp*new_cell_area*(k1*fA_tmp-kr1*cA);
    dB1 = 0.0;
    dC1 = 0.5*dtp*new_cell_area*(k3*fC_tmp-kr3*cC);

    // Update cell values
    cA += dA1/orig_cell_vol;
    cB += dB1/orig_cell_vol;
    cC += dC1/orig_cell_vol;

    // Adjusted fluid values. Grid vol represents the volume of fluid in the
    // grid cell not occupied by the particle
    Real dfA = dA1/fluid_vol;
    Real dfB = dB1/fluid_vol;
    Real dfC = dC1/fluid_vol;
    //if (dfA > fA_tmp) dfA = fA_tmp;
    //if (dfB > fB_tmp) dfB = fB_tmp;
    //if (dfC > fC_tmp) dfC = fC_tmp;
    fA_tmp -= dfA;
    fB_tmp -= dfB;
    fC_tmp -= dfC;
    if (cA < 0.0) {
      cA -= dA1/orig_cell_vol;
      fA_tmp += dfA;
      dA1 = -cA*orig_cell_vol;
      cA = 0.0;
      dfA = dA1/fluid_vol;
      fA_tmp -= dfA;
    }
    if (fA_tmp < 0.0) {
      cA -= dA1/orig_cell_vol;
      fA_tmp += dfA;
      dA1 = fA_tmp*fluid_vol;
      dfA = dA1/fluid_vol;
      fA_tmp = 0.0;
      cA += dA1/orig_cell_vol;
    }
    if (cB < 0.0) {
      cB -= dB1/orig_cell_vol;
      fB_tmp += dfB;
      dB1 = -cB*orig_cell_vol;
      cB = 0.0;
      dfB = dB1/fluid_vol;
      fB_tmp -= dfB;
    }
    if (fB_tmp < 0.0) {
      cB -= dB1/orig_cell_vol;
      fB_tmp += dfB;
      dB1 = fB_tmp*fluid_vol;
      dfB = dB1/fluid_vol;
      fB_tmp = 0.0;
      cB += dB1/orig_cell_vol;
    }
    if (cC < 0.0) {
      cC -= dC1/orig_cell_vol;
      fC_tmp += dfC;
      dC1 = -cC*orig_cell_vol;
      cC = 0.0;
      dfC = dC1/fluid_vol;
      fC_tmp -= dfC;
    }
    if (fC_tmp < 0.0) {
      cC -= dC1/orig_cell_vol;
      fC_tmp += dfC;
      dC1 = fC_tmp*fluid_vol;
      dfC = dC1/fluid_vol;
      fC_tmp = 0.0;
      cC += dC1/orig_cell_vol;
    }

#ifndef AMREX_USE_GPU
    if (cA < 0.0) printf("   Internal concentration A : %16.8e\n",cA);
    if (cB < 0.0) printf("   Internal concentration B : %16.8e\n",cB);
    if (cC < 0.0) printf("   Internal concentration C : %16.8e\n",cC);
    if (cA < 0.0) is_negative = true;
    if (cB < 0.0) is_negative = true;
    if (cC < 0.0) is_negative = true;
#endif

    if (fA_tmp < 0. || fB_tmp < 0. || fC_tmp < 0.) 
    {
#ifndef AMREX_USE_GPU
      std::cout << "Negative F_tmp values at first update " << fA_tmp << " " << fB_tmp << " " << fC_tmp << std::endl;
#endif
      amrex::Abort("Negative fluid concentration at first update");
    }

    // ********************************************************************
    // Compute changes due to reactions within the cell only 
    // Update cA, cB, cC
    // ********************************************************************
    Real rA, rB, rC, rV;
    rA = -k2*cA + kr2*cB*cC;
    rB =  k2*cA - kr2*cB*cC;
    rC =  k2*cA - kr2*cB*cC;
    printf("RV: %e KV: %e CB: %e V_cell: %e RB: %e\n",rV,kv,cB,orig_cell_vol,rB);

    // Increment concentrations
    Real dA3, dB3, dC3;
    dA3 = dtp*rA;
    dB3 = dtp*rB;
    dC3 = dtp*rC;
    cA += dA3;
    cB += dB3;
    cC += dC3;
    printf("CB at 1: %e\n",cB);
    if (cA < 0.0) {
      cA -= dA3;
      dA3 = -cA;
      dB3 = -dA3;
      dC3 = -dA3;
      cA += dA3;
      cB += dB3;
      cC += dC3;
    }
    if (cB < 0.0) {
      cB -= dB3;
      dB3 = -cB;
      dA3 = -dB3;
      dC3 = dB3;
      cA += dA3;
      cB += dB3;
      cC += dC3;
    }
    if (cC < 0.0) {
      cC -= dC3;
      dC3 = -cC;
      dA3 = -dC3;
      dB3 = dC3;
      cA += dA3;
      cB += dB3;
      cC += dC3;
    }

    // ********************************************************************
    // Increment B again to account for first order reaction. Handle this
    // reaction separately since it does not conserve A and B
    // ********************************************************************
    printf("CB at 2: %e\n",cB);
    rV =  kv*cB*orig_cell_vol;
    rB = -kg*cB;
    dB3 = dtp*rB;
    cB += dB3;
    Real dV = dtp*rV;
    if (cB < 0.0) {
      cB -= dB3;
      // Assume a shorter time step if cB goes to zero for a full time step
      Real ttmp = -cB/rB;
      dV = ttmp*rV;
      cB = 0.0;
    }
    printf("CB at 3: %e\n",cB);
    printf("DV: %e\n",dV);

#ifndef AMREX_USE_GPU
    if (cA < 0.0) printf("   Intermediate internal concentration A: %16.8e\n",cA);
    if (cB < 0.0) printf("   Intermediate internal concentration B: %16.8e\n",cB);
    if (cC < 0.0) printf("   Intermediate internal concentration C: %16.8e\n",cC);
    if (cA < 0.0) is_negative = true;
    if (cB < 0.0) is_negative = true;
    if (cC < 0.0) is_negative = true;
#endif

    // ********************************************************************
    // Increment volume 
    // ********************************************************************



    // Update the cell area, change in area, and sizes
    if (cell_type == cellType::YEAST) {
      // Yeast cells grow no matter what
      Real dvdt = rV;

      cell_par[realIdx::dvdt] = dvdt; 
      new_cell_vol += dV;
      cell_par[realIdx::vol] = new_cell_vol; 
      // Calculate new surface area, volume and radius
      Real radius = pow((3.0*new_cell_vol/(4.0*M_PI)),1.0/3.0);
      cell_par[realIdx::area] = 4.0*M_PI*radius*radius;
      new_cell_area = cell_par[realIdx::area];
      cell_par[realIdx::dadt] = 2.0*dvdt/radius;
      cell_par[realIdx::radius] = radius;
    } else if (cell_type == cellType::FUNGI) {
      // Fungi segments only grow if they are a growing tip
      // or any of the segment dimensions are below maximum
      // values
      if (is_tip) {
        Real dvdt = std::max(rV,0.0);
        cell_par[realIdx::dvdt] = rV; 
        new_cell_vol += dV;
        cell_par[realIdx::vol] = new_cell_vol; 
        Real radius, c_length, tau;
        if (cell_par[realIdx::radius] < radius_max
            && cell_par[realIdx::c_length] < length_max) {
          // grow segment proportionally
          tau = pow(new_cell_vol/orig_cell_vol,1.0/3.0);
          radius = tau*cell_par[realIdx::radius];
          c_length = tau*cell_par[realIdx::c_length];
          Real dtdt = 1.0/(3.0*pow(new_cell_vol,2.0/3.0))
            *rV/pow(orig_cell_vol,1.0/3.0);
          Real drdt = dtdt*cell_par[realIdx::radius];
          Real dldt = dtdt*cell_par[realIdx::c_length];
          cell_par[realIdx::dadt] = 2.0*M_PI*(drdt*(radius+c_length)+radius*(drdt+dldt));
          if (cell_ipar[intIdx::n_bnds] == 0) {
            cell_par[realIdx::area] = 2.0*M_PI*radius*(radius+c_length);
          } else {
            // Only one end of the cylinder is exposed to external fluid
            cell_par[realIdx::area] = M_PI*radius*(radius+2.0*c_length);
          }
//          printf("Final radius: %e length: %e id: %d cpu: %d\n",
//              cell_par[realIdx::radius],cell_par[realIdx::c_length],
//              cell_ipar[intIdx::id],cell_ipar[intIdx::cpu]);
        } else if (cell_par[realIdx::radius] < radius_max) {
          // increase radius only (this case probably never reached)
          tau = sqrt(new_cell_vol/orig_cell_vol);
          radius = tau*cell_par[realIdx::radius];
          c_length = cell_par[realIdx::c_length];
          Real dtdt = 1.0/(2.0*sqrt(new_cell_vol))*rV/sqrt(orig_cell_vol);
          Real drdt = dtdt*cell_par[realIdx::radius];
          cell_par[realIdx::dadt] = 2.0*M_PI*(2.0*radius+c_length)*drdt;
        } else {
          // increase length
          tau = new_cell_vol/orig_cell_vol;
          radius = cell_par[realIdx::radius];
          c_length = tau*cell_par[realIdx::c_length];
          Real dldt = rV*cell_par[realIdx::c_length]/orig_cell_vol;
          cell_par[realIdx::dadt] =  2.0*M_PI*radius*dldt;
        }
        cell_par[realIdx::radius] = radius;
        cell_par[realIdx::c_length] = c_length;
      } else {
        // Fungi segment can only grow radially
        if (cell_par[realIdx::radius] < radius_max) {
          cell_par[realIdx::dvdt] = rV; 
          new_cell_vol += dV;
          cell_par[realIdx::vol] = new_cell_vol; 
          Real tau = sqrt(new_cell_vol/orig_cell_vol);
          Real radius = tau*cell_par[realIdx::radius];
          Real c_length = cell_par[realIdx::c_length];
          Real dtdt = 1.0/(2.0*sqrt(new_cell_vol))*rV/sqrt(orig_cell_vol);
          Real drdt = dtdt*cell_par[realIdx::radius];
          cell_par[realIdx::radius] = radius;
          cell_par[realIdx::c_length] = c_length;
          // Both end caps of cylinder are hidden from fluid
          cell_par[realIdx::area] = 2.0*M_PI*radius*c_length;
          cell_par[realIdx::dadt] = 2.0*M_PI*(2.0*radius+c_length)*drdt;
        } else {
          // Cell volume doesn't change, reset everything to original values
          cell_par[realIdx::dvdt] = 0.0; 
          cell_par[realIdx::vol] = orig_cell_vol; 
          cell_par[realIdx::dadt] = 0.0;
        }
      }
    }

    // Adjust concentrations for change in volume
    Real ratio = orig_cell_vol/new_cell_vol;
    cA *= ratio;
    cB *= ratio;
    cC *= ratio;

    // cell parameters

    // ********************************************************************
    // Compute second exchange with fluid
    // Update cA, cB, cC
    // ********************************************************************
    Real dA2, dB2, dC2;
    dA2 = 0.5*dtp*new_cell_area*(k1*fA_tmp-kr1*cA);
    dB2 = 0.0;
    dC2 = 0.5*dtp*new_cell_area*(k3*fC_tmp-kr3*cC);

    cA += dA2/new_cell_vol;
    cB += dB2/new_cell_vol;
    cC += dC2/new_cell_vol;

    // ********************************************************************
    // Save fluid concentration increments - we will divide by grid volume 
    //      (and not multiply by dt) when we add it to the fluid 
    // ********************************************************************
    dfA = dA2/fluid_vol;
    dfB = dB2/fluid_vol;
    dfC = dC2/fluid_vol;
    //if (dfA > fA_tmp) dfA = fA_tmp;
    //if (dfB > fB_tmp) dfB = fB_tmp;
    //if (dfC > fC_tmp) dfC = fC_tmp;
    fA_tmp -= dfA;
    fB_tmp -= dfB;
    fC_tmp -= dfC;
    if (cA < 0.0) {
      cA -= dA2/new_cell_vol;
      fA_tmp += dfA;
      dA2 = -cA*new_cell_vol;
      cA = 0.0;
      dfA = dA2/fluid_vol;
      fA_tmp -= dfA;
    }
    if (fA_tmp < 0.0) {
      cA -= dA2/new_cell_vol;
      fA_tmp += dfA;
      dA2 = - fA_tmp*fluid_vol;
      dfA = -fA_tmp;
      fA_tmp = 0.0;
      cA += dA2/new_cell_vol;
    }
    if (cB < 0.0) {
      cB -= dB2/new_cell_vol;
      fB_tmp += dfB;
      dB2 = -cB*new_cell_vol;
      cA = 0.0;
      dfB = dB2/fluid_vol;
      fB_tmp -= dfB;
    }
    if (fB_tmp < 0.0) {
      cB -= dB2/new_cell_vol;
      fB_tmp += dfB;
      dB2 = - fB_tmp*fluid_vol;
      dfB = -fB_tmp;
      fB_tmp = 0.0;
      cB += dB2/new_cell_vol;
    }
    if (cC < 0.0) {
      cC -= dC2/new_cell_vol;
      fC_tmp += dfC;
      dC2 = -cC*new_cell_vol;
      cA = 0.0;
      dfC = dC2/fluid_vol;
      fC_tmp -= dfC;
    }
    if (fC_tmp < 0.0) {
      cC -= dC2/new_cell_vol;
      fC_tmp += dfC;
      dC2 = - fC_tmp*fluid_vol;
      dfC = -fC_tmp;
      fC_tmp = 0.0;
      cC += dC2/new_cell_vol;
    }

#ifndef AMREX_USE_GPU
    if (cA < 0.0) printf("   Final internal concentration A: %16.8e\n",cA);
    if (cB < 0.0) printf("   Final internal concentration B: %16.8e\n",cB);
    if (cC < 0.0) printf("   Final internal concentration C: %16.8e\n",cC);
    if (cA < 0.0) is_negative = true;
    if (cB < 0.0) is_negative = true;
    if (cC < 0.0) is_negative = true;
#endif

    if (fA_tmp < 0. || fB_tmp < 0. || fC_tmp < 0.) 
    {
#ifndef AMREX_USE_GPU
      std::cout << "Negative F_tmp values at second update " << fA_tmp << " " << fB_tmp << " " << fC_tmp << std::endl;
#endif
      amrex::Abort("Negative fluid concentration at second update");
    }

    p_vals[6] += -(dA1+dA2);
    p_vals[7] += -(dB1+dB2);
    p_vals[8] += -(dC1+dC2);

    orig_cell_vol = new_cell_vol;
  }
#if 0
  printf("   Final mesh concentration A    : %16.8e\n",fA_tmp);
  printf("   Final mesh concentration B    : %16.8e\n",fB_tmp);
  printf("   Final mesh concentration C    : %16.8e\n",fC_tmp);
  printf("   Final cell volume        :  %16.8e\n",new_cell_vol);
  printf("   Final cell concentration A: %16.8e\n",cA);
  printf("   Final cell concentration B: %16.8e\n",cB);
  printf("   Final cell concentration C: %16.8e\n",cC);
#endif
#if 0
  printf("   Increment for concentration A : %16.8e grid_vol: %16.8e\n",p_vals[6],grid_vol);
  printf("   Increment for concentration B : %16.8e grid_vol: %16.8e\n",p_vals[7],grid_vol);
  printf("   Increment for concentration C : %16.8e grid_vol: %16.8e\n",p_vals[8],grid_vol);
#endif
  if (is_negative) {
#ifndef AMREX_USE_GPU
    std::cout << "Negative concentrations found. Aborting..." << std::endl;
#endif
    amrex::Abort("Negative concentrations found. Aborting...");
  }
#if 0
  if (cell_type == cellType::FUNGI &&
      cell_ipar[intIdx::position] == siteLocation::TIP) {
    printf("Volume growth rate: %e length: %e id: %d cpu: %d\n",
        cell_par[realIdx::dvdt], cell_par[realIdx::c_length],
        cell_ipar[intIdx::id], cell_ipar[intIdx::cpu]);
  }
#endif

  // ********************************************************************
  // Now update the actual particle values 
  // ********************************************************************
  p_vals[3] = cA;
  p_vals[4] = cB;
  p_vals[5] = cC;
  p_vals[0] = p_vals[3];
  p_vals[1] = p_vals[4];
  p_vals[2] = p_vals[5];
#if 0
  printf("   Final concentration increment for A: %16.8e\n",p_vals[6]);
  printf("   Final concentration increment for B: %16.8e\n",p_vals[7]);
  printf("   Final concentration increment for C: %16.8e\n",p_vals[8]);
#endif
}

/**
 * Initialize exchange increments to zero
 * @param par real parameters describing particle
 * @param ipar integer parameters describing particle
 */
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void 
initExchange(Real *par, int *ipar)
{
  if (ipar[intIdx::cell_type] == cellType::FUNGI) {
    Real *cnc = &par[realIdx::first_data];
    cnc[3] = 0.0;
    cnc[4] = 0.0;
    cnc[5] = 0.0;
  }
}
/**
 * Calculate the exchange of material between two segments. r12 is computed
 * outside this routine and periodic boundary conditions, if applicable, have
 * already been applied
 * @param par1 real parameters describing particle 1
 * @param par2 real parameters describing particle 2
 * @param ipar1 integer parameters describing particle 1
 * @param ipar2 integer parameters describing particle 2
 * @param xpar array of params used in exchange calculation
 * @param dt magnitude of the time step
 */
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void 
evaluateExchange(Real *par1, Real *par2, int *ipar1, int *ipar2, Real *xpar, Real dt)
{
  if (ipar1[intIdx::cell_type] == cellType::FUNGI &&
      ipar2[intIdx::cell_type] == cellType::FUNGI) {
    // Mass exchange parameters
    Real mtA = xpar[0];
    Real mtB = xpar[1];
    Real mtC = xpar[2];
    int i_id[4];
    int i_cpu[4];
    int j_id[4];
    int j_cpu[4];

    int i_bnds = ipar1[intIdx::n_bnds];
    i_id[0] = ipar1[intIdx::seg1_id1];
    i_id[1] = ipar1[intIdx::seg2_id1];
    i_id[2] = ipar1[intIdx::seg3_id1];
    i_id[3] = ipar1[intIdx::seg4_id1];

    i_cpu[0] = ipar1[intIdx::seg1_id2];
    i_cpu[1] = ipar1[intIdx::seg2_id2];
    i_cpu[2] = ipar1[intIdx::seg3_id2];
    i_cpu[3] = ipar1[intIdx::seg4_id2];

    int j_bnds = ipar2[intIdx::n_bnds];
    j_id[0] = ipar2[intIdx::seg1_id1];
    j_id[1] = ipar2[intIdx::seg2_id1];
    j_id[2] = ipar2[intIdx::seg3_id1];
    j_id[3] = ipar2[intIdx::seg4_id1];

    j_cpu[0] = ipar2[intIdx::seg1_id2];
    j_cpu[1] = ipar2[intIdx::seg2_id2];
    j_cpu[2] = ipar2[intIdx::seg3_id2];
    j_cpu[3] = ipar2[intIdx::seg4_id2];

    int iid = ipar1[intIdx::id];
    int icpu = ipar1[intIdx::cpu];
    int jid = ipar2[intIdx::id];
    int jcpu = ipar2[intIdx::cpu];

    // Look for bonds connecting these two segments
    int ib, jb;
    for (ib=0; ib<i_bnds; ib++) {
      if (i_id[ib] == jid && i_cpu[ib] == jcpu) break;
    }
    for (jb=0; jb<j_bnds; jb++) {
      if (j_id[jb] == iid && j_cpu[jb] == icpu) break;
    }
    if (ib <i_bnds && jb < j_bnds) {
      // Found a connection so calculate the exchange increment
       Real *cnc1 = &par1[realIdx::first_data];
       Real *cnc2 = &par2[realIdx::first_data];
       Real cA_1 = cnc1[0];
       Real cB_1 = cnc1[1];
       Real cC_1 = cnc1[2];
       Real cA_2 = cnc2[0];
       Real cB_2 = cnc2[1];
       Real cC_2 = cnc2[2];
       
       // Calculate area of septum
       Real rad1 = par1[realIdx::radius];
       Real area1 = M_PI*rad1*rad1;
       Real rad2 = par2[realIdx::radius];
       Real area2 = M_PI*rad2*rad2;
       Real area = 0.5*(area1+area2);

       // Calculate exchange increment. Exchange is transfer from segment 1
       // to segment 2
       Real dA = (cA_1-cA_2)*area*mtA*dt;
       Real dB = (cB_1-cB_2)*area*mtB*dt;
       Real dC = (cC_1-cC_2)*area*mtC*dt;

       // Update increments
       cnc1[3] -= dA;
       cnc1[4] -= dB;
       cnc1[5] -= dC;
       cnc2[3] += dA;
       cnc2[4] += dB;
       cnc2[5] += dC;
    } // if (i<i_bnds && j_bnds) ...
  }
}

/**
 * apply increments from particle exchange calculation to concentrations
 * @param par real parameters describing particle
 * @param ipar integer parameters describing particle
 */
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void 
applyIncrements(Real *par, int *ipar)
{
  if (ipar[intIdx::cell_type] == cellType::FUNGI) {
    Real *cnc = &par[realIdx::first_data];
    cnc[0] += cnc[3];
    cnc[1] += cnc[4];
    cnc[2] += cnc[5];
  }
}
