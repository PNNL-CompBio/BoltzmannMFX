/**
 * Check to see if particle meets criteria for splitting
 * @param p_par real particle parameters
 * @param p_ipar integer particle parameters
 * @param max_vol maximum volume for YEAST splitting
 * @param max_len maximum segment length for splitting of FUNGI
 * @param max_rad maximum radius of FUNGI segments
 * @param engine object for implementing random number generator
 * @return true if particle should split
 */
AMREX_GPU_HOST_DEVICE AMREX_INLINE
bool 
checkSplit(Real *p_par, int *p_ipar, Real max_vol, Real max_len,
    Real max_rad, amrex::RandomEngine const& engine)
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
    } else if (position == siteLocation::SECOND_1 &&
               position == siteLocation::SECOND_2) {
      // Segment is immediately behind the growth tip. Generate side
      // branch with some probability
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
 * @param id, cpu ID and cpu of new segment
 * @param engine object for implementing random number generator
 */
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void 
setNewSegment(Real *pos_orig, Real *pos_new, Real *par_orig,
              Real *par_new, int *ipar_orig, int *ipar_new, int num_reals,
              int num_ints, Real split_len, int id, int cpu,
              amrex::RandomEngine const& engine)
{
  // Copy values from original particle to child
  copyChildParameters(par_orig, ipar_orig, par_new, ipar_new,
      num_reals, num_ints);

  // Get chemical gradient
  Real gx = par_orig[realIdx::gx];
  Real gy = par_orig[realIdx::gy];
  Real gz = par_orig[realIdx::gz];

  //Check out if this is an interior or terminal segment
  int nbnd = ipar_orig[intIdx::n_bnds];
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
    // Fix up topology of new and old segment
    //        1-----------2 parent
    //            /    \
    //           /      \
    //(child) 1----2  1----2 (parent)
    ipar_orig[intIdx::n_bnds] = 1;
    ipar_orig[intIdx::seg1_id1] = id;
    ipar_orig[intIdx::seg1_id2] = cpu;
    // Growth is towards site 2
    ipar_orig[intIdx::site1] = 1;
    ipar_orig[intIdx::position] = siteLocation::TIP;
    ipar_orig[intIdx::growth_dir] = 0;

    ipar_new[intIdx::n_bnds] = 1;
    ipar_new[intIdx::seg1_id1] = ipar_orig[intIdx::id];
    ipar_new[intIdx::seg1_id2] = ipar_orig[intIdx::cpu];
    ipar_new[intIdx::site1] = 2;
    ipar_new[intIdx::position] = siteLocation::TIP;
    // Growth is towards site 1
    ipar_new[intIdx::growth_dir] = 1;
    ipar_new[intIdx::id] = id;
    ipar_new[intIdx::cpu] = cpu;
  } else if (nbnd == 1) {
    // Terminal segment
    // Find direction of chemical gradient
    Real gx = par_orig[realIdx::gx];
    Real gy = par_orig[realIdx::gy];
    Real gz = par_orig[realIdx::gz];
    Real gn = sqrt(gx*gx+gy*gy+gz*gz);
    // Find out which end is bonded to a segment
    int bsite = ipar_orig[intIdx::site1];
    printf("NBNDS: %d SITE1: %d SITE2: %d SITE3: %d\n",
        ipar_orig[intIdx::n_bnds],ipar_orig[intIdx::site1],ipar_orig[intIdx::site2],
        ipar_orig[intIdx::site3]);
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
    printf("THETA: %f PHI: %f GX: %e GY: %e GZ: %e GN: %e\n",theta,phi,gx,gy,gz,gn);
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
    printf("GTHETA: %f GPHI: %f GX: %e GY: %e GZ: %e GN: %e\n",theta,phi,gx,gy,gz,gn);
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
      printf("RTHETA: %f RPHI: %f GNX: %f GNY: %f GNZ: %f NX: %f NY: %f NZ: %f\n",
          rtheta,rphi,gnx,gny,gnz,nx,ny,nz);
    // Now have orientation of old segment (nx,ny,nz) and new segment
    // (gnx,gny,gnz). Finish constructing new configuration.
    Real x = pos_orig[0];
    Real y = pos_orig[1];
    Real z = pos_orig[2];
    Real c_length = par_orig[realIdx::c_length];
    Real x1, y1, z1, x2, y2, z2;
    if (bsite == 1) {
      printf("Add to site 2\n");
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
      ipar_orig[intIdx::growth_dir] = 0;
      ipar_orig[intIdx::site2] = 2;
      ipar_orig[intIdx::position] = siteLocation::SECOND_2;
      Real seg_len = c_length-split_len;
      // find location of Site 2 on new segment (site 1 on new segment equals
      // site 2 of old segment)
      Real ax = x2 + gnx*seg_len;
      Real ay = y2 + gny*seg_len;
      Real az = z2 + gnz*seg_len;
      printf("B1 X1: %e Y1: %e Z1: %e X2: %e Y2: %e Z2: %e X3: %e Y3: %e Z3: %e\n",
          x1,y1,z1,x2,y2,z2,ax,ay,az);
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
      printf("SITE: %d GTHETA: %f GPHI: %f\n",bsite,gtheta,gphi);
      par_new[realIdx::theta] = gtheta;
      par_new[realIdx::phi] = gphi;
      ipar_new[intIdx::site1] = 1;
      ipar_new[intIdx::position] = siteLocation::TIP;
    } else {
      printf("Add to site 1\n");
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
      ipar_orig[intIdx::growth_dir] = 1;
      ipar_orig[intIdx::site2] = 1;
      ipar_orig[intIdx::position] = siteLocation::SECOND_1;
      Real seg_len = c_length-split_len;
      // find location of Site 1 on new segment (site 2 on new segment equals
      // site 1 of old segment)
      Real ax = x1 - gnx*seg_len;
      Real ay = y1 - gny*seg_len;
      Real az = z1 - gnz*seg_len;
      printf("B2 X1: %e Y1: %e Z1: %e X2: %e Y2: %e Z2: %e X3: %e Y3: %e Z3: %e\n",
          x1,y1,z1,x2,y2,z2,ax,ay,az);
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
      printf("SITE: %d GTHETA: %f GPHI: %f\n",bsite,gtheta,gphi);
      par_new[realIdx::theta] = gtheta;
      par_new[realIdx::phi] = gphi;
      ipar_new[intIdx::site1] = 2;
      ipar_new[intIdx::position] = siteLocation::TIP;
    }
//    printf("PARENT LENGTH: %e\n",par_orig[realIdx::c_length]);
//    printf("CHILD LENGTH: %e\n",par_new[realIdx::c_length]);
    ipar_orig[intIdx::n_bnds] = 2; 
    printf("Adding new segment ID: %d CPU: %d\n",id,cpu);
    ipar_orig[intIdx::seg2_id1] = id;
    ipar_orig[intIdx::seg2_id2] = cpu;
    ipar_new[intIdx::n_bnds] = 1; 
    ipar_new[intIdx::seg1_id1] = ipar_orig[intIdx::id];
    ipar_new[intIdx::seg1_id2] = ipar_orig[intIdx::cpu];
  } else if (nbnd == 2) {
    // Interior segment just behind growth tip
    if (ipar_orig[intIdx::position] == siteLocation::SECOND_1
        || ipar_orig[intIdx::position] == siteLocation::SECOND_2) {
    }
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
  Real k1, k2, k3, kr1, kr2, kr3, kg;
  Real fluid_vol = grid_vol/npart;
  int cell_type =  cell_ipar[intIdx::cell_type];
  bool is_tip = (cell_ipar[intIdx::n_bnds] <= 1);
  Real radius_max, length_max;

  k1 = chempar[0];
  k2 = chempar[1];
  k3 = chempar[2];
  kr1 = chempar[3];
  kr2 = chempar[4];
  kr3 = chempar[5];
  kg = chempar[6];
  radius_max = chempar[7];
  length_max = chempar[8];

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
  double dtp = dt/static_cast<double>(nloop);

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
    double dfA = dA1/fluid_vol;
    double dfB = dB1/fluid_vol;
    double dfC = dC1/fluid_vol;
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
    Real rA, rB, rC;
    rA = -k2*cA + kr2*cB*cC;
    rB =  k2*cA - kr2*cB*cC;
    rC =  k2*cA - kr2*cB*cC;

    // Increment concentrations
    Real dA3, dB3, dC3;
    dA3 = dtp*rA;
    dB3 = dtp*rB;
    dC3 = dtp*rC;
    cA += dA3;
    cB += dB3;
    cC += dC3;
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
      Real dvdt = std::max(kg*rB*orig_cell_vol,0.0);

      cell_par[realIdx::dvdt] = dvdt; 
      new_cell_vol += dtp*dvdt;
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
        Real dvdt = std::max(kg*rB*orig_cell_vol,0.0);
        cell_par[realIdx::dvdt] = dvdt; 
        new_cell_vol += dtp*dvdt;
        cell_par[realIdx::vol] = new_cell_vol; 
        Real radius, c_length, tau;
        if (cell_par[realIdx::radius] < radius_max
            && cell_par[realIdx::c_length] < length_max) {
          // grow segment proportionally
          tau = pow(new_cell_vol/orig_cell_vol,1.0/3.0);
          radius = tau*cell_par[realIdx::radius];
          c_length = tau*cell_par[realIdx::c_length];
          Real dtdt = 1.0/(3.0*pow(new_cell_vol,2.0/3.0))*dvdt/pow(orig_cell_vol,1.0/3.0);
          Real drdt = dtdt*cell_par[realIdx::radius];
          Real dldt = dtdt*cell_par[realIdx::c_length];
          cell_par[realIdx::dadt] = 2.0*M_PI*(drdt*(radius+c_length)+radius*(drdt+dldt));
          if (cell_ipar[intIdx::n_bnds] == 0) {
            cell_par[realIdx::area] = 2.0*M_PI*radius*(radius+c_length);
          } else {
            // Only one end of the cylinder is exposed to external fluid
            cell_par[realIdx::area] = M_PI*radius*(radius+2.0*c_length);
          }
        } else if (cell_par[realIdx::radius] < radius_max) {
          // increase radius only (this case probably never reached)
          tau = sqrt(new_cell_vol/orig_cell_vol);
          radius = tau*cell_par[realIdx::radius];
          c_length = cell_par[realIdx::c_length];
          Real dtdt = 1.0/(2.0*sqrt(new_cell_vol))*dvdt/sqrt(orig_cell_vol);
          Real drdt = dtdt*cell_par[realIdx::radius];
          cell_par[realIdx::dadt] = 2.0*M_PI*(2.0*radius+c_length)*drdt;
        } else {
          // increase length
          tau = new_cell_vol/orig_cell_vol;
          radius = cell_par[realIdx::radius];
          c_length = tau*cell_par[realIdx::c_length];
          Real dldt = dvdt*cell_par[realIdx::c_length]/orig_cell_vol;
          cell_par[realIdx::dadt] =  2.0*M_PI*radius*dldt;
        }
        cell_par[realIdx::radius] = radius;
        cell_par[realIdx::c_length] = c_length;
      } else {
        // Fungi segment can only grow radially
        if (cell_par[realIdx::radius] < radius_max) {
          Real dvdt = std::max(kg*rB*orig_cell_vol,0.0);
          cell_par[realIdx::dvdt] = dvdt; 
          new_cell_vol += dtp*dvdt;
          cell_par[realIdx::vol] = new_cell_vol; 
          Real tau = sqrt(new_cell_vol/orig_cell_vol);
          Real radius = tau*cell_par[realIdx::radius];
          Real c_length = cell_par[realIdx::c_length];
          Real dtdt = 1.0/(2.0*sqrt(new_cell_vol))*dvdt/sqrt(orig_cell_vol);
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

