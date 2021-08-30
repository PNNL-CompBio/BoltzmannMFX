/**
 * Check to see if particle meets criteria for splitting
 * @param p_par particle parameters
 * @param p_conc species concentrations
 * @return true if particle should split
 */
AMREX_GPU_HOST_DEVICE AMREX_INLINE
bool 
checkSplit(Real *p_par, Real* /*p_conc*/, Real max_vol)
{
  bool ret = false;
  if (p_par[realIdx::vol] > max_vol) ret = true;
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
setChildParameters(Real *p_real_orig, int *p_int_orig,
                   Real *p_real_child, int *p_int_child,
                   Real /*overlap*/, int num_reals, int num_ints) 
{
  int i;
  int nreals = num_reals + realIdx::count-1;
  int nints  = num_ints  +  intIdx::count-1;
  for (i=0; i<nreals; i++) p_real_child[i] = p_real_orig[i];
  for (i=0; i<nints; i++) p_int_child[i] = p_int_orig[i];

  // fix up values that need to be modified due to splitting
  Real volume = p_real_orig[realIdx::vol]/2.0;
  Real radius = pow((3.0*volume/(4.0*M_PI)),1.0/3.0);
  Real area = 4.0*M_PI*radius*radius;
  Real dvdt = p_real_orig[realIdx::dvdt];
  Real dadt = 2.0*dvdt/radius;
  p_real_orig[realIdx::vol] = volume;
  p_real_child[realIdx::vol] = volume;
  p_real_orig[realIdx::area] = area;
  p_real_child[realIdx::area] = area;
  p_real_orig[realIdx::dadt] = dadt;
  p_real_child[realIdx::dadt] = dadt;
  p_real_orig[realIdx::a_size] = radius;
  p_real_child[realIdx::a_size] = radius;
  p_real_orig[realIdx::b_size] = radius;
  p_real_child[realIdx::b_size] = radius;
  p_real_orig[realIdx::c_size] = radius;
  p_real_child[realIdx::c_size] = radius;
}

/**
 * Create positions and parameter values for new particle pair. This consist
 * of the original particle and a new particle
 * @param pos_orig, pos_new positions of original and new particles
 * @param par_orig, par_new real parameter values of original and new particles
 * @param ipar_orig, ipar_new integer parameter values of original and new particles
 */
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void 
setNewCell(Real *pos_orig, Real *pos_new, Real *par_orig,
           Real *par_new, int *ipar_orig, int *ipar_new,
           Real overlap, int num_reals, int num_ints, amrex::RandomEngine const& engine)
{
  Real x,y,z;
  x = pos_orig[0];
  y = pos_orig[1];
  z = pos_orig[2];

  // Copy values from original particle and modify some values
  // as appropriate
  setChildParameters(par_orig, ipar_orig, par_new, ipar_new,
                     overlap, num_reals, num_ints);

  // Find new locations for split particles
  Real radius = par_new[realIdx::a_size];
#if 0
  Real theta = M_PI ; // * amrex::Random();
  Real phi = 2.0 * M_PI ; // * amrex::Random();
#else
  Real theta = M_PI * amrex::Random(engine);
  Real phi = 2.0 * M_PI * amrex::Random(engine);
#endif
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
}

/**
 * Transfer mesh values to internal concentrations and evaluate chemistry
 * @param grid_vol volume of grid cell that contains biological cell
 * @param cell_par pointer to cell parameter values
 * @param mesh_vals values of concentrations on mesh
 * @param p_vals values of concentrations in particles
 * @param dt time step interval
 * @param nloop number of subintervals to use for chemistry
 */
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void xferMeshToParticleAndUpdateChem(Real grid_vol,
    Real *cell_par, Real *mesh_vals, Real *p_vals, Real dt, int nloop,
    Real k1, Real k2, Real k3, Real kr1, Real kr2, Real kr3,
    Real kg)
{
  // cell parameters
  Real orig_cell_vol  = cell_par[realIdx::vol];

  // ********************************************************************
  // Define original concentrations
  // ********************************************************************
  // Fluid concentrations
  Real fA, fB, fC;
  fA = mesh_vals[0];
  fB = mesh_vals[1];
  fC = mesh_vals[2];
  
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
#ifndef AMREX_USE_GPU
    if (cA < 0.0) printf("   Initial internal concentration A: %16.8e\n",cA);
    if (cB < 0.0) printf("   Initial internal concentration B: %16.8e\n",cB);
    if (cC < 0.0) printf("   Initial internal concentration C: %16.8e\n",cC);
    if (cA < 0.0) is_negative = true;
    if (cB < 0.0) is_negative = true;
    if (cC < 0.0) is_negative = true;
#endif

    // Adjusted fluid values -- these are just approximations because
    //    we are not taking into account the fluid volume fraction --
    fA_tmp += -dA1/grid_vol;
    fB_tmp += -dB1/grid_vol;
    fC_tmp += -dC1/grid_vol;

    // ********************************************************************
    // Compute changes due to reactions within the cell only 
    // Update cA, cB, cC
    // ********************************************************************
    Real rA, rB, rC;
    rA = -k2*cA + kr2*cB*cC;
    rB =  k2*cA - kr2*cB*cC;
    rC =  k2*cA - kr2*cB*cC;

    // Increment concentrations
    cA += dtp*rA;
    cB += dtp*rB;
    cC += dtp*rC;
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
    Real dvdt = std::max(kg*rB*new_cell_vol,0.0);

    cell_par[realIdx::dvdt] = dvdt; 
    orig_cell_vol = new_cell_vol;
    new_cell_vol += dtp*dvdt;
    cell_par[realIdx::vol] = new_cell_vol; 


    // Update the cell area, change in area, and sizes
    Real radius = pow((3.0*new_cell_vol/(4.0*M_PI)),1.0/3.0);
    cell_par[realIdx::area] = 4.0*M_PI*radius*radius;
    new_cell_area = cell_par[realIdx::area];
    cell_par[realIdx::dadt] = 2.0*dvdt/radius;
    cell_par[realIdx::a_size] = radius;
    cell_par[realIdx::b_size] = radius;
    cell_par[realIdx::c_size] = radius;

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

#ifndef AMREX_USE_GPU
    if (cA < 0.0) printf("   Final internal concentration A: %16.8e\n",cA);
    if (cB < 0.0) printf("   Final internal concentration B: %16.8e\n",cB);
    if (cC < 0.0) printf("   Final internal concentration C: %16.8e\n",cC);
    if (cA < 0.0) is_negative = true;
    if (cB < 0.0) is_negative = true;
    if (cC < 0.0) is_negative = true;
#endif

    // ********************************************************************
    // Save fluid concentration increments - we will divide by grid volume 
    //      (and not multiply by dt) when we add it to the fluid 
    // ********************************************************************
    fA_tmp += -dA2/grid_vol;
    fB_tmp += -dB2/grid_vol;
    fC_tmp += -dC2/grid_vol;

    p_vals[6] += -(dA1+dA2);
    p_vals[7] += -(dB1+dB2);
    p_vals[8] += -(dC1+dC2);
  }
  if (is_negative) {
    std::cout << "Negative concentrations found. Aborting..." << std::endl;
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
}

