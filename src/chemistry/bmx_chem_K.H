/**
 * Check to see if particle meets criteria for splitting
 * @param p_par particle parameters
 * @param p_conc species concentrations
 * @return true if particle should split
 */
AMREX_GPU_HOST_DEVICE AMREX_INLINE
bool 
checkSplit(Real *p_par, int cell_type, Real max_vol)
{
  bool ret = false;
  if (cell_type == cellType::YEAST) {
    if (p_par[realIdx::vol] > max_vol) ret = true;
  } else if (cell_type == cellType::FUNGI) {
    //TODO: Check fungi values
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
setChildParameters(Real *p_real_orig, int *p_int_orig,
                   Real *p_real_child, int *p_int_child,
                   int num_reals, int num_ints) 
{
  int i;
  int nreals = num_reals + realIdx::count-1;
  int nints  = num_ints  +  intIdx::count-1;
  for (i=0; i<nreals; i++) p_real_child[i] = p_real_orig[i];
  for (i=0; i<nints; i++) p_int_child[i] = p_int_orig[i];

  // fix up values that need to be modified due to splitting
  if (p_int_orig[intIdx::cell_type] == cellType::YEAST) {
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
    p_real_orig[realIdx::radius] = radius;
  } else if (p_int_orig[intIdx::cell_type] == cellType::FUNGI) {
    //TODO: set values for fungi
  }
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
  // Copy values from original particle and modify some values
  // as appropriate
  setChildParameters(par_orig, ipar_orig, par_new, ipar_new,
      num_reals, num_ints);

  if (ipar_orig[intIdx::cell_type] == cellType::YEAST) {
    Real x,y,z;
    x = pos_orig[0];
    y = pos_orig[1];
    z = pos_orig[2];

    // Find new locations for split particles
    Real radius = par_new[realIdx::radius];
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
  } else if (ipar_orig[intIdx::cell_type] == cellType::FUNGI) {
    //TODO: set fungi positions
  }
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

