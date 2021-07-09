/**
 * Transfer mesh values to internal concentrations and evaluate chemistry
 * @param grid_vol volume of grid cell that contains biological cell
 * @param cell_par pointer to cell parameter values
 * @param mesh_vals values of concentrations on mesh
 * @param p_vals values of concentrations in particles
 * @param dt time step interval
 */
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void xferMeshToParticleAndUpdateChem(Real grid_vol,
    Real *cell_par, Real *mesh_vals, Real *p_vals, Real dt,
    Real k1, Real k2, Real k3, Real kr1, Real kr2, Real kr3,
    Real kg)
{
  // cell parameters
  Real cell_vol = cell_par[realIdx::vol];
  Real cell_area = cell_par[realIdx::area];

  // fluid concentrations
  Real fA, fB, fC;
  fA = mesh_vals[0];
  fB = mesh_vals[1];
  fC = mesh_vals[2];
  
  // cell concentrations
  Real cA, cB, cC;
  cA = p_vals[0];
  cB = p_vals[1];
  cC = p_vals[2];
 
  // incremental changes
  Real dA, dB, dC;
  dA = 0.5*dt*cell_area*(k1*fA-kr1*cA);
  dB = 0.0;
  dC = 0.5*dt*cell_area*(k3*fC-kr3*cC);
  // adjusted cell values
  cA += dA/cell_vol;
  cB += dB/cell_vol;
  cC += dC/cell_vol;
  // adjusted fluid values
  fA = fA-dA/grid_vol;
  fB = fB-dB/grid_vol;
  fC = fC-dC/grid_vol;

  // Rate of change of A, B, C
  Real rA, rB, rC;
  rA = -k2*cA + kr2*cB*cC;
  rB = k2*cA - kr2*cB*cC;
  rC = k2*cA - kr2*cB*cC;
  // Increment volume 
  Real volume = cell_par[realIdx::vol];
  Real dvdt = kg*volume*rB;
  if (rB < 0.0) dvdt = 0.0;
  cell_par[realIdx::dvdt] = dvdt; 
  Real new_vol = volume + dt*dvdt;
  cell_par[realIdx::vol] = new_vol;
  Real radius = pow((3.0*volume/(4.0*M_PI)),1.0/3.0);
  cell_par[realIdx::area] = 4.0*M_PI*radius*radius;
  cell_par[realIdx::dadt] = 2.0*dvdt/radius;
  cell_par[realIdx::a_size] = radius;
  cell_par[realIdx::b_size] = radius;
  cell_par[realIdx::c_size] = radius;
  // Increment concentrations
  cA += dt*rA;
  cB += dt*rB;
  cC += dt*rC;
  // Adjust concentrations for change in volume
  Real ratio = volume/new_vol;
  cA *= ratio;
  cB *= ratio;
  cC *= ratio;

  // cell parameters
  cell_vol = cell_par[realIdx::vol];
  cell_area = cell_par[realIdx::area];
  // incremental changes
  dA = 0.5*dt*cell_area*(k1*fA-kr1*cA);
  dB = 0.0;
  dC = 0.5*dt*cell_area*(k3*fC-kr3*cC);
  // save adjusted cell concentrations
  p_vals[3] = cA+dA/cell_vol;
  p_vals[4] = cB+dB/cell_vol;
  p_vals[5] = cC+dC/cell_vol;
  // save fluid concentration increments

  p_vals[6] = -dA/(dt*grid_vol);
  p_vals[7] = -dB/(dt*grid_vol);
  p_vals[8] = -dC/(dt*grid_vol);

  // Now update the actual particle values 
  p_vals[0] = p_vals[3];
  p_vals[1] = p_vals[4];
  p_vals[2] = p_vals[5];
}

/**
 * Transfer mesh values to internal concentrations
 * @param grid_vol volume of grid cell that contains biological cell
 * @param cell_par pointer to cell parameter values
 * @param mesh_vals values of concentrations on mesh
 * @param p_vals values of concentrations in particles
 * @param dt time step interval
 */
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void xferMeshToParticle(Real grid_vol, Real *cell_par,
    Real *mesh_vals, Real *p_vals, Real dt,
    Real k1, Real k2, Real k3, Real kr1, Real kr2, Real kr3,
    Real kg)
{
  // cell parameters
  Real cell_vol = cell_par[realIdx::vol];
  Real cell_area = cell_par[realIdx::area];
  // fluid concentrations
  Real fA, fB, fC;
  fA = mesh_vals[0];
  fB = mesh_vals[1];
  fC = mesh_vals[2];
  // cell concentrations
  Real cA, cB, cC;
  cA = p_vals[0];
  cB = p_vals[1];
  cC = p_vals[2];
  // incremental changes
  Real dA, dB, dC;
  dA = dt*cell_area*(k1*fA-kr1*cA);
  dB = 0.0;
  dC = dt*cell_area*(k3*fC-kr3*cC);
  // adjusted cell values
  p_vals[0] = p_vals[0]+dA/cell_vol;
  p_vals[1] = p_vals[1]+dB/cell_vol;
  p_vals[2] = p_vals[2]+dC/cell_vol;
  // adjusted fluid values
  p_vals[3] = fA-dA/grid_vol;
  p_vals[4] = fB-dB/grid_vol;
  p_vals[5] = fC-dC/grid_vol;
}

/**
 * Evaluate chemical rate of change inside chemistry module
 * @param pval current values of concentrations in particles
 * @param cell_par pointer to cell parameter values
 * @param dt time step interval
 */
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void updateChemistry(Real *p_vals, Real *cell_par, Real dt,
    Real k1, Real k2, Real k3, Real kr1, Real kr2, Real kr3,
    Real kg)
{
  // Concentrations of A, B, C
  Real A, B, C;
  A = p_vals[0];
  B = p_vals[1];
  C = p_vals[2];
  // Rate of change of A, B, C
  Real fA, fB, fC;
  fA = -k2*A + kr2*B*C;
  fB = k2*A - kr2*B*C;
  fC = k2*A - kr2*B*C;
  // Increment volume 
  Real volume = cell_par[realIdx::vol];
  Real dvdt = kg*volume*fB;
  if (fB < 0.0) dvdt = 0.0;
  cell_par[realIdx::dvdt] = dvdt; 
  Real new_vol = volume + dt*dvdt;
  cell_par[realIdx::vol] = new_vol;
  Real radius = pow((3.0*volume/(4.0*M_PI)),1.0/3.0);
  cell_par[realIdx::area] = 4.0*M_PI*radius*radius;
  cell_par[realIdx::dadt] = 2.0*dvdt/radius;
  cell_par[realIdx::a_size] = radius;
  cell_par[realIdx::b_size] = radius;
  cell_par[realIdx::c_size] = radius;
  // Increment concentrations
  p_vals[0] += dt*fA;
  p_vals[1] += dt*fB;
  p_vals[2] += dt*fC;
  // Adjust concentrations for change in volume
  Real ratio = volume/new_vol;
  p_vals[0] *= ratio;
  p_vals[1] *= ratio;
  p_vals[2] *= ratio;
}

