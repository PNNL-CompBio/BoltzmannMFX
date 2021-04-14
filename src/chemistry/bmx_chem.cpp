#include <stdio.h>
#include <math.h>
#include <bmx.H>
#include <bmx_fluid_parms.H>
#include <bmx_chem.H>

BMXChemistry *BMXChemistry::p_instance = NULL;

/**
 * Retrieve instance of BMXChemistry object
 */
BMXChemistry* BMXChemistry::instance()
{
  if (p_instance == NULL) {
    p_instance = new BMXChemistry();
  }
  return p_instance;
}

/**
 * Constructor
 */
BMXChemistry::BMXChemistry()
{
  p_num_species = 3;
  p_num_ivals = 0;
  p_num_reals = 2*p_num_species;
  p_num_ints = 0;
  p_inc_offset = p_num_species;
  p_pi = 4.0*atan(1.0);
}

/**
 * Destructor
 */
BMXChemistry::~BMXChemistry()
{
  delete p_instance;
}

/**
 * Set values of particle integer data
 * @param ipar pointer to internal integer values for each particle
 */
void BMXChemistry::setIntegers(int *ipar)
{
  ipar[0] = p_num_species;
  ipar[1] = p_num_ivals;
  ipar[2] = p_num_reals;
  ipar[3] = p_num_ints;
  ipar[4] = p_inc_offset;
}

/**
 * Setup chemistry model by reading a parameter file
 * @param file name of parameter file used by chemistry model
 */
void BMXChemistry::setParams(const char *file)
{
  k1 = 0.01;
  kr1 = 0.01;
  k2 = 0.005;
  kr2 = 0.005;
  k3 = 0.01;
  kr3 = 0.01;
  kg = 2000.0;
  p_overlap = 0.2;
}

/**
 * Return the number of Real and int variables
 * @param num_ints number of integer variables in cell model
 * @param num_reals number of Real variables in cell model
 * @param tot_ints total number of ints needed by BMXChemistry (may include
 *        utility data)
 * @param tot_reals total number of Reals needed by BMXChemistry (may include
 *        utility data)
 */
void BMXChemistry::getVarArraySizes(int *num_ints, int *num_reals, int *tot_ints, int *tot_reals)
{
  *num_reals = p_num_species;
  *num_ints = p_num_ivals;
  *tot_reals = p_num_reals;
  *tot_ints = p_num_ints;
}

/**
 * Transfer mesh values to internal concentrations
 * @param grid_vol volume of grid cell that contains biological cell
 * @param cell_par pointer to cell parameter values
 * @param mesh_vals values of concentrations on mesh
 * @param p_vals values of concentrations in particles
 * @param dt time step interval
 */
void BMXChemistry::xferMeshToParticle(Real grid_vol, Real *cell_par,
    Real *mesh_vals, Real *p_vals, Real dt)
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
void BMXChemistry::updateChemistry(Real *p_vals, Real *cell_par, Real dt)
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
  Real radius = pow((3.0*volume/(4*p_pi)),1.0/3.0);
  cell_par[realIdx::area] = 4.0*p_pi*radius*radius;
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

/**
 * Calculate transfer increments based on current concentrations in biological
 * cell and in grid cell
 * @param grid_vol volume of grid cell that contains biological cell
 * @param cell_par pointer to cell parameter values
 * @param mesh_inc values to increment concentrations on mesh
 * @param p_vals values of concentrations in particles
 * @param dt time step interval
 */
void BMXChemistry::xferParticleToMesh(Real grid_vol, Real *cell_par,
    Real *mesh_vals, Real *p_vals, Real dt)
{
  // cell parameters
  Real cell_vol = cell_par[realIdx::vol];
  Real cell_area = cell_par[realIdx::area];
  // fluid concentrations
  Real fA, fB, fC;
  fA = p_vals[3];
  fB = p_vals[4];
  fC = p_vals[5];
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
  // adjust cell concentrations
  p_vals[0] = p_vals[0]+dA/cell_vol;
  p_vals[1] = p_vals[1]+dB/cell_vol;
  p_vals[2] = p_vals[2]+dC/cell_vol;
  // fluid concentration increments
  mesh_vals[0] = -dA/(dt*grid_vol);
  mesh_vals[1] = -dB/(dt*grid_vol);
  mesh_vals[2] = -dC/(dt*grid_vol);
}

/**
 * Print concentrations of chemical species in cell
 * @param p_vals values of concentrations in particles
 * @param p_par values of particle parameters
 */
void BMXChemistry::printCellConcentrations(Real *p_vals, Real *p_par)
{
  printf("\n");
  printf("        Concentration A: %18.6f\n",p_vals[0]);
  printf("        Concentration B: %18.6f\n",p_vals[1]);
  printf("        Concentration C: %18.6f\n",p_vals[2]);
  printf("        Cell volume    : %18.6f\n",p_par[realIdx::vol]);
  printf("        X velocity     : %18.6f\n",p_par[realIdx::velx]);
  printf("        Y velocity     : %18.6f\n",p_par[realIdx::vely]);
  printf("        Z velocity     : %18.6f\n",p_par[realIdx::velz]);
}

/**
 * Duplicate data from original particle to child when splitting
 * @param p_real_orig pointer to real values from original particle
 * @param p_int_orig pointer to integer values from original particle
 * @param p_real_child pointer to real values on child  particle
 * @param p_int_child pointer to integer values on child particle
 */
void BMXChemistry::setChildParameters(Real *p_real_orig, int *p_int_orig,
    Real *p_real_child, int *p_int_child)
{
  int i;
  int nreals = p_num_reals + realIdx::count-1;
  int nints = p_num_ints + intIdx::count-1;
  for (i=0; i<nreals; i++) p_real_child[i] = p_real_orig[i];
  for (i=0; i<nints; i++) p_int_child[i] = p_int_orig[i];

  // fix up values that need to be modified due to splitting
  Real volume = p_real_orig[realIdx::vol]/2.0;
  Real radius = pow((3.0*volume/(4*p_pi)),1.0/3.0);
  Real area = 4.0*p_pi*radius*radius;
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
 * Check to see if particle meets criteria for splitting
 * @param p_par particle parameters
 * @param p_conc species concentrations
 * @return true if particle should split
 */
bool BMXChemistry::checkSplit(Real *p_par, Real *p_conc)
{
  bool ret = false;
  if (p_par[realIdx::vol] > FLUID::max_vol) ret = true;
  return ret;
}

/**
 * Create positions and parameter values for new particle pair. This consist
 * of the original particle and a new particle
 * @param pos_orig, pos_new positions of original and new particles
 * @param par_orig, par_new real parameter values of original and new particles
 * @param ipar_orig, ipar_new integer parameter values of original and new particles
 */
void BMXChemistry::setNewCell(Real *pos_orig, Real *pos_new, Real *par_orig,
    Real *par_new, int *ipar_orig, int *ipar_new)
{
  Real x,y,z;
  x = pos_orig[0];
  y = pos_orig[1];
  z = pos_orig[2];

  // Copy values from original particle and modify some values
  // as appropriate
  setChildParameters(par_orig, ipar_orig, par_new,ipar_new);
  // Find new locations for split particles
  Real radius = par_new[realIdx::a_size];
  Real theta = p_pi * amrex::Random();
  Real phi = 2.0 * p_pi * amrex::Random();
  Real nx = cos(theta)*cos(phi);
  Real ny = cos(theta)*sin(phi);
  Real nz = sin(theta);
  Real scale = 1.0 - p_overlap;
  pos_new[0] = x + scale*0.5*nx*radius;
  pos_new[1] = y + scale*0.5*ny*radius;
  pos_new[2] = z + scale*0.5*nz*radius;
  pos_orig[0] = x - scale*0.5*nx*radius;
  pos_orig[1] = y - scale*0.5*ny*radius;
  pos_orig[2] = z - scale*0.5*nz*radius;
}
