#include <stdio.h>
#include <math.h>
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
 * Set the number of integer and real variables used by the chemistry model
 * @param num_ints number of integer variables used by chemistry model
 * @param num_reals number of real variables used by chemistry models
 */
void BMXChemistry::setParams(int num_ints, int num_reals)
{
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
 * @param cell_vol volume of biological cell
 * @param cell_area surface area of biological cell
 * @param mesh_vals values of concentrations on mesh
 * @param p_vals values of concentrations in particles
 * @param dt time step interval
 */
void BMXChemistry::xferMeshToParticle(Real grid_vol, Real cell_vol, Real cell_area,
    Real *mesh_vals, Real *p_vals, Real dt)
{
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
 * @param dt time step interval
 */
void BMXChemistry::updateChemistry(Real *p_vals, Real dt)
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
  // Increment concentrations
  p_vals[0] += dt*fA;
  p_vals[1] += dt*fB;
  p_vals[2] += dt*fC;
}

/**
 * Calculate transfer increments based on current concentrations in biological
 * cell and in grid cell
 * @param grid_vol volume of grid cell that contains biological cell
 * @param cell_vol volume of biological cell
 * @param cell_area surface area of biological cell
 * @param mesh_inc values to increment concentrations on mesh
 * @param p_vals values of concentrations in particles
 * @param dt time step interval
 */
void BMXChemistry::xferParticleToMesh(Real grid_vol, Real cell_vol,
    Real cell_area, Real *mesh_vals, Real *p_vals, Real dt)
{
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
 */
void BMXChemistry::printCellConcentrations(Real *p_vals)
{
  printf("        Concentration A: %18.6f\n",p_vals[0]);
  printf("        Concentration B: %18.6f\n",p_vals[1]);
  printf("        Concentration C: %18.6f\n",p_vals[2]);
}
