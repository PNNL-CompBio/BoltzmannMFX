#include <stdio.h>
#include <bmx.H>
#include <bmx_fluid_parms.H>
#include <bmx_chem.H>

amrex::Real BMXChemistry::k1 = 0.0;
amrex::Real BMXChemistry::k2 = 0.0;
amrex::Real BMXChemistry::k3 = 0.0;
amrex::Real BMXChemistry::kr1 = 0.0;
amrex::Real BMXChemistry::kr2 = 0.0;
amrex::Real BMXChemistry::kr3 = 0.0;
amrex::Real BMXChemistry::kg = 0.0;

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
  p_verbose = false;
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
  ParmParse pp("chem_species");
  pp.get("k1",k1);
  pp.get("kr1",kr1);
  pp.get("k2",k2);
  pp.get("kr2",kr2);
  pp.get("k3",k3);
  pp.get("kr3",kr3);
  pp.get("kg",kg);
  p_overlap = 0.2;
  /* figure out cutoff for neighbor list */
  Real vol = SPECIES::max_vol;
  Real radius = pow((3.0*vol/(4.0*M_PI)),1.0/3.0);
  ParmParse ppF("cell_force");
  Real width;
  ppF.get("boundary_width",width);
  DEM::neighborhood = 1.1*(2.0*radius+width);
  ParmParse ppV("bmx");
  int verbose = 0;
  ppV.query("verbose",verbose);
  if (verbose != 0) {
    p_verbose = true;
  } else {
    p_verbose = false;
  }

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
 * Print concentrations of chemical species in cell
 * @param p_vals values of concentrations in particles
 * @param p_par values of particle parameters
 */
void BMXChemistry::printCellConcentrations(Real *p_vals, Real *p_par)
{
  if (p_verbose) {
    printf("\n");
    printf("        Concentration A: %18.6e\n",p_vals[0]);
    printf("        Concentration B: %18.6e\n",p_vals[1]);
    printf("        Concentration C: %18.6e\n",p_vals[2]);
    printf("        Cell volume    : %18.6e\n",p_par[realIdx::vol]);
    printf("        X velocity     : %18.6e\n",p_par[realIdx::velx]);
    printf("        Y velocity     : %18.6e\n",p_par[realIdx::vely]);
    printf("        Z velocity     : %18.6e\n",p_par[realIdx::velz]);
  }
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
 * Check to see if particle meets criteria for splitting
 * @param p_par particle parameters
 * @param p_conc species concentrations
 * @return true if particle should split
 */
bool BMXChemistry::checkSplit(Real *p_par, Real *p_conc)
{
  bool ret = false;
  if (p_par[realIdx::vol] > SPECIES::max_vol) ret = true;
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
  setChildParameters(par_orig, ipar_orig, par_new, ipar_new);
  // Find new locations for split particles
  Real radius = par_new[realIdx::a_size];
  Real theta = M_PI * amrex::Random();
  Real phi = 2.0 * M_PI * amrex::Random();
  Real nx = sin(theta)*cos(phi);
  Real ny = sin(theta)*sin(phi);
  Real nz = cos(theta);
  Real scale = 1.0 - p_overlap;
  pos_new[0] = x + scale*0.5*nx*radius;
  pos_new[1] = y + scale*0.5*ny*radius;
  pos_new[2] = z + scale*0.5*nz*radius;
  pos_orig[0] = x - scale*0.5*nx*radius;
  pos_orig[1] = y - scale*0.5*ny*radius;
  pos_orig[2] = z - scale*0.5*nz*radius;
}
