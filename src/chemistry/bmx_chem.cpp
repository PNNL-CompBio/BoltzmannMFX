//
//     Copyright (c) 2013 Battelle Memorial Institute
//     Licensed under modified BSD License. A copy of this license can be found
//     in the LICENSE file in the top level directory of this distribution.
//
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
amrex::Real BMXChemistry::kv = 0.0;
amrex::Real BMXChemistry::mtA = 0.0;
amrex::Real BMXChemistry::mtB = 0.0;
amrex::Real BMXChemistry::mtC = 0.0;
amrex::Real BMXChemistry::radius_max = 0.0;
amrex::Real BMXChemistry::length_max = 0.0;
amrex::Real BMXChemistry::p_overlap = 0.0;
amrex::Real BMXChemistry::fusion_prob= 0.0;
amrex::Real BMXChemistry::max_fusion_separation= 0.0;

int BMXChemistry::p_num_reals = 0;
int BMXChemistry::p_num_ints = 0;

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
  p_num_reals = 3*p_num_species;
  p_num_ints = 0;
  p_inc_offset = 2*p_num_species;
  p_verbose = false;
  p_int_vals.push_back(p_num_species);
  p_int_vals.push_back(p_num_ivals);
  p_int_vals.push_back(p_num_reals);
  p_int_vals.push_back(p_num_ints);
  p_int_vals.push_back(p_inc_offset);
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
 * Setup chemistry model by reading a parameter file (a file is not currently
 * being used but probably will be at some point).
 * @param file name of parameter file used by chemistry model
 */
void BMXChemistry::setParams(const char* /*file*/)
{
  ParmParse pp("chem_species");
  pp.get("k1",k1);
  pp.get("kr1",kr1);
  pp.get("k2",k2);
  pp.get("kr2",kr2);
  pp.get("k3",k3);
  pp.get("kr3",kr3);
  pp.get("kg",kg);
  pp.get("kv",kv);
  mtA = 0.0;
  pp.query("mass_transfer_A",mtA);
  mtB = 0.0;
  pp.query("mass_transfer_B",mtB);
  mtC = 0.0;
  pp.query("mass_transfer_C",mtC);
  fusion_prob = 0.0;
  pp.query("fusion_probability",fusion_prob);
  max_fusion_separation = 5.0e-4;
  pp.query("max_fusion_separation",max_fusion_separation);
  p_overlap = 0.2;
  /* figure out cutoff for neighbor list */
  Real vol = SPECIES::max_vol;
  radius_max = SPECIES::max_rad;
  length_max = SPECIES::max_len;
  Real radius = pow((3.0*vol/(4.0*M_PI)),1.0/3.0);
  ParmParse ppF("cell_force");
  Real width;
  ppF.get("neighbor_width",width);
  // TODO: Come up with correct neighborhood value base on what types of cells
  //       are being simulated
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
 * @param id particle index
 * @param p_vals values of concentrations in particles
 * @param p_par values of particle parameters
 */
void BMXChemistry::printCellConcentrations(int id, Real *p_vals, Real *p_par)
{
  if (p_verbose) {
    printf("\n");
    printf("        Particle: %d\n",id);
    if (p_vals[0] < 0) {
      printf(" XXXX   Concentration A: %18.6e\n",p_vals[0]);
    } else {
      printf("        Concentration A: %18.6e\n",p_vals[0]);
    }
    if (p_vals[1] < 0) {
      printf(" XXXX   Concentration B: %18.6e\n",p_vals[1]);
    } else {
      printf("        Concentration B: %18.6e\n",p_vals[1]);
    }
    if (p_vals[2] < 0) {
      printf(" XXXX   Concentration C: %18.6e\n",p_vals[2]);
    } else {
      printf("        Concentration C: %18.6e\n",p_vals[2]);
    }
   // printf("        Concentration B: %18.6e\n",p_vals[1]);
   // printf("        Concentration C: %18.6e\n",p_vals[2]);
    printf("        Cell volume    : %18.6e\n",p_par[realIdx::vol]);
   // printf("        X velocity     : %18.6e\n",p_par[realIdx::velx]);
   // printf("        Y velocity     : %18.6e\n",p_par[realIdx::vely]);
   // printf("        Z velocity     : %18.6e\n",p_par[realIdx::velz]);
  }
}

/**
 * Return values that are also stored in particle idata fields
 * @param idx intIdx parameter representing index of desired value
 * @return integer value at index location
 */
int BMXChemistry::getIntData(int idx)
{
  return p_int_vals[idx];
}

/**
 * Return a vector containing all chemical paramters needed to integrate
 * chemistry
 * @param chempar vector of parameters used in chemistry
 */
void BMXChemistry::getChemParams(std::vector<Real> &chempar)
{
  chempar.clear();
  chempar.push_back(k1);
  chempar.push_back(k2);
  chempar.push_back(k3);
  chempar.push_back(kr1);
  chempar.push_back(kr2);
  chempar.push_back(kr3);
  chempar.push_back(kg);
  chempar.push_back(kv);
  chempar.push_back(radius_max);
  chempar.push_back(length_max);
}

/** Returen a vector containing exchange parameters
 * @return vector of exchange parameters
 */
std::vector<Real> BMXChemistry::getExchangeParameters()
{
  std::vector<Real> ret;
  ret.push_back(mtA);
  ret.push_back(mtB);
  ret.push_back(mtC);
  return ret;
}

/**
 * Return a vector containing parameters for segment fusion
 * @return vector of fusion parameters
 */
std::vector<Real> BMXChemistry::getFusionParameters()
{
  std::vector<Real> ret;
  ret.push_back(fusion_prob);
  ret.push_back(max_fusion_separation);
  return ret;
}
