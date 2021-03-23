#ifndef _BMX_CHEM_H_
#define _BMX_CHEM_H_

#include <AMReX.H>
#include <AMReX_REAL.H>
#include <bmx_pc.H>

using amrex::Real;

#define MAX_CHEM_REAL_VAR 31
#define MAX_CHEM_INT_VAR 5

class BMXChemistry
{
  public:

  /**
   * Retrieve instance of BMXChemistry object
   */
  static BMXChemistry* instance();

  /**
   * Set the number of integer and real varaibles used by the chemistry model
   * @param num_ints number of integer variables used by chemistry model
   * @param num_reals number of real variables used by chemistry models
   */
  void setParams(int num_ints, int num_reals);

  /**
   * Setup chemistry model by reading a parameter file
   * @param file name of parameter file used by chemistry model
   */
  void setParams(const char *file);

  /**
   * Set values of particle integer data
   * @param ipar pointer to internal integer values for each particle
   */
  void setIntegers(int *ipar);

  /**
   * Return the number of Real and int variables
   * @param num_ints number of integer variables in cell model
   * @param num_reals number of Real variables in cell model
   * @param tot_ints total number of ints needed by BMXChemistry (may include
   *        utility data)
   * @param tot_reals total number of Reals needed by BMXChemistry (may include
   *        utility data)
   */
  void getVarArraySizes(int *num_ints, int *num_reals, int *tot_ints, int *tot_reals);

  /**
   * Transfer mesh values to internal concentrations
   * @param grid_vol volume of grid cell that contains biological cell
   * @param cell_vol volume of biological cell
   * @param cell_area surface area of biological cell
   * @param mesh_vals values of concentrations on mesh
   * @param p_vals values of concentrations in particles
   * @param dt time step interval
   */
  void xferMeshToParticle(Real grid_vol, Real cell_vol, Real cell_area,
      Real *mesh_vals, Real *p_vals, Real dt);

  /**
   * Evaluate chemical rate of change inside chemistry module
   * @param p_vals current values of concentrations in particles
   * @param dt time step interval
   */
  void updateChemistry(Real *p_vals, Real dt);

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
  void xferParticleToMesh(Real grid_vol, Real cell_vol, Real cell_area,
      Real *mesh_vals, Real *p_vals, Real dt);

  /**
   * Print concentrations of chemical species in cell
   * @param p_vals values of concentrations in particles
   */
  void printCellConcentrations(Real *p_vals);

  protected:
  /**
   * Constructor
   */
  BMXChemistry();

  /**
   * Destructor
   */
  ~BMXChemistry();

  private:

  int p_num_species; // Total number chemical species in model
  int p_num_ivals;   // Total number integer settings in model
  int p_num_reals;   // Total number of real values allocated in model
  int p_num_ints;    // Total number of int values allocated in model
  int p_inc_offset;  // Offset of increments relative to start of species data

  // Hardwired variables for ABC model
  Real k1, kr1;   // rate coefficients for species A into and out of biological cell
  Real k2, kr2;   // rate coefficients for converting A into B and C inside biological cell
  Real k3, kr3;   // rate coefficients for species C into and out of biological cell

  // Value of pi
  Real p_pi;
  static BMXChemistry *p_instance;
};
#endif
