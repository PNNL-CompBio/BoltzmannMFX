#ifndef _BMX_CHEM_H_
#define _BMX_CHEM_H_

#include <AMReX.H>
#include <AMReX_REAL.H>

#define MAX_CHEM_REAL_VAR 29
#define MAX_CHEM_INT_VAR 16

using amrex::Real;

/**
 * Note: if the enumeration starts at 0, then count will represent the total
 * number of values defined prior to the count variable.
 */

/**
 * Generic set of indices for variables characterizing a cell. Specific real
 * properties of a cell are enumerated first, followed by the first index of an
 * uncharacterized block of real variables that can be fed to other routines.
 * The indices below assume that the most generic cell is an ellipsoid with axes
 * a, b, c oriented in space by three Euler angles and acted upon by a linear
 * force on the center of mass and a torque. The cell is also characterized by a
 * volume V. The cell has linear velocity and angular velocity and the cell
 * volume is also growing or shrinking with some rate. For spherical cells,
 * a=b=c=radius and the Euler angles, rotation velocity and torques are all
 * zero.
 */
struct realIdx
{
  enum { // Locations of real particle data stored in particle type
    radius,      // radius of cylinder
    c_length,    // length of cylinder
    theta,       // polar angle
    phi,         // azimuthal angle
    area,        // segment surface area                                 5
    vol,         // segment volume
    velx,        // x-component of linear velocity
    vely,        // y-component of linear velocity
    velz,        // z-component of linear velocity
    wx,          // x-component of angular velocity                     10
    wy,          // y-component of angular velocity
    wz,          // z-component of angular velocity
    fx,          // x-component of linear force
    fy,          // y-component of linear force
    fz,          // z-component of linear force                         15
    taux,        // x-component of torque
    tauy,        // y-component of torque
    tauz,        // z-component of torque
    dadt,        // rate of change of the segment area
    dvdt,        // rate of change of the segment volume                20
    first_data,  // index of first value in uncharacterized block 22
    count        // total number of elements in this enumeration
  };
};

/**
 * Generic set of indices for integer variables characterizing a cell. Specific
 * integers are listed first followed by an uncharacterized block of integers
 * that can be fed to other routines. The total number of real and integer
 * variables in this block and the readIdx block are stored in this struct.
 */
struct intIdx
{
  enum { // Locations of integer data stored in particle type
    num_reals,      // Number unique reals in uncharacterized block of real data
    num_ints,       // Number unique integers in uncharacterized block of intger data
    real_tot,       // Total number reals (unique plus temporary values)
    int_tot,        // Total number of integers (unique plus temporary values)
    first_real_inc, // Location of first value of concentration increments        5
    cell_type,      // Enumeration to distinguish between different cell types
    n_bnds,         // Number of segments bonded to this segment
    seg1_id1,       // unique identifier of bonded segment 1
    seg1_id2,       // original process of bonded segment 1
    seg2_id1,       // unique identifier of bonded segment 2                     10
    seg2_id2,       // original process of bonded segment 2
    seg3_id1,       // unique identifier of bonded segment 3
    seg3_id2,       // original process of bonded segment 3
    site1,          // bonding site of segment 1
    site2,          // bonding site of segment 2                                 15
    site3,          // bonding site of segment 3                                 16
    first_data,     // index of first value of uncharacterized block
    count           // total number of elements in this enumeration
  };
};

struct cellType
{
  enum { // Different cell types in simulation
    YEAST = 0, // Spherical, unconnected cells
    FUNGI      // Segment of filamentous fungi
  };
};

class BMXChemistry
{
  public:

  /**
   * Retrieve instance of BMXChemistry object
   */
  static BMXChemistry* instance();

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

#if 0
  /**
   * Transfer mesh values to internal concentrations and evaluate chemistry
   * @param grid_vol volume of grid cell that contains biological cell
   * @param cell_par pointer to cell parameter values
   * @param mesh_vals values of concentrations on mesh
   * @param p_vals values of concentrations in particles
   * @param dt time step interval
   */
  AMREX_GPU_HOST_DEVICE AMREX_INLINE
  void xferMeshToParticleAndUpdateChem(Real grid_vol, Real *cell_par,
      Real *mesh_vals, Real *p_vals, Real dt);
#endif

  /**
   * Print concentrations of chemical species in cell
   * @param id particle index
   * @param p_vals values of concentrations in particles
   * @param p_par values of particle parameters
   */
  void printCellConcentrations(int id, Real *p_vals, Real *p_par);

  /**
   * Return values that are also stored in particle idata fields
   * @param idx intIdx parameter representing index of desired value
   * @return integer value at index location
   */
  int getIntData(int idx);

  /**
   * Return a vector containing all chemical paramters needed to integrate
   * chemistry
   * @param chempar vector of parameters used in chemistry
   */
  void getChemParams(std::vector<Real> &chempar);

  // Hardwired variables for ABC model
  static Real k1, kr1;   // rate coefficients for species A into and out of biological cell
  static Real k2, kr2;   // rate coefficients for converting A into B and C inside biological cell
  static Real k3, kr3;   // rate coefficients for species C into and out of biological cell
  static Real kg;        // cell growth coefficient

  // Fraction of overlap
  static Real p_overlap;

  static int p_num_reals;   // Total number of real values allocated in model
  static int p_num_ints;    // Total number of int values allocated in model

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
  int p_inc_offset;  // Offset of increments relative to start of species data

  std::vector<int> p_int_vals; // stored values of integers

  // Flag for verbose output
  bool p_verbose;

  static BMXChemistry *p_instance;
};
#endif
