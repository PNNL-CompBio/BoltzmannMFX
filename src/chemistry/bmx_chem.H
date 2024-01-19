//
//     Copyright (c) 2013 Battelle Memorial Institute
//     Licensed under modified BSD License. A copy of this license can be found
//     in the LICENSE file in the top level directory of this distribution.
//
#ifndef _BMX_CHEM_H_
#define _BMX_CHEM_H_

#include <AMReX.H>
#include <AMReX_REAL.H>

/**
 * These two values need to be coordinated with the number of variables
 * defined in the realIdx and intIdx enumerations, plus the number of
 * variables in the undefined blocks */
#define NUM_CHEM_COMPONENTS 5
#define MAX_CHEM_REAL_VAR 28 + 3*NUM_CHEM_COMPONENTS
#define MAX_CHEM_INT_VAR 35

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
 * a=b=c=radius and the orientation angles, rotation velocity and torques are
 * all zero.
 */
struct realIdx
{
  enum { // Locations of real particle data stored in particle type
    radius,      // radius of particle or cylinder
    c_length,    // length of cylinder
    theta,       // polar angle
    phi,         // azimuthal angle
    area,        // surface area of sphere or segment                   5
    vol,         // volume of sphere or segment
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
    gx,          // x-component of chemical gradient
    gy,          // y-component of chemical gradient                    20
    gz,          // z-component of chemical gradient
    dadt,        // rate of change of the area
    dvdt,        // rate of change of the volume
    tau_split,   // value of tau to use if fusion occurs
    bond_scale,  // parameter to reduce value of bond force             25
    ranvx,       // random x velocity
    ranvy,       // random y velocity
    ranvz,       // random z velocity                                   28
    first_data,  // index of first value in uncharacterized block
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
    seg2_id1,       // unique identifier of bonded segment 2
    seg3_id1,       // unique identifier of bonded segment 3                     10
    seg4_id1,       // unique identifier of bonded segment 4
    seg1_id2,       // original process of bonded segment 1
    seg2_id2,       // original process of bonded segment 2
    seg3_id2,       // original process of bonded segment 3
    seg4_id2,       // original process of bonded segment 4                      15
    site1,          // bonding site of segment 1
    site2,          // bonding site of segment 2
    site3,          // bonding site of segment 3
    site4,          // bonding site of segment 4
    position,       // identify tip, second position, interior segment           20
    fuse_flag,      // flag that this particle is fusing to another particle
    split_flag,     // flag to identify splitting segment in a fusion event
    new_flag,       // flag to identify new particle in fusion event
    fuse_id,        // id of particle that this particle is fusing to
    fuse_cpu,       // cpu of particle that this particle is fusing to           25
    fuse_tip,       // site that fuses                                           
    dltd_site1,     // deleted site 1
    dltd_site2,     // deleted site 2
    dltd_id1,       // deleted id 1
    dltd_id2,       // deleted id 2                                              30
    dltd_cpu1,      // deleted cpu 1
    dltd_cpu2,      // deleted cpu 2
    fix_site,       // parameter to make sure that bonds fixed correctly
    id,             // unique identifier of segment
    cpu,            // original process of segment                               35
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

struct siteLocation
{
  enum { // Different locations of fungi segments
    TIP = 0,   // Growth tip filament
    SECOND_1,  // Segment behind growth tip, tip is bound to site 1
    SECOND_2,  // Segment behind growth tip, tip is bound to site 2
    INTERIOR   // 
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
  void getChemParams(amrex::Gpu::DeviceVector<Real> &chempar);

  /** Return a vector containing exchange parameters
   * @return vector of exchange parameters
   */
  amrex::Gpu::DeviceVector<Real> getExchangeParameters();

  /**
   * Return a vector containing parameters for segment fusion
   * @return vector of fusion parameters
   */
  amrex::Gpu::DeviceVector<Real> getFusionParameters();

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

  // Hardwired variables for ABC model
  static Real k1, kr1;   // rate coefficients for species A into and out of fungi
  static Real k2, kr2;   // rate coefficients for converting A into B and C inside fungi
  static Real k3, kr3;   // rate coefficients for species C into and out of fungi
  static Real k4, kr4;   // rate coefficients for species D into and out of bacterium
  static Real k5, kr5;   // rate coefficients for species C and D into E and A
  static Real k6, kr6;   // rate coefficients for species A into and out of bacterium
  static Real k7, kr7;   // rate coefficients for species C into and out of bacterium

  static Real kg, kv;        // fungi growth coefficients
  static Real kb, kbv;       // bacterium growth coefficients

  static Real bacteria_radius_max;  // maximum value of radius for bacterium
  static Real radius_max;  // maximum value of cylinder radius for fungi segment
  static Real length_max;  // maximum value of cylinder length for fungi segment
  static Real brnch_prob;  // probability of side branching
  static Real split_prob;  // probability of tip splitting

  static Real mtA;  // Septum mass transfer coefficient for A
  static Real mtB;  // Septum mass transfer coefficient for B
  static Real mtC;  // Septum mass transfer coefficient for C

  int p_num_species; // Total number chemical species in model
  int p_num_ivals;   // Total number integer settings in model
  int p_inc_offset;  // Offset of increments relative to start of species data

  static Real fusion_prob; // probability of segments fusing
  static Real max_fusion_separation; // maximum separtion when fusion can occur
  static Real scale_inc; // increment for increasing bond scaling parameter

  std::vector<int> p_int_vals; // stored values of integers

  // Flag for verbose output
  bool p_verbose;

  static BMXChemistry *p_instance;
};
#endif
