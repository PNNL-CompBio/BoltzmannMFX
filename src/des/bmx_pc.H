#ifndef _BMX_PC_H_
#define _BMX_PC_H_

#include <AMReX.H>
#include <AMReX_AmrParGDB.H>
#include <AMReX_IArrayBox.H>
#include <AMReX_Particles.H>
#include <AMReX_NeighborParticles.H>

#include <AMReX_FabArray.H>
#include <AMReX_FabArray.H>
#include <AMReX_IArrayBox.H>
#include <AMReX_Particles.H>
#include <AMReX_NeighborParticles.H>
#include <bmx_bc_list.H>
#ifdef NEW_CHEM
#define MAX_CHEM_REAL_VAR 31
#define MAX_CHEM_INT_VAR 5
#include <bmx_chem.H>
#endif

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
    a_size,      // first ellipsoid dimension
    b_size,      // second ellipsoid dimension
    c_size,      // third ellipsoid dimension
    psi,         // first Euler angle
    theta,       // second Euler angle
    phi,         // third Euler angle
    area,        // cell surface area
    vol,         // cell volume
    velx,        // x-component of linear velocity
    vely,        // y-component of linear velocity
    velz,        // z-component of linear velocity
    wx,          // x-component of angular velocity
    wy,          // y-component of angular velocity
    wz,          // z-component of angular velocity
    fx,          // x-component of linear force
    fy,          // y-component of linear force
    fz,          // z-component of linear force
    taux,        // x-component of torque
    tauy,        // y-component of torque
    tauz,        // z-component of torque
    dadt,        // rate of change of the cell area
    dvdt,        // rate of change of the cell volume
    first_data,  // index of first value in uncharacterized block
    count        // total number of elements in this enumeration
  };
};

struct realData
{
  enum { // Real particle data stored in p.m_data
    velx,
    vely,
    velz,
    radius,
    volume,
    fluid_A,
    fluid_B,
    consume_A,
    consume_B,
    count
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
    first_real_inc, // Location of first value of concentration increments
    first_data,     // index of first value of uncharacterized block
    count           // total number of elements in this enumeration
  };
};

struct intData
{
  enum { // Integer particle data stored in p.m_data
    state,
    phase,
    count
  };
};

///
/// We add a particle to the neighbor list if its center is within 3 particle radii
///
struct BMXCheckPair
{
  amrex::Real neighborhood;
  BMXCheckPair (amrex::Real neighborhood_) : neighborhood(neighborhood_) {}

  template <class P>
  AMREX_GPU_DEVICE AMREX_INLINE
  bool operator() (const P* particles, const int i, const int j) const
  {
    const P& p1 = particles[i];
    const P& p2 = particles[j];

    if (i < j)
      if ((p1.pos() - p2.pos()).radSquared() <= neighborhood)
        return true;  

    return false;
  }
};

// Cannot inherit from this class?
#ifdef NEW_CHEM
class BMXParticleContainer final
    : public amrex::NeighborParticleContainer<MAX_CHEM_REAL_VAR,MAX_CHEM_INT_VAR>
#else
class BMXParticleContainer final
    : public amrex::NeighborParticleContainer<realData::count,intData::count>
#endif
{
public:
#ifdef NEW_CHEM
    using BMXParIter = amrex::ParIter<MAX_CHEM_REAL_VAR,MAX_CHEM_INT_VAR,0,0>;
    using BMXParConstIter = amrex::ParConstIter<MAX_CHEM_REAL_VAR,MAX_CHEM_INT_VAR,0,0>;
#else
    using BMXParIter = amrex::ParIter<realData::count,intData::count,0,0>;
    using BMXParConstIter = amrex::ParConstIter<realData::count,intData::count,0,0>;
#endif

    // Get number of particles
    int NumberOfParticles (BMXParIter& pti)
    { return pti.numParticles(); }

    BMXParticleContainer (amrex::AmrCore* amr_core);

    virtual ~BMXParticleContainer () {};

    void InitData ();

    void InitParticlesAscii (const std::string& file);

    void printParticles ();

    void AllocData ();

    void PrintParticleCounts ();

    void EvolveParticles (int lev,
                          int nstep,
                          amrex::Real dt,
                          amrex::Real time,
                          amrex::MultiFab * cost,
                          std::string & knapsack_weight_type,
                          int & nsubsteps);

    ///
    /// Write out all particles in the format of the particle_input.dat files so we
    ///       can initialize a run with them
    ///
    void WriteAsciiFileForInit (const std::string& filename);

    ///
    /// Write out all particles at a level
    ///
    void writeAllAtLevel (int lev = 0);

    void SolidsVolumeDeposition (int lev,
                                 amrex::MultiFab& mf);

    template <typename FB>
    void SolidsVolumeDeposition (FB WeightFunc,
                                 int lev,
                                 amrex::MultiFab& mf);

    void InterphaseTxfrDeposition (int lev,
                                   amrex::MultiFab& txfr_mf, amrex::Real dt);
 
    template <typename F>
    void InterphaseTxfrDeposition (F WeightFunc,
                                   int lev,
                                   amrex::MultiFab& txfr_mf, amrex::Real dt);

  private:
    int nlev;

    amrex::IntVect m_sorting_bin = amrex::IntVect::TheZeroVector();

    static int domain_bc[6];

    static void ReadStaticParameters ();
};
#endif
