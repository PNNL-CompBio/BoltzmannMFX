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

class BMXParticleContainer final
    : public amrex::NeighborParticleContainer<realData::count,intData::count>
{
public:
    using BMXParIter = amrex::ParIter<realData::count,intData::count,0,0>;
    using BMXParConstIter = amrex::ParConstIter<realData::count,intData::count,0,0>;

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
                                   amrex::MultiFab& txfr_mf);
 
    template <typename F>
    void InterphaseTxfrDeposition (F WeightFunc,
                                   int lev,
                                   amrex::MultiFab& txfr_mf);

  private:
    int nlev;

    amrex::IntVect m_sorting_bin = amrex::IntVect::TheZeroVector();

    static int domain_bc[6];

    static void ReadStaticParameters ();
};
#endif
