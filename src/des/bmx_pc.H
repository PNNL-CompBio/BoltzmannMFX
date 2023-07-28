//
//     Copyright (c) 2013 Battelle Memorial Institute
//     Licensed under modified BSD License. A copy of this license can be found
//     in the LICENSE file in the top level directory of this distribution.
//
#ifndef _BMX_PC_H_
#define _BMX_PC_H_

#include <AMReX_AmrParGDB.H>
#include <AMReX_Particles.H>
#include <AMReX_NeighborParticles.H>

#include <bmx_bc_list.H>
#ifdef NEW_CHEM
#include <bmx_chem.H>
#endif

/**
 * Note: if the enumeration starts at 0, then count will represent the total
 * number of values defined prior to the count variable.
 */

#ifndef NEW_CHEM
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
#endif

#ifndef NEW_CHEM
struct intData
{
  enum { // Integer particle data stored in p.m_data
    state,
    phase,
    count
  };
};
#endif

///
/// We add a particle to the neighbor list if its center is within 3 particle radii
///
struct BMXCheckPair
{
  amrex::Real neighborhood;
  bool single_pass;
  BMXCheckPair (amrex::Real neighborhood_) : neighborhood(neighborhood_), single_pass(true) {}
  BMXCheckPair (amrex::Real neighborhood_, bool single_pass_) : neighborhood(neighborhood_),
      single_pass(single_pass_) {}

  template <class P>
  AMREX_GPU_DEVICE AMREX_INLINE
  bool operator() (const P* particles, const int i, const int j) const
  {
    const P& p1 = particles[i];
    const P& p2 = particles[j];

    if (single_pass) {
      if (i < j) {
        if ((p1.pos() - p2.pos()).radSquared() <= neighborhood)
          return true;  
      }
    } else { 
      if ((p1.pos() - p2.pos()).radSquared() <= neighborhood)
        return true;  
    }

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

    void WriteToAscii (const std::string& file,int nstep,amrex::Real time);
    void countParticleTypes(int &nparticles, int &nsegments);
    void WriteSegmentsToAscii (const std::string& file,int nstep,
        amrex::Real time, int nsegments);
    void WriteParticlesToAscii (const std::string& file,int nstep,
        amrex::Real time, int nparticles);

    void printParticles ();

    Real computeParticleVolume () const;

    Real computeParticleContent (int comp) const;

    void AllocData ();

    void PrintParticleCounts ();

    void EvolveParticles (amrex::Real dt,
                          const amrex::Vector<amrex::MultiFab *> cost,
                          std::string & knapsack_weight_type,
                          int & nsubsteps);
    void ParticleExchange (Real dt,
                           const amrex::Vector<amrex::MultiFab*> cost,
                           std::string& knapsack_weight_type,
                           int& nsubsteps);

    bool EvaluateTipFusion(const amrex::Vector<amrex::MultiFab*> cost,
                                   std::string& knapsack_weight_type);

    void EvaluateInteriorFusion(const amrex::Vector<amrex::MultiFab*> cost,
                                   std::string& knapsack_weight_type);

    void CleanupFusion(const amrex::Vector<amrex::MultiFab*> cost,
                                   std::string& knapsack_weight_type);

    void InitBonds(const amrex::Vector<amrex::MultiFab*> cost,
                                   std::string& knapsack_weight_type);

    void PrintConnectivity(const amrex::Vector<amrex::MultiFab*> cost,
                                   std::string& knapsack_weight_type);

    ///
    /// Write out all particles in the format of the particle_input.dat files so we
    ///       can initialize a run with them
    ///
    void WriteAsciiFileForInit (const std::string& filename);

    ///
    /// Write out all particles at a level
    ///
    void writeAllAtLevel (int lev = 0);

    void split_particles (Real time);

  private:
    int nlev;
    int finest_level;

    amrex::IntVect m_sorting_bin = amrex::IntVect::TheZeroVector();

    static int domain_bc[6];

    static void ReadStaticParameters ();

    bool p_verbose;
};
#endif
