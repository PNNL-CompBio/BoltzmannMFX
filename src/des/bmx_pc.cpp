//
//     Copyright (c) 2013 Battelle Memorial Institute
//     Licensed under modified BSD License. A copy of this license can be found
//     in the LICENSE file in the top level directory of this distribution.
//
#include <bmx_des_K.H>

#include <bmx_dem_parms.H>
#include <bmx_bc_parms.H>
#include <bmx_pc.H>

using namespace amrex;

int  BMXParticleContainer::domain_bc[6] {0};

#ifdef NEW_CHEM
BMXParticleContainer::BMXParticleContainer (AmrCore* amr_core)
    : NeighborParticleContainer<MAX_CHEM_REAL_VAR,MAX_CHEM_INT_VAR>
#else
BMXParticleContainer::BMXParticleContainer (AmrCore* amr_core)
    : NeighborParticleContainer<realData::count,intData::count>
#endif
      (amr_core->GetParGDB(), 1)
{
    ReadStaticParameters();

    {
      this->SetVerbose(0);
      ParmParse pp("bmx");
      int verbose = 0; 
      pp.query("verbose",verbose);
      if (verbose != 0) {
        p_verbose = true;
      } else {
        p_verbose = false;
      }
    }

    nlev         = amr_core->finestLevel()+1;
    finest_level = amr_core->finestLevel();
}

void BMXParticleContainer::AllocData ()
{
    reserveData();
    resizeData();
}

void BMXParticleContainer::PrintParticleCounts ()
{
  const int lev = 0;
  amrex::AllPrintToFile("load_balance") << "Particles on each box: \n";
  long local_count = 0;
  for (BMXParIter pti(*this, lev); pti.isValid(); ++pti)
    {
      long np = pti.numParticles();
      local_count += np;
      amrex::AllPrintToFile("load_balance") << "Box:" << pti.index() << ", count: " << np << std::endl;
    }
  amrex::AllPrintToFile("load_balance") << "Total for this process: " << local_count << std::endl << std::endl;
}

void BMXParticleContainer::printParticles ()
{
    const int lev = 0;
    auto& plevel = GetParticles(lev);

    for (auto& kv : plevel)
    {
       const auto& particles = kv.second.GetArrayOfStructs();
       auto& soa = kv.second.GetStructOfArrays();
       auto p_realarray = soa.realarray();

       for (int i = 0; i < particles.numParticles(); ++i)
       {
          std::cout << "Particle ID  = " << i << " " << std::endl;
          std::cout << "X            = " << particles[i].pos(0) << " " << std::endl;
          std::cout << "Y            = " << particles[i].pos(1) << " " << std::endl;
          std::cout << "Z            = " << particles[i].pos(2) << " " << std::endl;
#ifdef NEW_CHEM
          for (int j = realIdx::count-1; j < realIdx::count + particles[i].idata(intIdx::num_reals)-1; j++)
#else
          std::cout << "phase        = " << p_intarray[intData::phase][i] << " " << std::endl;
          std::cout << "Real properties = " << std::endl;

          for (int j = 0; j < realData::count; j++)
#endif
            std::cout << "property " << j << "  = " << p_realarray[j][i] << " " << std::endl;

          std::cout << std::endl;
       }
    }
}

Real 
BMXParticleContainer::computeParticleVolume () const
{
    auto r = amrex::ReduceSum(*this, [=]
       AMREX_GPU_HOST_DEVICE (const ParticleType& p) noexcept -> amrex::Real
       {
           return p.rdata(realIdx::vol);
       });
    ParallelDescriptor::ReduceRealSum(r);

    return r;
}

Real 
BMXParticleContainer::computeParticleContent (int comp) const
{
    auto r = amrex::ReduceSum(*this, [=]
       AMREX_GPU_HOST_DEVICE (const ParticleType& p) noexcept -> amrex::Real
       {
           return (p.rdata(realIdx::vol) * p.rdata(comp));
       });

    ParallelDescriptor::ReduceRealSum(r);
    return r;
}

void BMXParticleContainer::ReadStaticParameters ()
{
    static bool initialized = false;

    if (!initialized)
        initialized = true;
}

void BMXParticleContainer::writeAllAtLevel (int lev)
{
    // Not threaded because its print to terminal
    for (BMXParIter pti(*this, lev); pti.isValid(); ++pti)
    {
        auto& particles = pti.GetArrayOfStructs();
        int np = pti.numParticles();
        Gpu::HostVector<ParticleType> host_particles(np);
        Gpu::copy(Gpu::deviceToHost, particles.begin(), particles.end(), host_particles.begin());

        for (const auto& p: host_particles)
        {
           const IntVect& iv = Index(p, lev);

           RealVect xyz(p.pos(0), p.pos(1), p.pos(2));
           std::cout << " id " << p.id()
                << " index " << iv
                << " position " << xyz << std::endl;
       }
    }
}

