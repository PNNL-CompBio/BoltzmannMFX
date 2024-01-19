//
//     Copyright (c) 2013 Battelle Memorial Institute
//     Licensed under modified BSD License. A copy of this license can be found
//     in the LICENSE file in the top level directory of this distribution.
//
#include <bmx_pc.H>
#include <bmx_dem_parms.H>
#include <bmx_chem_species_parms.H>
#include <bmx_fluid_parms.H>
#ifdef NEW_CHEM
#include <bmx_chem.H>
#endif
// #include <bmx_ic_parms.H>

using namespace amrex;

void BMXParticleContainer::InitParticlesAscii (const std::string& file)
{
  // only read the file on the IO proc
  if (ParallelDescriptor::IOProcessor())
  {
    std::ifstream ifs;
    ifs.open(file.c_str(), std::ios::in);

#ifdef NEW_CHEM
    BMXChemistry *bmxchem = BMXChemistry::instance();
#endif

    if (!ifs.good())
      amrex::FileOpenFailed(file);

    // read in number of particles from input file
    int np = -1;
    ifs >> np >> std::ws;

    amrex::Print() << "Now reading " << np << " particles from " << file << std::endl;

    // Issue an error if nparticles = 0 is specified
    if ( np == -1 ){
      Abort("\nCannot read number of particles from particle_input.dat: file is corrupt.\
                   \nPerhaps you forgot to specify the number of particles on the first line??? ");
    }

    // we add all the particles to grid 0 and tile 0 and let
    // Redistribute() put them in the right places.
    const int lev  = 0;
    const int grid = 0;
    const int tile = 0;

    auto& particles = DefineAndReturnParticleTile(lev,grid,tile);
    particles.resize(np);

    Gpu::HostVector<ParticleType> host_particles(np);

    std::vector<amrex::Real> init_conc = FLUID::init_conc;

    for (int i = 0; i < np; i++)
    {
      // Read from input file
      ifs >> host_particles[i].pos(0);
      ifs >> host_particles[i].pos(1);
      ifs >> host_particles[i].pos(2);
      ifs >> host_particles[i].rdata(realIdx::radius);
      ifs >> host_particles[i].rdata(realIdx::c_length);   // 5
      ifs >> host_particles[i].rdata(realIdx::theta);
      ifs >> host_particles[i].rdata(realIdx::phi);
      ifs >> host_particles[i].rdata(realIdx::area);
      ifs >> host_particles[i].rdata(realIdx::vol);
      ifs >> host_particles[i].rdata(realIdx::velx);       // 10
      ifs >> host_particles[i].rdata(realIdx::vely);
      ifs >> host_particles[i].rdata(realIdx::velz);
      ifs >> host_particles[i].rdata(realIdx::wx);
      ifs >> host_particles[i].rdata(realIdx::wy);
      ifs >> host_particles[i].rdata(realIdx::wz);         // 15
      ifs >> host_particles[i].rdata(realIdx::fx);
      ifs >> host_particles[i].rdata(realIdx::fy);
      ifs >> host_particles[i].rdata(realIdx::fz);
      ifs >> host_particles[i].rdata(realIdx::taux);
      ifs >> host_particles[i].rdata(realIdx::tauy);       // 20
      ifs >> host_particles[i].rdata(realIdx::tauz);
      ifs >> host_particles[i].rdata(realIdx::gx);
      ifs >> host_particles[i].rdata(realIdx::gy);
      ifs >> host_particles[i].rdata(realIdx::gz);
      ifs >> host_particles[i].rdata(realIdx::dadt);       // 25
      ifs >> host_particles[i].rdata(realIdx::dvdt);
      ifs >> host_particles[i].idata(intIdx::cell_type);   // 27
      for (int c=0; c<FLUID::nchem_species; c++) 
      {
       host_particles[i].rdata(realIdx::first_data+c) = init_conc[c];
        //host_particles[i].rdata(realIdx::first_data+c) = tmp[c];
      }
      bmxchem->setIntegers(&host_particles[i].idata(0));
      host_particles[i].rdata(realIdx::bond_scale) = 1.0;
      host_particles[i].idata(intIdx::fuse_tip) = 0;
      
      // Set id and cpu for this particle
      host_particles[i].id()  = ParticleType::NextID();
      host_particles[i].cpu() = ParallelDescriptor::MyProc();
      host_particles[i].idata(intIdx::id) = static_cast<int>(host_particles[i].id());
      host_particles[i].idata(intIdx::cpu) = static_cast<int>(host_particles[i].cpu());
      host_particles[i].idata(intIdx::position) = siteLocation::TIP;
      host_particles[i].idata(intIdx::n_bnds) = 0;
      host_particles[i].idata(intIdx::fuse_flag) = 0;

      if (!ifs.good())
          amrex::Abort("Error initializing particles from Ascii file. \n");
    }
    // For testing inter-segment exchange
    // host_particles[3].rdata(realIdx::first_data+1) = 1.0e-5;
    for (int i = 0; i < np; i++)
    {
      printf("PARTICLE: %d:%d x: %f y: %f z: %f theta: %f phi: %f [A]: %e [B]: %e [C]: %e type: %d\n",
          (int)host_particles[i].id(),
          (int)host_particles[i].cpu(),
          host_particles[i].pos(0),
          host_particles[i].pos(1),
          host_particles[i].pos(2),
          host_particles[i].rdata(realIdx::theta),
          host_particles[i].rdata(realIdx::phi),
          host_particles[i].rdata(realIdx::first_data),
          host_particles[i].rdata(realIdx::first_data+1),
          host_particles[i].rdata(realIdx::first_data+2),
          host_particles[i].idata(intIdx::cell_type)
          );
    }

    auto& aos = particles.GetArrayOfStructs();
    Gpu::DeviceVector<ParticleType>& gpu_particles = aos();

    // Copy particles from host to device
    Gpu::copyAsync(Gpu::hostToDevice, host_particles.begin(), host_particles.end(), gpu_particles.begin());
  }

  Redistribute();
}

