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

    amrex::Print() << "Now reading " << np << " particles from particle_input.dat" << std::endl;

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
#ifdef NEW_CHEM
      // Read from input file
      ifs >> host_particles[i].pos(0);
      ifs >> host_particles[i].pos(1);
      ifs >> host_particles[i].pos(2);
      ifs >> host_particles[i].rdata(realIdx::a_size);
      ifs >> host_particles[i].rdata(realIdx::b_size);
      ifs >> host_particles[i].rdata(realIdx::c_size);
      ifs >> host_particles[i].rdata(realIdx::psi);
      ifs >> host_particles[i].rdata(realIdx::theta);
      ifs >> host_particles[i].rdata(realIdx::phi);
      ifs >> host_particles[i].rdata(realIdx::area);
      ifs >> host_particles[i].rdata(realIdx::vol);
      ifs >> host_particles[i].rdata(realIdx::velx);
      ifs >> host_particles[i].rdata(realIdx::vely);
      ifs >> host_particles[i].rdata(realIdx::velz);
      ifs >> host_particles[i].rdata(realIdx::wx);
      ifs >> host_particles[i].rdata(realIdx::wy);
      ifs >> host_particles[i].rdata(realIdx::wz);
      ifs >> host_particles[i].rdata(realIdx::fx);
      ifs >> host_particles[i].rdata(realIdx::fy);
      ifs >> host_particles[i].rdata(realIdx::fz);
      ifs >> host_particles[i].rdata(realIdx::taux);
      ifs >> host_particles[i].rdata(realIdx::tauy);
      ifs >> host_particles[i].rdata(realIdx::tauz);
      ifs >> host_particles[i].rdata(realIdx::dadt);
      ifs >> host_particles[i].rdata(realIdx::dvdt);
      /*
      Real tmp[3];
      tmp[0] = 2.0e-05;
      tmp[1] = 0.0;
      tmp[2] = 2.0e-06;
      */
      for (int c=0; c<FLUID::nchem_species; c++) 
      {
       host_particles[i].rdata(realIdx::first_data+c) = init_conc[c];
        //host_particles[i].rdata(realIdx::first_data+c) = tmp[c];
      }
      bmxchem->setIntegers(&host_particles[i].idata(0));
      
#else
      ifs >> host_particles[i].rdata(realData::velx);
      ifs >> host_particles[i].rdata(realData::vely);
      ifs >> host_particles[i].rdata(realData::velz);
      ifs >> host_particles[i].rdata(realData::radius);
      ifs >> host_particles[i].rdata(realData::volume);
      ifs >> host_particles[i].idata(intData::phase);
      ifs >> host_particles[i].idata(intData::state);

      // These will hold the values interpolated from the fluid
      host_particles[i].rdata(realData::fluid_A) = 0.;
      host_particles[i].rdata(realData::fluid_B) = 0.;

      // These will hold the values that the particle is going to consume from the fluid -- 
      //  these will be deposited onto the grid to change X_A and X_B
      host_particles[i].rdata(realData::consume_A) = 0.;
      host_particles[i].rdata(realData::consume_B) = 0.;
#endif

      // Set id and cpu for this particle
      host_particles[i].id()  = ParticleType::NextID();
      host_particles[i].cpu() = ParallelDescriptor::MyProc();

      if (!ifs.good())
          amrex::Abort("Error initializing particles from Ascii file. \n");
    }

    auto& aos = particles.GetArrayOfStructs();
    Gpu::DeviceVector<ParticleType>& gpu_particles = aos();

    // Copy particles from host to device
    Gpu::copyAsync(Gpu::hostToDevice, host_particles.begin(), host_particles.end(), gpu_particles.begin());
  }

  Redistribute();
}

