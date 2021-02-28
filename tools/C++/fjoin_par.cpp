//  Author: Jordan Musser
//  Date  : May, 2020
//
#include <iostream>
#include <cstdlib>
#include <stdlib.h>
#include <string.h>
#include <sstream>
#include <iomanip>
#include <AMReX.H>
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>
#include <iostream>

#include <iostream>     // std::cout
#include <fstream>      // std::ifstream

using namespace amrex;
using namespace std;

//
// Prototypes
//
void get_input_arguments ( const int argc, char** argv,
        std::string& a_fbase, std::string& a_fjoin,
        int& a_id, amrex::Vector<int>& a_var,
        int& a_start, int& a_end, amrex::Real& a_dt,
        int& a_format, bool& a_verbose);

void help ();

struct particle_t {
  amrex::Vector<int>         idata;
  amrex::Vector<amrex::Real> rdata;
};


amrex::Real calc_granular_temperature (amrex::Vector<particle_t> a_particles);
amrex::Real cleaned_value (amrex::Real value_in);



int main ( int argc, char* argv[] )
{
  std::string  fbase;
  std::string  fjoin;
  int id(-1);
  int istart(0);
  int iend(0);
  int iformat(10);
  amrex::Real dt(0.0);
  amrex::Vector<int> var;
  bool verbose(false);

  get_input_arguments(argc, argv, fbase, fjoin, id, var, istart, iend, dt,
      iformat, verbose);

  const size_t var_count = var.size();

  if(verbose)
  {
    // Print summary
    std::cout << "\n\nfjoin_par input summary: **"
              << "\n   File base name " << fbase
              << "\n   Output File    " << fjoin
              << "\n   Verbose        " << verbose
              << "\n   Particle ID    " << id
              << "\n   Simulation dt  " << dt
              << "\n   Interval start " << istart
              << "\n   Interval end   " << iend
              << "\n   Format         " << iformat
              << "\n   Variable count " << var_count
              << "\n\n";
  }

  int fcount(0);
  int lc2(0);
  int err(0);

  int npo(-1), nro(-1), nio(-1);

  ofstream  outfile;
  outfile.open( fjoin.c_str(), ofstream::out | ofstream::trunc );

  while( lc2 + istart <= iend && err == 0 )
  {
    std::stringstream clc2;
    clc2 << std::setw(5) << std::setfill('0') << istart + lc2;
    std::string lfile = fbase + clc2.str();

    std::ifstream ifs;
    ifs.open(lfile.c_str(), std::ios::in);

    if (ifs.good()) {
      //if(verbose) std::cout << lfile << "     found!" << std::endl;

      int np(-1), AoS_nr(-1), AoS_ni(-1), SoA_nr(-1), SoA_ni(-1);

      ifs >> np;
      ifs >> AoS_nr;
      ifs >> AoS_ni;
      ifs >> SoA_nr;
      ifs >> SoA_ni;

      AMREX_ALWAYS_ASSERT_WITH_MESSAGE((AoS_nr + SoA_nr) == 19,
          "Number of reals does not equal 20. Need to update fjoin par.");

      AMREX_ALWAYS_ASSERT_WITH_MESSAGE((AoS_ni + SoA_ni) ==  2,
          "Number of ints does not equal 2. Need to update fjoin par.");

      if(lc2 == 0){
        npo = np;
        nro = AoS_nr + SoA_nr;
        nio = AoS_ni + SoA_ni;
      } else {
        err += (np == npo) ? 0 : 1;
        err += ((AoS_nr+SoA_nr) == nro) ? 0 : 1;
        err += ((AoS_ni+SoA_ni) == nio) ? 0 : 1;
      }

      if( dt > 0.0)
      {
        outfile << fixed << uppercase << setprecision(6);
        outfile << setw(18) << dt*lc2;
      }

      // Read and store particle data.
      amrex::Vector<particle_t> particles;

      for (int lc(0); lc < np; lc++)
      {
        particle_t p;

        amrex::Real rtmp;
        int itmp;

        ifs >> rtmp;  p.rdata.push_back(rtmp);  //   : x-position (AoS)
        ifs >> rtmp;  p.rdata.push_back(rtmp);  //   : y-position (AoS)
        ifs >> rtmp;  p.rdata.push_back(rtmp);  //   : z-position (AoS)

        for (int AoS_rcount(0); AoS_rcount < AoS_nr; ++AoS_rcount) {
          ifs >> rtmp;  p.rdata.push_back(rtmp);  // : AoS Real variables
        }

        ifs >> itmp;  p.idata.push_back(itmp);  //   : id  (AoS)
        ifs >> itmp;  p.idata.push_back(itmp);  //   : cpu (AoS)

        for (int AoS_icount(0); AoS_icount < AoS_ni; ++AoS_icount) {
          ifs >> itmp;  p.idata.push_back(itmp);  // : AoS int variables
        }

        for (int SoA_rcount(0); SoA_rcount < SoA_nr; ++SoA_rcount) {
          ifs >> rtmp;  p.rdata.push_back(rtmp);  // : SoA Real variables
        }

        for (int SoA_icount(0); SoA_icount < SoA_ni; ++SoA_icount) {
          ifs >> itmp;  p.idata.push_back(itmp);  // : SoA int variables
        }

        particles.push_back(p);

        if (!ifs.good()){
          amrex::Abort("Error reading particles from Ascii file. \n");
        }
      }

      outfile << fixed << uppercase << setprecision(iformat);

      for(size_t lc=0; lc < var_count; lc++)
      {
        const bool last = (lc == var_count-1);
        const int var_id = var[lc];

        if( var_id == 100 )
        {
          outfile << setw(24) << cleaned_value(calc_granular_temperature(particles));
        }
        else 
        {
          if( id == -1 )
          {
            if(verbose)
              std::cout << "Writing all particles!" << std::endl;

            for(int llc(0); llc<np; llc++)
            {
              outfile << setw(24) << cleaned_value(particles[llc].rdata[var_id-1]);
            }
          }
          else
          {
            outfile << setw(24) << cleaned_value(particles[id-1].rdata[var_id-1]);
          }
        }

        if(last) 
          outfile << std::endl;
      }

    }
    else {
      if(verbose)
        std::cout << lfile << " not found!" << std::endl;
    }

    lc2++;
  }

  outfile.close();
  return err;
};


//
// Parse command line arguments
//
void get_input_arguments ( const int argc, char** argv,
        std::string& a_fbase, std::string& a_fjoin,
        int& a_id, amrex::Vector<int>& a_var,
        int& a_start, int& a_end, amrex::Real& a_dt,
        int& a_format, bool& a_verbose )
{

    int i = 1; // skip program name

    while ( i < argc )
    {

      if ( !strcmp(argv[i],"-f") || !strcmp(argv[i],"--file") ) {
          a_fbase = argv[++i];

      } else if ( !strcmp(argv[i],"-j") || !strcmp(argv[i],"--join") ) {
          a_fjoin = argv[++i];

      } else if ( !strcmp(argv[i],"--verbose") ) {
          a_verbose = true;

      } else if ( !strcmp(argv[i],"--id") ) {
          a_id = atoi(argv[++i]);
          if ( ( a_id < -1 ) or ( a_id == 0 ) ) {
            std::cout << "\nParticle ID must be positive! " << a_id << std::endl;
            help ();
            exit (EXIT_FAILURE);
          }

      } else if ( !strcmp(argv[i],"--format") ) {
          a_format = atoi(argv[++i]);

      } else if ( !strcmp(argv[i],"--var") ) {
          a_var.push_back(atoi(argv[++i]));

      } else if ( !strcmp(argv[i],"--start") ) {
          a_start = atoi(argv[++i]);

      } else if ( !strcmp(argv[i],"--end") ) {
          a_end = atoi(argv[++i]);

      } else if ( !strcmp(argv[i],"--dt") ) {
          a_dt = atof(argv[++i]);

      } else {
          std::cout << "\n\nOption " << argv[i] << " not recognized" << std::endl;
          help ();
          exit ( EXIT_FAILURE );
      }

      // Go to the next parameter name
      ++i;
    }

    if ( a_fbase.empty () ) {
      std::cout << "\n\nRequired option [--file] is missing" << std::endl;
      help ();
      exit (EXIT_FAILURE);
    }

    if ( a_fjoin.empty () ) {
      a_fjoin = a_fbase + ".join";
    }

}



//
// Print usage info
//
void help ()
{
    std::cout << "\n\nProgram to compare particles data written in ASCII file  by AMReX."
        << "\n "
        << "\nUsage:"
        << "\n   fjoin_par [--file STRING] [--id INTEGER]"
        << "\n        [--var INTEGER] [--dt REAL] [--verbose]"
        << "\n "
        << "\nargs --file, -f  : The base name used to generate ASCII particle output"
        << "\n                   files. This is the same as `amr.par_ascii_file` in"
        << "\n                   the input deck."
        << "\n     --join, -j  : File name of output file"
        << "\n     --verbose   : Output format, Number of sig-figs to write {10 (default)} "
        << "\n     --id        : ID of particle for data extraction."
        << "\n     --dt        : This is the simulation dt separating the ascii output"
        << "\n                   files. When this is passed, the output tries to reconstruct"
        << "\n                   the simulation output time."
        << "\n     --start     : "
        << "\n     --end       : "
        << "\n     --var       : Index of particle property for extraction. This input"
        << "\n                   may be passed multiple times."
        << "\n                      1 - position-x"
        << "\n                      2 - position-y"
        << "\n                      3 - position-z"
        << "\n                      4 - radius"
        << "\n                      5 - volume"
        << "\n                      6 - mass"
        << "\n                      7 - density"
        << "\n                      8 - oneOverI"
        << "\n                      9 - velocity-x"
        << "\n                     10 - velocity-y"
        << "\n                     11 - velocity-z"
        << "\n                     12 - omega-x"
        << "\n                     13 - omega-y"
        << "\n                     14 - omega-z"
        << "\n                     15 - drag-x"
        << "\n                     16 - drag-y"
        << "\n                     17 - drag-z"
        << "\n "
        << "\n                    100 - kinetic energy"
        << "\n\n";

}


amrex::Real calc_granular_temperature (amrex::Vector<particle_t> a_particles)
{
  amrex::Real gtmp(0.0);
  amrex::Real np = a_particles.size();

  for(int lc(0); lc < np; lc++){
    gtmp += a_particles[lc].rdata[ 8]*a_particles[lc].rdata[ 8]
         +  a_particles[lc].rdata[ 9]*a_particles[lc].rdata[ 9]
         +  a_particles[lc].rdata[10]*a_particles[lc].rdata[10];
  }
  amrex::Real myval = gtmp / (3.0 *np);
  return gtmp / (3.0 * np);
}


amrex::Real cleaned_value (amrex::Real value_in)
{
  const Real tolerance = std::numeric_limits<Real>::epsilon();
  return (std::abs(value_in) > tolerance) ? value_in : 0.0;
}
