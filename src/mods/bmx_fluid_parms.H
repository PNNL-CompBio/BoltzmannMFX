#ifndef BMX_FLUID_PARMS_H_
#define BMX_FLUID_PARMS_H_

#include <AMReX_REAL.H>

#include <bmx_species_parms.H>

namespace FLUID
{
  // Flag to solve fluid equations
  extern int solve;

  // Flag to solve species fluid equations
  extern int solve_species;

  // Fluid phase species names
  extern amrex::Vector<std::string> species;

  // Species unique identifying code
  extern std::vector<int> species_id;

  // Total number of fluid species
  extern int nspecies;

  // Specified constant gas phase species diffusion coefficients
  extern std::vector<amrex::Real> D_gk0;

  // Name to later reference when building input regions
  extern std::string name;

  struct FLUID_t {

    amrex::Vector<SPECIES::SPECIES_t> species;

    FLUID_t() :
      species(0)
    {}

  };

  void Initialize ();

}

#endif
