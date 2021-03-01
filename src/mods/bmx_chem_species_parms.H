#ifndef BMX_SPECIES_PARMS_H_
#define BMX_SPECIES_PARMS_H_

#include <AMReX_REAL.H>

namespace SPECIES
{
  constexpr int NMAX = 15;

  // List of recognized chem_species diffusivity models
  struct DIFFUSIVITYMODEL { enum {Invalid=-1, Constant}; };

  extern int DiffusivityModel;

  // List of recognized chem_species specific heat models
  struct SPECIFICHEATMODEL { enum {Invalid=-1, Constant}; };

  extern int SpecificHeatModel;

  // Flag to solve chem_species equations
  extern int solve;

  // Total number of chem_species
  extern int nchem_species;

  // Specified chem_species names
  extern std::vector<std::string> chem_species;

  // ChemSpecies unique identifying code (at the moment = their index in the input
  // entries)
  extern std::vector<int> chem_species_id;

  // Specified chem_species constant diffusion coefficients
  extern std::vector<amrex::Real> D_k0;

  struct SPECIES_t {
    // Constructor
    SPECIES_t() :
      mass_fraction(0.)
    {}

    // ChemSpecies mass fraction, used for ICs and BCs
    amrex::Real mass_fraction;
  };

  void Initialize ();

}

#endif
