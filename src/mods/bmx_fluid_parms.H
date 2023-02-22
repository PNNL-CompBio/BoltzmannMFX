//
//     Copyright (c) 2013 Battelle Memorial Institute
//     Licensed under modified BSD License. A copy of this license can be found
//     in the LICENSE file in the top level directory of this distribution.
//
#ifndef BMX_FLUID_PARMS_H_
#define BMX_FLUID_PARMS_H_

#include <AMReX_REAL.H>

#include <bmx_chem_species_parms.H>

namespace FLUID
{
  // Flag to solve fluid equations
  extern int solve;

  // Flag to solve chem_species fluid equations
  extern int solve_chem_species;

  // Fluid phase chem_species names
  extern amrex::Vector<std::string> chem_species;

  // ChemSpecies unique identifying code
  extern std::vector<int> chem_species_id;

  // Total number of fluid chem_species
  extern int nchem_species;

#ifdef NEW_CHEM
  // Initial fluid concentrations of species
  extern std::vector<amrex::Real> init_conc;

  // Location of top of support layer
  extern amrex::Real surface_location;
#endif

  // Specified constant gas phase chem_species diffusion coefficients
  extern std::vector<amrex::Real> D_k0;

  // Name to later reference when building input regions
  extern std::string name;

  struct FLUID_t {

    amrex::Vector<SPECIES::SPECIES_t> chem_species;

    FLUID_t() :
      chem_species(0)
    {}

  };

  void Initialize ();

}

#endif
