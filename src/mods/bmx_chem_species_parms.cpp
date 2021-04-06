#include <AMReX.H>
#include <AMReX_Arena.H>
#include <AMReX_Print.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Utility.H>
#include <AMReX_ParallelDescriptor.H>

#include <bmx_chem_species_parms.H>

namespace SPECIES
{
  // Flag to solve chem_species equations
  int solve;

  // Total number of chem_species
  int nchem_species(0);

  // ChemSpecies names
  std::vector<std::string> chem_species(0);

  // ChemSpecies unique identifying code (at the moment = their index in the input
  // entries)
  std::vector<int> chem_species_id(0);

  // Specified chem_species diffusion coefficients
  std::vector<amrex::Real> D_k0(0);

  // Specified constant chem_species specific heat
  std::vector<amrex::Real> cp_k0(0);

  // Maximum cell volume
  amrex::Real max_vol;


  void Initialize ()
  {
    amrex::ParmParse pp("chem_species");

    if (pp.contains("solve"))
    {
      pp.getarr("solve", chem_species);

      pp.get("max_solve", max_vol);

      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(chem_species.size() > 0,
          "No input provided for chem_species.solve");

      // Disable the chem_species solver if the chem_species are defined as "None" (case
      // insensitive) or 0
      if (amrex::toLower(chem_species[0]).compare("none") == 0 or
          (chem_species[0]).compare("0") == 0)
      {
        solve = 0;
      }
      else
      {
        solve = 1;
        nchem_species = chem_species.size();

        chem_species_id.resize(nchem_species);
        for (int n(0); n < nchem_species; n++) {
          chem_species_id[n] = n;
        }

        D_k0.resize(nchem_species);
      }

      if(solve)
      {
        // Get molecular weights input --------------------------------//
        for (int n(0); n < nchem_species; n++) {
          std::string name = "chem_species." + chem_species[n];
          amrex::ParmParse ppChemSpecies(name.c_str());
        }

      }
    }
  }

}
