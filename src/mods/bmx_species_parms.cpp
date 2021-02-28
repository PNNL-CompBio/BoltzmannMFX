#include <AMReX.H>
#include <AMReX_Arena.H>
#include <AMReX_Print.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Utility.H>
#include <AMReX_ParallelDescriptor.H>

#include <bmx_species_parms.H>

namespace SPECIES
{
  // Flag to solve species equations
  int solve;

  // Total number of species
  int nspecies(0);

  // Species names
  std::vector<std::string> species(0);

  // Species unique identifying code (at the moment = their index in the input
  // entries)
  std::vector<int> species_id(0);

  // Specified species diffusion coefficients
  std::vector<amrex::Real> D_k0(0);

  // Specified constant species specific heat
  std::vector<amrex::Real> cp_k0(0);


  void Initialize ()
  {
    amrex::ParmParse pp("species");

    if (pp.contains("solve"))
    {
      pp.getarr("solve", species);

      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(species.size() > 0,
          "No input provided for species.solve");

      // Disable the species solver if the species are defined as "None" (case
      // insensitive) or 0
      if (amrex::toLower(species[0]).compare("none") == 0 or
          (species[0]).compare("0") == 0)
      {
        solve = 0;
      }
      else
      {
        solve = 1;
        nspecies = species.size();

        species_id.resize(nspecies);
        for (int n(0); n < nspecies; n++) {
          species_id[n] = n;
        }

        D_k0.resize(nspecies);
      }

      if(solve)
      {
        // Get molecular weights input --------------------------------//
        for (int n(0); n < nspecies; n++) {
          std::string name = "species." + species[n];
          amrex::ParmParse ppSpecies(name.c_str());
        }

      }
    }
  }

}
