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

  // Maximum segment length
  amrex::Real max_len;

  // Post-split length
  amrex::Real split_len;

  // Maximum segment radius
  amrex::Real max_rad;

  // Maximum branching probability
  amrex::Real brnch_prob;

  void Initialize ()
  {
    amrex::ParmParse pp("chem_species");

    pp.get("max_vol", max_vol);
    max_len = 0.0035;
    pp.query("max_seg_length", max_len);
    split_len = 0.0030;
    pp.query("seg_split_length", split_len);
    max_rad = 0.00025;
    pp.query("max_seg_radius", max_rad);
    brnch_prob = 0.0001;
    pp.query("branching_probability", brnch_prob);

    if (pp.contains("solve"))
    {
      pp.getarr("solve", chem_species);

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
