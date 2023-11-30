//
//     Copyright (c) 2013 Battelle Memorial Institute
//     Licensed under modified BSD License. A copy of this license can be found
//     in the LICENSE file in the top level directory of this distribution.
//
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

  // Branching probability
  amrex::Real brnch_prob;

  // Tip-splitting probability
  amrex::Real split_prob;

  // Frequency to compute radius of gyration and center of mass
  int rg_frequency;

  // Number of bins in density profile
  int dens_prof_bins;

  // Maximum radius of density profile
  amrex::Real dens_prof_max;

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
    split_prob = 0.1;
    pp.query("splitting_probability", split_prob);
    rg_frequency = 100;
    pp.query("rg_frequency", rg_frequency);
    dens_prof_bins = 50;
    pp.query("density_profile_bins", dens_prof_bins);
    dens_prof_max = 0.0150;
    pp.query("density_profile_maximum", dens_prof_max);

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
