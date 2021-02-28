#include <AMReX.H>
#include <AMReX_Arena.H>
#include <AMReX_Print.H>

#include <AMReX_ParmParse.H>
#include <AMReX_Utility.H>
#include <AMReX_ParallelDescriptor.H>

#include <bmx_fluid_parms.H>
#include <bmx_species_parms.H>

namespace FLUID
{
  // Flag to solve fluid equations
  int solve;

  // Flag to solve species fluid equations
  int solve_species(0);

  // Fluid phase species names
  amrex::Vector<std::string> species;

  // Species unique identifying code
  std::vector<int> species_id;

  // Total number of fluid species
  int nspecies(0);

  // Specified constant gas phase species diffusion coefficients
  std::vector<amrex::Real> D_gk0(0);

  // Name to later reference when building inputs for IC/BC regions
  std::string name;


  void Initialize ()
  {
    amrex::ParmParse pp("fluid");

    std::vector<std::string> fluid_name;
    pp.queryarr("solve", fluid_name);

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(fluid_name.size() == 1,
       "Fluid solver not specified. fluid.sove = ? ");

    // Disable the fluid solver if the fluid is defined as "None"
    if (amrex::toLower(fluid_name[0]).compare("none") == 0 or fluid_name[0] == "0" ) {
      solve = 0;
    } else {
      solve = 1;
      name = fluid_name[0];
    }

    if (solve)
    {
      amrex::ParmParse ppFluid(name.c_str());

      // Fluid species inputs
      if (ppFluid.contains("species"))
      {
        solve_species = 1;

        ppFluid.getarr("species", species);

        nspecies = species.size();

        amrex::Print() << " " << std::endl;

        amrex::Print() << " Reading in " << species.size() << " species from the inputs file" << std::endl;
        for (int n = 0; n < nspecies; n++)
           amrex::Print() << "species " << n << " is " << species[n] << std::endl;

        D_gk0.resize(nspecies);

        ppFluid.getarr("species_diff", D_gk0);
        for (int n = 0; n < nspecies; n++)
           amrex::Print() << "diff coeffs for species " << n << " is " << D_gk0[n] << std::endl;

        amrex::Print() << " " << std::endl;


      }
    }
  }
}
