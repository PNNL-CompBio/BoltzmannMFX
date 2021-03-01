#include <AMReX.H>
#include <AMReX_Arena.H>
#include <AMReX_Print.H>

#include <AMReX_ParmParse.H>
#include <AMReX_Utility.H>
#include <AMReX_ParallelDescriptor.H>

#include <bmx_fluid_parms.H>
#include <bmx_chem_species_parms.H>

namespace FLUID
{
  // Flag to solve fluid equations
  int solve;

  // Flag to solve chem_species fluid equations
  int solve_chem_species(0);

  // Fluid phase chem_species names
  amrex::Vector<std::string> chem_species;

  // ChemSpecies unique identifying code
  std::vector<int> chem_species_id;

  // Total number of fluid chem_species
  int nchem_species(0);

  // Specified constant gas phase chem_species diffusion coefficients
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

      // Fluid chem_species inputs
      if (ppFluid.contains("chem_species"))
      {
        solve_chem_species = 1;

        ppFluid.getarr("chem_species", chem_species);

        nchem_species = chem_species.size();

        amrex::Print() << " " << std::endl;

        amrex::Print() << " Reading in " << chem_species.size() << " chem_species from the inputs file" << std::endl;
        for (int n = 0; n < nchem_species; n++)
           amrex::Print() << "chem_species " << n << " is " << chem_species[n] << std::endl;

        D_gk0.resize(nchem_species);

        ppFluid.getarr("chem_species_diff", D_gk0);
        for (int n = 0; n < nchem_species; n++)
           amrex::Print() << "diff coeffs for chem_species " << n << " is " << D_gk0[n] << std::endl;

        amrex::Print() << " " << std::endl;


      }
    }
  }
}
