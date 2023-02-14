//
//     Copyright (c) 2013 Battelle Memorial Institute
//     Licensed under modified BSD License. A copy of this license can be found
//     in the LICENSE file in the top level directory of this distribution.
//
#include <AMReX.H>
#include <AMReX_Arena.H>
#include <AMReX_Print.H>

#include <AMReX_Vector.H>
#include <AMReX_ParmParse.H>
#include <bmx_dem_parms.H>

namespace DEM
{
    int solve = 0;

    int NPHASE;

    amrex::Real dtsolid;

    amrex::Real small_number = 1.0e-15;
    amrex::Real large_number = 1.0e32;
    amrex::Real eps = std::numeric_limits<amrex::Real>::epsilon();

    amrex::Real neighborhood = 1.0;

    void Initialize ()
    {

      amrex::ExecOnFinalize(Finalize);

      amrex::ParmParse ppDEM("dem");

      // Names of the solids used to build input regions.
      amrex::Vector<std::string> names;

      ppDEM.queryarr("solve", names);

      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(names.size() >= 1,
           "DEM solver not specified: Input dem.solve is undefined!");

      solve = 1;
      for(int lc=0; lc < names.size(); ++lc){
        if (amrex::toLower(names[0]).compare("none") == 0 or
          (names[0]).compare("0") == 0) solve = 0;
      }

      // You can't name a solids "None" or "0" -- you just can't
      if( solve == 0 && names.size() > 1 ){
        amrex::Abort("Invalid input: One or more DEM solids defined"
                     "but, the solver is disabled!");
      }

      if( solve ) {

        // Store the total number of solids
        NPHASE = names.size();

      }
    }

    void Finalize ()
    {
    }
}
