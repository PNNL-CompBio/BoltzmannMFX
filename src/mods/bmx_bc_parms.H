//
//     Copyright (c) 2013 Battelle Memorial Institute
//     Licensed under modified BSD License. A copy of this license can be found
//     in the LICENSE file in the top level directory of this distribution.
//
#ifndef BMX_BC_PARMS_H_
#define BMX_BC_PARMS_H_

#include <AMReX_REAL.H>
#include <AMReX_RealBox.H>
#include <AMReX_Geometry.H>

#include <AMReX_Vector.H>
#include <AMReX_LO_BCTYPES.H>

#include <bmx_fluid_parms.H>

namespace BC
{
  struct BC_t {
    // BC type (mi, po, nsw)
    int type;
  };

  // Lists of BCs applied to the domain extent
  extern amrex::Vector<int> bc_xlo, bc_xhi;
  extern amrex::Vector<int> bc_ylo, bc_yhi;
  extern amrex::Vector<int> bc_zlo, bc_zhi;

  // Flags for DEM domain extent reinforcement.
  extern int domain_bc[6];

  // Use AMReX LinOpBCType to define the domain extents
  extern std::array<amrex::LinOpBCType,3> diff_chem_species_lobc;
  extern std::array<amrex::LinOpBCType,3> diff_chem_species_hibc;

  // Data structure for BC data
  extern amrex::Vector<BC_t> bc;

  void Initialize (amrex::Geometry& geom);
}

#endif
