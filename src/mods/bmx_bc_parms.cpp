#include <AMReX.H>
#include <AMReX_Arena.H>
#include <AMReX_Print.H>

#include <AMReX_IntVect.H>
#include <AMReX_RealVect.H>
#include <AMReX_Geometry.H>

#include <AMReX_ParmParse.H>

#include <bmx_bc_parms.H>
#include <bmx_dem_parms.H>
#include <bmx_chem_species_parms.H>
#include <bmx_bc_list.H>

namespace BC
{

  // Direction of pressure drop (0:x, 1:y, 2:z)
  int delp_dir = -1;
  amrex::Real delp[3];

  int domain_bc[6];

  // Lists of BCs applied to the domain extent
  amrex::Vector<int> bc_xlo, bc_xhi;
  amrex::Vector<int> bc_ylo, bc_yhi;
  amrex::Vector<int> bc_zlo, bc_zhi;

  std::array<amrex::LinOpBCType,3> diff_chem_species_lobc;
  std::array<amrex::LinOpBCType,3> diff_chem_species_hibc;


  // Data structure storing individual BC information
  amrex::Vector<BC_t> bc;


  void Initialize (amrex::Geometry& geom)
  {

    BCList bc_mask;

    // Integer ids for BCs
    const int pout_ = bc_mask.get_pout();
    const int pinf_ = bc_mask.get_pinf();
    const int minf_ = bc_mask.get_minf();
    const int nsw_  = bc_mask.get_nsw();
    const int eb_   = bc_mask.get_eb();

    // Set flag to keep particles from leaving unless periodic.
    for (int dir(0); dir<3; ++dir) {
      if (geom.isPeriodic(dir)) {
        domain_bc[2*dir  ] = 0;
        domain_bc[2*dir+1] = 0;
      } else {
        domain_bc[2*dir  ] = 1;
        domain_bc[2*dir+1] = 1;
      }
    }


    // Default all sides of the domain to Neumann
    for (int dir=0; dir < 3; dir ++ ){
      if (geom.isPeriodic(dir))
      {

        diff_chem_species_lobc[dir] = amrex::LinOpBCType::Periodic;
        diff_chem_species_hibc[dir] = amrex::LinOpBCType::Periodic;

      } else {

        diff_chem_species_lobc[dir] = amrex::LinOpBCType::Dirichlet;
        diff_chem_species_hibc[dir] = amrex::LinOpBCType::Dirichlet;

      }
    }
  }// END Initialize

}
