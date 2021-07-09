#include <stdio.h>
#include <math.h>
#include <bmx.H>
#include <bmx_fluid_parms.H>
#include <bmx_cell_interaction.H>

BMXCellInteraction *BMXCellInteraction::p_instance = NULL;
amrex::Real BMXCellInteraction::p_z_gravity = 0.0;
amrex::Real BMXCellInteraction::p_bndry_width = 0.0;
amrex::Real BMXCellInteraction::p_stiffness = 0.0;
amrex::Real BMXCellInteraction::p_z_bndry_width = 0.0;
amrex::Real BMXCellInteraction::p_z_stiffness = 0.0;
amrex::Real BMXCellInteraction::p_z_wall = 0.0;

/**
 * Retrieve instance of BMXCellInteraction object
 */
BMXCellInteraction* BMXCellInteraction::instance()
{
  if (p_instance == NULL) {
    p_instance = new BMXCellInteraction();
  }
  return p_instance;
}

/**
 * Constructor
 */
BMXCellInteraction::BMXCellInteraction()
{
}

/**
 * Destructor
 */
BMXCellInteraction::~BMXCellInteraction()
{
  delete p_instance;
}

/**
 * Setup chemistry model by reading a parameter file
 * @param file name of parameter file used by chemistry model
 */
void BMXCellInteraction::setParams(const char *file)
{
  ParmParse pp("cell_force");

  pp.get("boundary_width",p_bndry_width);
  p_z_gravity = 0.0;
  pp.query("gravity",p_z_gravity);
#if 1
  pp.get("stiffness",p_stiffness);
#else
  pp.get("repulsion_stiffness",p_r_stiffness);
  pp.get("adhesion_stiffness",p_a_stiffness);
#endif
  pp.get("wall_boundary_width",p_z_bndry_width);
#if 1
  pp.get("wall_stiffness",p_z_stiffness);
#else
  pp.get("wall_repulsion_stiffness",p_z_r_stiffness);
  pp.get("wall_adhesion_stiffness",p_z_a_stiffness);
#endif

  p_z_wall = FLUID::surface_location;
  amrex::Print() << "SURFACE LOCATION: "<<p_z_wall<<'\n';
}
