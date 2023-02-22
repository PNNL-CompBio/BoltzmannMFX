//
//     Copyright (c) 2013 Battelle Memorial Institute
//     Licensed under modified BSD License. A copy of this license can be found
//     in the LICENSE file in the top level directory of this distribution.
//

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
amrex::Real BMXCellInteraction::p_bond_strength = 0.0;
amrex::Real BMXCellInteraction::p_bond_cutoff = 5.0;
amrex::Real BMXCellInteraction::p_viscosity = 20.0;
std::vector<amrex::Real> BMXCellInteraction::p_force_params;

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
void BMXCellInteraction::setParams(const char* /*file*/)
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
  pp.get("viscous_drag",p_viscosity);

  p_bond_strength = 5.0;
  pp.get("bond_strength",p_bond_strength);

  p_bond_cutoff = 0.00001;
  pp.get("bond_cutoff",p_bond_cutoff);

  p_z_wall = FLUID::surface_location;

  p_force_params.clear();
  p_force_params.push_back(p_bndry_width);
  p_force_params.push_back(p_stiffness);
  p_force_params.push_back(p_z_wall);
  p_force_params.push_back(p_z_bndry_width);
  p_force_params.push_back(p_z_stiffness);
  p_force_params.push_back(p_z_gravity);
  p_force_params.push_back(p_bond_strength);
  p_force_params.push_back(p_bond_cutoff);
  p_force_params.push_back(p_viscosity);
  amrex::Print() << "SURFACE LOCATION: "<<p_z_wall<<'\n';
}

/**
 * Return a vector containing force parameters
 * @return force parameters
 */
std::vector<Real> BMXCellInteraction::getForceParams()
{
  return p_force_params;
}
