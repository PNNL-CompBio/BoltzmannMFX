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
amrex::Real BMXCellInteraction::p_fungi_bndry_width = 0.0;
amrex::Real BMXCellInteraction::p_fungi_stiffness = 0.0;
amrex::Real BMXCellInteraction::p_cell_bndry_width = 0.0;
amrex::Real BMXCellInteraction::p_cell_stiffness = 0.0;
amrex::Real BMXCellInteraction::p_fungi_z_bndry_width = 0.0;
amrex::Real BMXCellInteraction::p_fungi_z_stiffness = 0.0;
amrex::Real BMXCellInteraction::p_cell_z_bndry_width = 0.0;
amrex::Real BMXCellInteraction::p_cell_z_stiffness = 0.0;
amrex::Real BMXCellInteraction::p_z_wall = 0.0;
amrex::Real BMXCellInteraction::p_bond_strength = 0.0;
amrex::Real BMXCellInteraction::p_bond_cutoff = 5.0;
amrex::Real BMXCellInteraction::p_viscosity = 20.0;
amrex::Real BMXCellInteraction::p_ran_scale = 0.0;
amrex::Real BMXCellInteraction::p_fluc_mix = 0.0;
amrex::Gpu::DeviceVector<amrex::Real> BMXCellInteraction::p_force_params;

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

  pp.get("fungi_boundary_width",p_fungi_bndry_width);
  pp.get("fungi_stiffness",p_fungi_stiffness);
  pp.get("cell_boundary_width",p_cell_bndry_width);
  pp.get("cell_stiffness",p_cell_stiffness);
  pp.get("fungi_wall_boundary_width",p_fungi_z_bndry_width);
  pp.get("fungi_wall_stiffness",p_fungi_z_stiffness);
  pp.get("cell_wall_boundary_width",p_cell_z_bndry_width);
  pp.get("cell_wall_stiffness",p_cell_z_stiffness);

  p_z_gravity = 0.0;
  pp.query("gravity",p_z_gravity);
  pp.get("viscous_drag",p_viscosity);

  p_bond_strength = 5.0;
  pp.get("bond_strength",p_bond_strength);

  p_bond_cutoff = 0.00001;
  pp.get("bond_cutoff",p_bond_cutoff);
  p_ran_scale = 0.0;
  pp.query("fluctuation_scale",p_ran_scale);
  p_fluc_mix = 0.0;
  pp.query("fluctuation_mixing",p_fluc_mix);

  p_z_wall = FLUID::surface_location;

  p_force_params.clear();
  p_force_params.push_back(p_fungi_bndry_width);
  p_force_params.push_back(p_fungi_stiffness);
  p_force_params.push_back(p_cell_bndry_width);
  p_force_params.push_back(p_cell_stiffness);
  p_force_params.push_back(p_fungi_z_bndry_width);     // 5
  p_force_params.push_back(p_fungi_z_stiffness);
  p_force_params.push_back(p_cell_z_bndry_width);
  p_force_params.push_back(p_cell_z_stiffness);
  p_force_params.push_back(p_z_wall);
  p_force_params.push_back(p_z_gravity);               // 10
  p_force_params.push_back(p_bond_strength);
  p_force_params.push_back(p_bond_cutoff);
  p_force_params.push_back(p_viscosity);
  p_force_params.push_back(p_ran_scale);
  p_force_params.push_back(p_fluc_mix);                // 15
  amrex::Print() << "SURFACE LOCATION: "<<p_z_wall<<'\n';
}

/**
 * Return a vector containing force parameters
 * @return force parameters
 */
amrex::Gpu::DeviceVector<Real> BMXCellInteraction::getForceParams()
{
  return p_force_params;
}
