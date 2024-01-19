//
//     Copyright (c) 2013 Battelle Memorial Institute
//     Licensed under modified BSD License. A copy of this license can be found
//     in the LICENSE file in the top level directory of this distribution.
//
#ifndef _BMX_CELL_INTERACTION_H_
#define _BMX_CELL_INTERACTION_H_

#include <AMReX.H>
#include <AMReX_REAL.H>
#include <bmx_chem.H>

using amrex::Real;

class BMXCellInteraction
{
  public:

  /**
   * Retrieve instance of BMXCellInteraction object
   */
  static BMXCellInteraction* instance();

  /**
   * Setup chemistry model by reading a parameter file
   * @param file name of parameter file used by chemistry model
   */
  void setParams(const char *file);

  /**
   * Return a vector containing force parameters
   * @return force parameters
   */
  amrex::Gpu::DeviceVector<Real> getForceParams();

  protected:
  /**
   * Constructor
   */
  BMXCellInteraction();

  /**
   * Destructor
   */
  ~BMXCellInteraction();
  
  public:

  static Real p_fungi_bndry_width;   // fungi surface interaction parameter
  static Real p_fungi_stiffness;     // strength of fungi interaction.
  static Real p_cell_bndry_width;   // cell surface interaction parameter
  static Real p_cell_stiffness;     // strength of cell interaction.

  static Real p_z_wall;        // z-location of surface

  static Real p_fungi_z_bndry_width; // fungi-growth surface interaction parameter
  static Real p_fungi_z_stiffness;   // strength of fungi-growth surface interaction
  static Real p_cell_z_bndry_width; // cell-growth surface interaction parameter
  static Real p_cell_z_stiffness;   // strength of cell-growth surface interaction

  static Real p_z_gravity;     // factor for gravity

  static Real p_bond_strength; // strength of harmonic bond linking fungi segments

  static Real p_bond_cutoff;   // cutoff radius for detecting a bond

  static Real p_viscosity;     // drag parameter for segments

  static Real p_ran_scale;     // magnitude of random velocity fluctuations

  static Real p_fluc_mix;      // fluctuation mixing coefficient

  static amrex::Gpu::DeviceVector<Real> p_force_params;  // vector containing all force params

  private:

  static BMXCellInteraction *p_instance;
};
#endif
