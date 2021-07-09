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

  static Real p_bndry_width;   // cell surface interaction parameter
#if 1
  static Real p_stiffness;     // strength of cell interaction. This parameter assumes
                               // viscous drag is incorporated into force
#else
  static Real p_r_stiffness;   // strength of repulsive cell interaction. This parameter assumes
                            // viscous drag is incorporated into force
  static Real p_a_stiffness;   // strength of adhesive cell interaction. This parameter assumes
                               // viscous drag is incorporated into force
#endif

  static Real p_z_wall;        // z-location of surface
  static Real p_z_bndry_width; // cell-growth surface interaction parameter

#if 1
  static Real p_z_stiffness;   // strength of cell-growth surface interaction
#else
  static Real p_z_r_stiffness; // strength of repulsive cell-growth surface interaction
  static Real p_z_a_stiffness; // strength of adhesive cell-growth surface interaction
#endif

  static Real p_z_gravity;     // factor for gravity

  private:

  static BMXCellInteraction *p_instance;
};
#endif
