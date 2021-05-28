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
   * Calculate the maximum separation distance for which there is any
   * interaction
   * @param p1 parameters describing particle 1
   * @param p2 parameters describing particle 2
   * @return maximum interatction distance between particles 1 and 2
   */
  Real maxInteractionDistance(const Real *p1, const Real *p2);

  /**
   * Calculate the force between 2 particles
   * Calculate the force between 2 particles. r12 is computed outside this routine
   * and periodic boundary conditions, if applicable, have already been applied
   * @param r12 position of particle 1 - position of particle 2
   * @param par1 parameters describing particle 1
   * @param par2 parameters describing particle 2
   * @param frc force between particle 1 and particle 2
   */
  void evaluateForce(const Real *r12, Real *par1, Real *par2, Real *frc);

  /**
   * Calculate force between growth surface and particle
   * @param pos position of particle
   * @param par parameters describing particle
   * @param frc force on particle from wall
   */
  void evaluateSurfaceForce(const Real *pos, Real *par, Real *frc);
  protected:
  /**
   * Constructor
   */
  BMXCellInteraction();

  /**
   * Destructor
   */
  ~BMXCellInteraction();

  private:

  Real p_bndry_width;   // cell surface interaction parameter
#if 1
  Real p_stiffness;     // strength of cell interaction. This parameter assumes
                        // viscous drag is incorporated into force
#else
  Real p_r_stiffness;   // strength of repulsive cell interaction. This parameter assumes
                        // viscous drag is incorporated into force
  Real p_a_stiffness;   // strength of adhesive cell interaction. This parameter assumes
                        // viscous drag is incorporated into force
#endif

  Real p_zwall;         // z-location of surface
  Real p_z_bndry_width; // cell-growth surface interaction parameter
#if 1
  Real p_z_stiffness;   // strength of cell-growth surface interaction
#else
  Real p_z_r_stiffness; // strength of repulsive cell-growth surface interaction
  Real p_z_a_stiffness; // strength of adhesive cell-growth surface interaction
#endif

  static BMXCellInteraction *p_instance;
};
#endif
