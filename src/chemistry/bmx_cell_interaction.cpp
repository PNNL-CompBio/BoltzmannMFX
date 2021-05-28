#include <stdio.h>
#include <math.h>
#include <bmx.H>
#include <bmx_fluid_parms.H>
#include <bmx_cell_interaction.H>

BMXCellInteraction *BMXCellInteraction::p_instance = NULL;

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

  p_zwall = FLUID::surface_location;
  amrex::Print() << "SURFACE LOCATION: "<<p_zwall<<'\n';
}

/**
 * Calculate the maximum separation distance for which there is any
 * interaction
 * @param p1 parameters describing particle 1
 * @param p2 parameters describing particle 2
 * @return maximum interatction distance between particles 1 and 2
 */
Real BMXCellInteraction::maxInteractionDistance(const Real *p1, const Real *p2)
{
  Real r1 = p1[realIdx::a_size];
  Real r2 = p2[realIdx::a_size];
  return r1+r2+p_bndry_width;
}
/**
 * Calculate the force between 2 particles. r12 is computed outside this routine
 * and periodic boundary conditions, if applicable, have already been applied
 * @param r12 position of particle 2 - position of particle 1
 * @param par1 parameters describing particle 1
 * @param par2 parameters describing particle 2
 * @param frc force between particle 1 and particle 2
 */
void BMXCellInteraction::evaluateForce(const Real *r12, Real *par1,
    Real *par2, Real *frc)
{
  Real rS = par1[realIdx::a_size] + par2[realIdx::a_size];
  Real rA = rS + p_bndry_width;
  Real rx, ry, rz, rn;
  rx = r12[0];
  ry = r12[1];
  rz = r12[2];
  rn = sqrt(rx*rx+ry*ry+rz*rz);
  rx /= rn;
  ry /= rn;
  rz /= rn;
#if 1
  Real F = p_stiffness*(rn-rA)*(rn-rA)*(rn-rS);
#else
  Real F;
  if (rn < rS) {
    F = p_r_stiffness*(rn-rA)*(rn-rA)*(rn-rS);
  } else {
    F = p_a_stiffness*(rn-rS)*(rA-rn);
  }
#endif
  frc[0] = F*rx;
  frc[1] = F*ry;
  frc[2] = F*rz;
}

/**
 * Calculate force between growth surface and particle
 * @param pos position of particle
 * @param par parameters describing particle
 * @param frc force on particle from wall
 */
void BMXCellInteraction::evaluateSurfaceForce(const Real *pos, Real *par, Real *frc)
{
  frc[0] = 0.0;
  frc[1] = 0.0;
  Real z = pos[2];
  Real dz = z - p_zwall; 
  Real rS = par[realIdx::a_size];
  Real rA = rS + p_z_bndry_width;
  if (dz > 0.0 && dz < rA) {
#if 1
    frc[2] = -p_z_stiffness*(dz-rA)*(dz-rA)*(dz-rS);
#else
    if (dz < rS) {
      frc[2] = -p_z_r_stiffness*(dz-rA)*(dz-rA)*(dz-rS);
    } else {
      frc[2] = -p_z_a_stiffness*(rA-dz)*(dz-rS);
    }
#endif
  } else if (dz <= 0) {
#if 1
    frc[2] = p_z_stiffness*rA*rA*rS;
#else
    frc[2] = p_z_r_stiffness*rA*rA*rS;
#endif
  } else {
    frc[2] = 0.0;
  }
}
