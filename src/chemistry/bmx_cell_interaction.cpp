#include <stdio.h>
#include <math.h>
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
  p_bndry_width = 0.05;
  p_stiffness = 10.0;

  p_z_bndry_width = 0.03;
  p_z_stiffness = 10.0;

  p_zwall = FLUID::surface_location;
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
 * @param r12 position of particle 1 - position of particle 2
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
  Real F = p_stiffness*(rn-rA)*(rn-rA)*(rn-rS);
  frc[0] = F*rx;
  frc[1] = F*ry;
  frc[2] = F*rz;
  //printf("fx: %12.6f fy: %12.6f fz: %12.6f\n",frc[0],frc[1],frc[2]);
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
  if (dz > 0.0) {
    frc[2] = -p_z_stiffness*(dz-rA)*(dz-rA)*(dz-rS);
  } else {
    frc[2] = 0.0;
  }
}
