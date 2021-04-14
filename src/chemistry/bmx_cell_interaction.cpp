#include <stdio.h>
#include <math.h>
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
  p_stiffness = 100.0;
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
  Real F = p_stiffness*(rn-rS)*(rn-rS)*(rn-rS);
  frc[0] = F*rx;
  frc[1] = F*ry;
  frc[2] = F*rz;
}
