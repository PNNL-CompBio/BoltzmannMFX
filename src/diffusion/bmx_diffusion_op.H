#ifndef DIFFUSION_EQUATION_H_
#define DIFFUSION_EQUATION_H_

#include <AMReX_AmrCore.H>
#include <AMReX_MLMG.H>
#include <AMReX_MLABecLaplacian.H>

//
// Solver for the implicit part of the diffusion equation:
//
// ( 1 - dt / rho * div ( eta grad ) ) u* = rhs
//
// Input "vel" holds the right hand side of this equation:
//
//  rhs = u + dt * ( conv + div ( eta (grad u)^T ) / rho + g - grad(p + p0) / rho )
//
// This is replaced by the solution of the matrix solve.
// Output "vel" therefore holds u*
//
// Note: we actually solve the above equation multiplied by the density ro.
//

class DiffusionOp
{
  public:
    // Constructor, with initialisation
    DiffusionOp (amrex::AmrCore* _amrcore,
                 std::array<amrex::LinOpBCType,AMREX_SPACEDIM> a_chem_specieslbc_lo,
                 std::array<amrex::LinOpBCType,AMREX_SPACEDIM> a_chem_specieslbc_hi,
                 int _nghost);

    // Destructor
    ~DiffusionOp ();

    // Actually do the setup work here
    void setup (amrex::AmrCore* _amrcore);

    // Read input from ParmParse (solver settings)
    void readParameters ();

    // Set user-supplied solver settings (must be done every time step)
    void setSolverSettings (amrex::MLMG& solver);

    // Species mass fraction diffusion
    void diffuse_chem_species (amrex::Vector< amrex::MultiFab* > X_gk,
                          const amrex::Vector< amrex::MultiFab* > D_gk,
                          amrex::Real theta, amrex::Real dt);

    // Explicit tensor temperature term
    void ComputeLapX (const amrex::Vector< amrex::MultiFab* >& lapX_out,
                      const amrex::Vector< amrex::MultiFab* >& X_gk_in,
                      const amrex::Vector< amrex::MultiFab const*>& D_gk);

  private:
    // AmrCore data
    amrex::AmrCore* amrcore;

    amrex::Vector<amrex::Geometry> geom;
    amrex::Vector<amrex::BoxArray> grids;
    amrex::Vector<amrex::DistributionMapping> dmap;

    int nghost = 0;

    // Boundary conditions
    std::array<amrex::LinOpBCType,3> m_chem_speciesbc_lo;
    std::array<amrex::LinOpBCType,3> m_chem_speciesbc_hi;

    // Internal data used in the matrix solve
    //
    // MultiLevel EmbeddedBoundary cell-centered Laplacian:
    //
    // ( alpha a - beta div ( b grad ) ) phi = rhs
    //
    amrex::Vector<amrex::Array<std::unique_ptr<amrex::MultiFab>, AMREX_SPACEDIM>> b;

    amrex::Vector<std::unique_ptr<amrex::MultiFab>> phi;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> rhs;

    amrex::Vector<std::unique_ptr<amrex::MultiFab>> chem_species_phi;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> chem_species_rhs;
    //
    // MultiLevel cell-centered Laplacian:
    //
    // ( alpha a - beta div ( b grad ) ) phi = rhs
    //
    std::unique_ptr<amrex::MLABecLaplacian> chem_species_matrix;

    amrex::Vector<amrex::Array<std::unique_ptr<amrex::MultiFab>, AMREX_SPACEDIM>> chem_species_b;

    // DiffusionOp verbosity
    int verbose = 0;

    // Options to control MLMG behavior
    int mg_verbose = 0;
    int mg_bottom_verbose = 0;
    int mg_maxiter = 100;
    int mg_bottom_maxiter = 100;
    int mg_max_fmg_iter = 0;
    int mg_max_coarsening_level = 100;
    amrex::Real mg_rtol = 1.0e-11;
    amrex::Real mg_atol = 1.0e-14;
    std::string bottom_solver_type = "bicgstab";
};

#endif
