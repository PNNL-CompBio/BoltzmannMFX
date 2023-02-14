//
//     Copyright (c) 2013 Battelle Memorial Institute
//     Licensed under modified BSD License. A copy of this license can be found
//     in the LICENSE file in the top level directory of this distribution.
//
#include <AMReX_ParmParse.H>
#include <AMReX_Vector.H>

#include <bmx_diffusion_op.H>
#include <bmx_bc_parms.H>
#include <bmx_chem_species_parms.H>
#include <bmx_fluid_parms.H>

using namespace amrex;

//
// Constructor:
// We set up everything which doesn't change between timesteps here
//
DiffusionOp::DiffusionOp (AmrCore* _amrcore,
                          std::array<amrex::LinOpBCType,3> a_chem_speciesbc_lo,
                          std::array<amrex::LinOpBCType,3> a_chem_speciesbc_hi,
                          int _nghost)
{
    if(verbose > 0)
        amrex::Print() << "Constructing DiffusionOp class" << std::endl;

    nghost = _nghost;

    m_chem_speciesbc_lo = a_chem_speciesbc_lo;
    m_chem_speciesbc_hi = a_chem_speciesbc_hi;

    // Get inputs from ParmParse
    readParameters();

    // Actually do the setup work here
    setup(_amrcore);
}

void DiffusionOp::setup (AmrCore* _amrcore)
{
    // The amrcore boxArray and DistributionMap change when we regrid so we must
    // pass the new object in here.
    amrcore = _amrcore;

    grids = amrcore->boxArray();
    dmap  = amrcore->DistributionMap();

    int finest_level = amrcore->finestLevel();

    const int nchem_species = FLUID::nchem_species;

    chem_species_phi.resize(finest_level + 1);
    chem_species_rhs.resize(finest_level + 1);
    chem_species_b.resize  (finest_level + 1);

    for(int lev = 0; lev <= finest_level; lev++)
    {
        // One ghost cell needed for solution "chem_species_phi"
        chem_species_phi[lev].reset(new MultiFab(grids[lev], dmap[lev], nchem_species, 1, MFInfo()));

        // No ghost cells needed for rhs
        chem_species_rhs[lev].reset(new MultiFab(grids[lev], dmap[lev], nchem_species, 0, MFInfo()));

        // No ghost cells needed for face-based coeff array
        for(int dir = 0; dir < AMREX_SPACEDIM; dir++)
        {
              BoxArray edge_ba = grids[lev];
              edge_ba.surroundingNodes(dir);
              chem_species_b[lev][dir].reset(new MultiFab(edge_ba, dmap[lev], nchem_species, 0, MFInfo()));
        }
    }

    LPInfo info;
    info.setMaxCoarseningLevel(mg_max_coarsening_level);

    amrex::Print() << "Initializing solver with " << nchem_species << " components " << std::endl;
    chem_species_matrix.reset(new MLABecLaplacian(amrcore->Geom(0,finest_level), grids, dmap, info, {}, nchem_species));

    chem_species_matrix->setMaxOrder(2);
    chem_species_matrix->setDomainBC(m_chem_speciesbc_lo, m_chem_speciesbc_hi);
}

DiffusionOp::~DiffusionOp ()
{}

void DiffusionOp::readParameters ()
{
    ParmParse pp("diffusion");

    pp.query("verbose_solver", verbose);
    pp.query("verbose", mg_verbose);
    pp.query("bottom_verbose", mg_bottom_verbose);
    pp.query("maxiter", mg_maxiter);
    pp.query("bottom_maxiter", mg_bottom_maxiter);
    pp.query("mg_max_fmg_iter", mg_max_fmg_iter);
    pp.query("mg_max_coarsening_level", mg_max_coarsening_level);
    pp.query("rtol", mg_rtol);
    pp.query("atol", mg_atol);
    pp.query("bottom_solver", bottom_solver_type);
}


//
// Set the user-supplied settings for the MLMG solver
// (this must be done every time step, since MLMG is created after updating matrix
//
void DiffusionOp::setSolverSettings (MLMG& solver)
{
    // The default bottom solver is BiCG
    if(bottom_solver_type == "smoother")
    {
        solver.setBottomSolver(MLMG::BottomSolver::smoother);
    }
    else if(bottom_solver_type == "hypre")
    {
        solver.setBottomSolver(MLMG::BottomSolver::hypre);
    }
        // Maximum iterations for MultiGrid / ConjugateGradients
        solver.setMaxIter(mg_maxiter);
        solver.setMaxFmgIter(mg_max_fmg_iter);
        solver.setBottomMaxIter(mg_bottom_maxiter);

        // Verbosity for MultiGrid / ConjugateGradients
        solver.setVerbose(mg_verbose);
        solver.setBottomVerbose(mg_bottom_verbose);

        // This ensures that ghost cells of phi are correctly filled when
        // returned from the solver
        solver.setFinalFillBC(true);
}

void DiffusionOp::ComputeLapX (const Vector< MultiFab* >& lapX_out,
                               const Vector< MultiFab* >& X_k_in,
                               const Vector< MultiFab const*>& D_k_in,
                               const Vector< MultiFab const*>& vf_in)
{
  BL_PROFILE("DiffusionOp::ComputeLapX");

  int finest_level = amrcore->finestLevel();

  // We want to return div (D_k grad)) phi
  chem_species_matrix->setScalars(0.0, -1.0);

  define_coeffs_on_faces(D_k_in, vf_in);

  for (int lev = 0; lev <= finest_level; lev++)
  {
    chem_species_matrix->setBCoeffs(lev, GetArrOfConstPtrs(chem_species_b[lev]));

    chem_species_matrix->setLevelBC(lev, GetVecOfConstPtrs(X_k_in)[lev]);
  }

  MLMG solver(*chem_species_matrix);

  solver.apply(lapX_out, X_k_in);
}
