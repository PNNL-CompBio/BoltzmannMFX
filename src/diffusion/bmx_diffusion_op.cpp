#include <AMReX_MultiFabUtil.H>
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

    geom  = amrcore->Geom();
    grids = amrcore->boxArray();
    dmap  = amrcore->DistributionMap();

    int max_level = amrcore->maxLevel();

    const int nchem_species_g = FLUID::nchem_species;

    // Resize and reset data
    b.resize(max_level + 1);

    phi.resize(max_level + 1);
    rhs.resize(max_level + 1);

    chem_species_phi.resize(max_level + 1);
    chem_species_rhs.resize(max_level + 1);

    for(int lev = 0; lev <= max_level; lev++)
    {
        for(int dir = 0; dir < AMREX_SPACEDIM; dir++)
        {
            BoxArray edge_ba = grids[lev];
            edge_ba.surroundingNodes(dir);
            b[lev][dir].reset(new MultiFab(edge_ba, dmap[lev], 1, nghost, MFInfo()));
        }

        phi[lev].reset(new MultiFab(grids[lev], dmap[lev], 3, 1, MFInfo()));

        // No ghost cells needed for rhs
        rhs[lev].reset(new MultiFab(grids[lev], dmap[lev], 3, 0, MFInfo()));

        chem_species_phi[lev].reset(new MultiFab(grids[lev], dmap[lev], nchem_species_g, 1, MFInfo()));

        // No ghost cells needed for rhs
        chem_species_rhs[lev].reset(new MultiFab(grids[lev], dmap[lev], nchem_species_g, 0, MFInfo()));
    }

    LPInfo info;
    info.setMaxCoarseningLevel(mg_max_coarsening_level);

    amrex::Print() << "Initializing solver with " << nchem_species_g << " components " << std::endl;
    chem_species_matrix.reset(new MLABecLaplacian(geom, grids, dmap, info, {}, nchem_species_g));

    chem_species_matrix->setMaxOrder(2);

    chem_species_matrix->setDomainBC(m_chem_speciesbc_lo, m_chem_speciesbc_hi);

    chem_species_b.resize(max_level + 1);

    for(int lev = 0; lev <= max_level; lev++)
    {
        for(int dir = 0; dir < AMREX_SPACEDIM; dir++)
        {
              BoxArray edge_ba = grids[lev];
              edge_ba.surroundingNodes(dir);
              chem_species_b[lev][dir].reset(new MultiFab(edge_ba, dmap[lev], nchem_species_g, nghost, MFInfo()));
        }
    }
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
                               const Vector< MultiFab* >& X_gk_in,
                               const Vector< MultiFab const*>& D_gk_in)
{
  BL_PROFILE("DiffusionOp::ComputeLapX");

  int finest_level = amrcore->finestLevel();

  // Number of fluid chem_species
  const int nchem_species_g = FLUID::nchem_species;

  // We want to return div (D_gk grad)) phi
  chem_species_matrix->setScalars(0.0, -1.0);

  Vector<BCRec> bcs_X; 
  bcs_X.resize(3*nchem_species_g);

  // Zero out the coefficients in the high-z part
  for (int lev = 0; lev <= finest_level; lev++)
  {
    MultiFab b_coeffs(X_gk_in[lev]->boxArray(), X_gk_in[lev]->DistributionMap(), nchem_species_g, 1, MFInfo());
    b_coeffs.setVal(0.);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*X_gk_in[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
      Box const& bx = mfi.growntilebox(IntVect(1,1,1));

      if (bx.ok())
      {
        Array4<Real const> const& D_gk_arr     = D_gk_in[lev]->const_array(mfi);
        Array4<Real      > const& b_coeffs_arr = b_coeffs.array(mfi);

        amrex::ParallelFor(bx, nchem_species_g, [=]
          AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            b_coeffs_arr(i,j,k,n) = D_gk_arr(i,j,k,n);
        });
      }
    }

    average_cellcenter_to_face( GetArrOfPtrs(chem_species_b[lev]), b_coeffs, geom[lev], nchem_species_g );

    // Zero out the coefficients in the high-z part
    const Box& domain = geom[lev].Domain();
    const int zhi = domain.bigEnd()[2];

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*X_gk_in[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
      Box const& bx = mfi.tilebox();

      if (bx.ok())
      {
        Array4<Real> const& bx_arr = chem_species_b[lev][0]->array(mfi);
        Array4<Real> const& by_arr = chem_species_b[lev][1]->array(mfi);
        Array4<Real> const& bz_arr = chem_species_b[lev][2]->array(mfi);

        Box const& xbx = mfi.growntilebox(IntVect(1,0,0));
        amrex::ParallelFor(xbx, nchem_species_g, [=]
          AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            if (k > zhi/2)
                bx_arr(i,j,k,n) = 0.;
        });

        Box const& ybx = mfi.growntilebox(IntVect(0,1,0));
        amrex::ParallelFor(ybx, nchem_species_g, [=]
          AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            if (k > zhi/2)
                by_arr(i,j,k,n) = 0.;
        });

        Box const& zbx = mfi.growntilebox(IntVect(0,0,1));
        amrex::ParallelFor(zbx, nchem_species_g, [=]
          AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            if (k > zhi/2)
                bz_arr(i,j,k,n) = 0.;
        });
      }

      Box grown_bx(bx); grown_bx.growLo(2,1); grown_bx.setBig(2,0);

      // Set Dirichlet conditions on solution
      if (grown_bx.ok())
      {
        Array4<Real> const& soln_arr = X_gk_in[lev]->array(mfi);
        amrex::ParallelFor(grown_bx, nchem_species_g, [=]
          AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            soln_arr(i,j,k,n) = 1.;
        });
      }
    }

    chem_species_matrix->setBCoeffs(lev, GetArrOfConstPtrs(chem_species_b[lev]));

    chem_species_matrix->setLevelBC(lev, GetVecOfConstPtrs(X_gk_in)[lev]);
  }

  MLMG solver(*chem_species_matrix);

  solver.apply(lapX_out, X_gk_in);
}
