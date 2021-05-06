#include <bmx.H>

#include <AMReX_BC_TYPES.H>
#include <AMReX_Box.H>

#include <bmx_fluid_parms.H>
#include <bmx_chem_species_parms.H>

std::string      bmx::particle_init_type   = "AsciiFile";
std::string      bmx::load_balance_type    = "KnapSack";
std::string      bmx::knapsack_weight_type = "RunTimeCosts";
int              bmx::load_balance_fluid   = 1;
int              bmx::knapsack_nmax        = 128;
DepositionScheme bmx::m_deposition_scheme;
amrex::Real      bmx::m_deposition_scale_factor = 1.0;

int bmx::nlev  = 1;

int  bmx::plot_int        = -1;
Real bmx::plot_per_approx = -1.;
Real bmx::plot_per_exact  = -1.;

// Destructor
bmx::~bmx ()
{
  if (DEM::solve)
    delete pc;

  for (int lev(0); lev < nlev; ++lev)
  {
    // Face-based coefficients b in MAC projection and implicit diffusion solve
    delete bcoeff[lev][0];
    delete bcoeff[lev][1];
    delete bcoeff[lev][2];

    // Boundary conditions types
    delete bc_ilo[lev];
    delete bc_ihi[lev];
    delete bc_jlo[lev];
    delete bc_jhi[lev];
    delete bc_klo[lev];
    delete bc_khi[lev];
  }

  // used if load_balance_type == "KnapSack"
  for (int lev = 0; lev < particle_cost.size(); ++lev)
    delete particle_cost[lev];

  for (int lev = 0; lev < fluid_cost.size(); ++lev)
    delete fluid_cost[lev];
} 

// Constructor
bmx::bmx ()
{
    // NOTE: Geometry on all levels has just been defined in the AmrCore
    // constructor. No valid BoxArray and DistributionMapping have been defined.
    // But the arrays for them have been resized.

    /****************************************************************************
     *                                                                          *
     * Set max number of levels (nlev)                                         *
     *                                                                          *
     ***************************************************************************/

    nlev = maxLevel() + 1;
    amrex::Print() << "Number of levels: " << nlev << std::endl;

    /****************************************************************************
     *                                                                          *
     * Initialize time steps                                                    *
     *                                                                          *
     ***************************************************************************/

    t_old.resize(nlev,-1.e100);
    t_new.resize(nlev,0.0);


    /****************************************************************************
     *                                                                          *
     * Initialize boundary conditions (used by fill-patch)                      *
     *                                                                          *
     ***************************************************************************/

    bcs_X.resize(0); // X_k
    bcs_D.resize(0); // D_k

    //___________________________________________________________________________
    // Boundary conditions used for level-sets

    // walls (Neumann)
    int bc_lo[] = {FOEXTRAP, FOEXTRAP, FOEXTRAP};
    int bc_hi[] = {FOEXTRAP, FOEXTRAP, FOEXTRAP};

    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        // lo-side BCs
        if (bc_lo[idim] == BCType::int_dir  ||  // periodic uses "internal Dirichlet"
            bc_lo[idim] == BCType::foextrap ||  // first-order extrapolation
            bc_lo[idim] == BCType::ext_dir ) {  // external Dirichlet
        }
        else {
            amrex::Abort("Invalid level-set bc_lo");
        }

        // hi-side BCSs
        if (bc_hi[idim] == BCType::int_dir  ||  // periodic uses "internal Dirichlet"
            bc_hi[idim] == BCType::foextrap ||  // first-order extrapolation
            bc_hi[idim] == BCType::ext_dir ) {  // external Dirichlet
        }
        else {
            amrex::Abort("Invalid level-set bc_hi");
        }
    }

    m_X_k_bc_types["Dirichlet"] = {bc_list.get_minf()};

    Gpu::synchronize();
}

void
bmx::avgDown (int crse_lev, const MultiFab& S_fine, MultiFab& S_crse)
{
    BL_PROFILE("bmx::avgDown()");

    average_down(S_fine, S_crse, 0, S_fine.nComp(), refRatio(crse_lev));
}

Vector< MultiFab* > bmx::get_X_k () noexcept
{
  Vector<MultiFab*> r;
  r.reserve(m_leveldata.size());
  for (int lev = 0; lev < m_leveldata.size(); ++lev) {
    r.push_back(m_leveldata[lev]->X_k);
  }
  return r;
}

Vector< MultiFab* > bmx::get_X_k_old () noexcept
{
  Vector<MultiFab*> r;
  r.reserve(m_leveldata.size());
  for (int lev = 0; lev < m_leveldata.size(); ++lev) {
    r.push_back(m_leveldata[lev]->X_ko);
  }
  return r;
}

Vector< MultiFab* > bmx::get_D_k () noexcept
{
  Vector<MultiFab*> r;
  r.reserve(m_leveldata.size());
  for (int lev = 0; lev < m_leveldata.size(); ++lev) {
    r.push_back(m_leveldata[lev]->D_k);
  }
  return r;
}

Vector< MultiFab* > bmx::get_vf () noexcept
{
  Vector<MultiFab*> r;
  r.reserve(m_leveldata.size());
  for (int lev = 0; lev < m_leveldata.size(); ++lev) {
    r.push_back(m_leveldata[lev]->vf);
  }
  return r;
}

Vector< MultiFab const*> bmx::get_X_k_const () const noexcept
{
  Vector<MultiFab const*> r;
  r.reserve(m_leveldata.size());
  for (int lev = 0; lev < m_leveldata.size(); ++lev) {
    r.push_back(m_leveldata[lev]->X_k);
  }
  return r;
}

Vector< MultiFab const*> bmx::get_X_k_old_const () const noexcept
{
  Vector<MultiFab const*> r;
  r.reserve(m_leveldata.size());
  for (int lev = 0; lev < m_leveldata.size(); ++lev) {
    r.push_back(m_leveldata[lev]->X_ko);
  }
  return r;
}

Vector< MultiFab const*> bmx::get_D_k_const () const noexcept
{
  Vector<MultiFab const*> r;
  r.reserve(m_leveldata.size());
  for (int lev = 0; lev < m_leveldata.size(); ++lev) {
    r.push_back(m_leveldata[lev]->D_k);
  }
  return r;
}

Vector< MultiFab const*> bmx::get_vf_const () const noexcept
{
  Vector<MultiFab const*> r;
  r.reserve(m_leveldata.size());
  for (int lev = 0; lev < m_leveldata.size(); ++lev) {
    r.push_back(m_leveldata[lev]->vf);
  }
  return r;
}
