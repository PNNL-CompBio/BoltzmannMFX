//
//     Copyright (c) 2013 Battelle Memorial Institute
//     Licensed under modified BSD License. A copy of this license can be found
//     in the LICENSE file in the top level directory of this distribution.
//
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
DepositionScheme bmx::m_cnc_deposition_scheme;
DepositionScheme bmx::m_vf_deposition_scheme;
amrex::Real      bmx::m_deposition_scale_factor = 1.0;

int  bmx::plot_int        = -1;
Real bmx::plot_per_approx = -1.;
Real bmx::plot_per_exact  = -1.;

// Destructor
bmx::~bmx ()
{
  if (DEM::solve)
    delete pc;

  for (int lev(0); lev <= finest_level; ++lev)
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

    SetNProper(2);

    /****************************************************************************
     *                                                                          *
     * Initialize time steps                                                    *
     *                                                                          *
     ***************************************************************************/

    t_old.resize(max_level+1,-1.e100);
    t_new.resize(max_level+1,0.0);


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

    fine_mask = 0;

    Gpu::synchronize();
}

void
bmx::ComputeAndPrintSums() 
{
  return;
    BL_PROFILE("bmx::ComputeAndPrintSums()");

    const auto p_lo = Geom(0).ProbLoArray();
    const auto p_hi = Geom(0).ProbHiArray();

    Real domain_vol = (p_hi[2]-p_lo[2])*(p_hi[1]-p_lo[1])*(p_hi[0]-p_lo[0]);

    Real fluid_vol = volSum();

    Real particle_vol = pc->computeParticleVolume();

    Real A_in_fluid     = volWgtSum(get_X_k_const(), 0);
    Real A_in_particles = pc->computeParticleContent(22);

    Real B_in_fluid     = volWgtSum(get_X_k_const(), 1);
    Real B_in_particles = pc->computeParticleContent(23);

    Real C_in_fluid     = volWgtSum(get_X_k_const(), 2);
    Real C_in_particles = pc->computeParticleContent(24);

    amrex::Print() << "Domain   volume : " << domain_vol << std::endl;
    amrex::Print() << "Particle + Fluid: " << fluid_vol+particle_vol << std::endl;
    amrex::Print() << "Fluid    volume : " << fluid_vol  << std::endl;
    amrex::Print() << "Particle volume : " << particle_vol << std::endl;

    if (std::abs(domain_vol - (fluid_vol+particle_vol)) > 1.e-12 * domain_vol)
       amrex::Abort("Volumes don't match!");

#if 1
    amrex::Print() << " A in fluid        : " << A_in_fluid << std::endl;
    amrex::Print() << " B in fluid        : " << B_in_fluid << std::endl;
    amrex::Print() << " C in fluid        : " << C_in_fluid << std::endl;

    amrex::Print() << " A in particles    : " << A_in_particles << std::endl;
    amrex::Print() << " B in particles    : " << B_in_particles << std::endl;
    amrex::Print() << " C in particles    : " << C_in_particles << std::endl;

    amrex::Print() << " A in fluid+part   : " << A_in_fluid + A_in_particles << std::endl;
    amrex::Print() << " C in fluid+part   : " << C_in_fluid + C_in_particles << std::endl;

    amrex::Print() << " A+C   in fluid    : " << A_in_fluid + C_in_fluid << std::endl;
    amrex::Print() << " A+C   in particles: " << A_in_particles + C_in_particles << std::endl;
#endif
    amrex::Print() << " Total A + C       : " << A_in_fluid + A_in_particles +
                                                 C_in_fluid + C_in_particles << std::endl;

}
void
bmx::compute_grad_X(int lev, Real time, MultiFab& gradx_X_k, MultiFab& grady_X_k, MultiFab& gradz_X_k) 
{
    fillpatch_Xk(get_X_k(), time);

    auto const dxInv = geom[lev].InvCellSizeArray();
 
    MultiFab* X_k = get_X_k()[lev];
    for (MFIter mfi(*X_k,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Box const& bx = mfi.tilebox();

        const Array4<const Real>  X_k_arr = X_k->const_array(mfi);
        const Array4<      Real> gx_k_arr = gradx_X_k.array(mfi);
        const Array4<      Real> gy_k_arr = grady_X_k.array(mfi);
        const Array4<      Real> gz_k_arr = gradz_X_k.array(mfi);

        ParallelFor(bx, FLUID::nchem_species, [=]
          AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
          {
              gx_k_arr(i,j,k,n) = (X_k_arr(i+1,j,k,n) - X_k_arr(i-1,j,k,n)) * 0.5 * dxInv[0];
              gy_k_arr(i,j,k,n) = (X_k_arr(i,j+1,k,n) - X_k_arr(i,j-1,k,n)) * 0.5 * dxInv[1];
              gz_k_arr(i,j,k,n) = (X_k_arr(i,j,k+1,n) - X_k_arr(i,j,k-1,n)) * 0.5 * dxInv[2];
          });
    }
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

Vector< MultiFab* > bmx::get_vf_old () noexcept
{
  Vector<MultiFab*> r;
  r.reserve(m_leveldata.size());
  for (int lev = 0; lev < m_leveldata.size(); ++lev) {
    r.push_back(m_leveldata[lev]->vf_o);
  }
  return r;
}

Vector< MultiFab* > bmx::get_vf_new () noexcept
{
  Vector<MultiFab*> r;
  r.reserve(m_leveldata.size());
  for (int lev = 0; lev < m_leveldata.size(); ++lev) {
    r.push_back(m_leveldata[lev]->vf_n);
  }
  return r;
}

Vector< MultiFab* > bmx::get_X_rhs () noexcept
{
  Vector<MultiFab*> r;
  r.reserve(m_leveldata.size());
  for (int lev = 0; lev < m_leveldata.size(); ++lev) {
    r.push_back(m_leveldata[lev]->X_rhs);
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

Vector< MultiFab const*> bmx::get_vf_old_const () const noexcept
{
  Vector<MultiFab const*> r;
  r.reserve(m_leveldata.size());
  for (int lev = 0; lev < m_leveldata.size(); ++lev) {
    r.push_back(m_leveldata[lev]->vf_o);
  }
  return r;
}

Vector< MultiFab const*> bmx::get_vf_new_const () const noexcept
{
  Vector<MultiFab const*> r;
  r.reserve(m_leveldata.size());
  for (int lev = 0; lev < m_leveldata.size(); ++lev) {
    r.push_back(m_leveldata[lev]->vf_n);
  }
  return r;
}
