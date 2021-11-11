#include <bmx.H>
#include <bmx_fluid_parms.H>

// Location is just a marker used to determine at which location in the code
// this function is being called.
void bmx::print_mesh(int location)
{
#if 0
  const int nchem_species = FLUID::nchem_species;
  for (int lev = 1; lev <= finest_level; lev++)
  {
    auto& ld = *m_leveldata[lev];

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*ld.vf_n,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
      Box const& bx = mfi.tilebox();

      Array4<Real const> const& vf_n     = ld.vf_n->const_array(mfi);
      Array4<Real const> const& X_k_arr = ld.X_k->const_array(mfi);
      Array4<Real const> const& X_rhs_arr = ld.X_rhs->const_array(mfi);

      int ix, iy, iz;
      ix = -1;
      iy = -1;
      iz = -1;
      ParallelFor(bx, nchem_species, [&ix,&iy,&iz,vf_n]
          AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
          {
          if (vf_n(i,j,k) != 1.0)
          {
          ix = i;
          iy = j;
          iz = k;
          }
          });
      if (ix > -1 && iy > -1 && iz > -1) {
        std::cout << "VF("<<location<<") " << vf_n(ix,iy,iz,0) << std::endl;
        std::cout << "XK("<<location<<") " << X_k_arr(ix,iy,iz,0) << std::endl;
        std::cout << "XK("<<location<<") " << X_k_arr(ix,iy,iz,1) << std::endl;
        std::cout << "XK("<<location<<") " << X_k_arr(ix,iy,iz,2) << std::endl;
        std::cout << "RHS("<<location<<") " << X_rhs_arr(ix,iy,iz,0) << std::endl;
        std::cout << "RHS("<<location<<") " << X_rhs_arr(ix,iy,iz,1) << std::endl;
        std::cout << "RHS("<<location<<") " << X_rhs_arr(ix,iy,iz,2) << std::endl;
      }
    } // mfi
  } // lev
#endif
}

void bmx::check_mesh_values()
{
#ifndef AMREX_USE_GPU
  const int nchem_species = FLUID::nchem_species;
  for (int lev = 1; lev <= finest_level; lev++)
  {
    auto& ld = *m_leveldata[lev];

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*ld.vf_n,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
      Box const& bx = mfi.tilebox();

      Array4<Real const> const& X_k_arr = ld.X_k->const_array(mfi);

      ParallelFor(bx, nchem_species, [&X_k_arr]
          AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
          {
          if (X_k_arr(i,j,k,n) < 0.0)
          {
          printf("NEGATIVE Concentration X_k(%d,%d,%d,%d): %e\n",i,j,k,n,X_k_arr(i,j,k,n));
          }
          });
    } // mfi
  } // lev
#endif
}

void
bmx::EvolveFluid (int nstep,
                   Real& dt,
                   Real& /*prev_dt*/,
                   Real& time,
                   Real /*stop_time*/,
                   Real& /*coupling_timing*/)
{
    BL_PROFILE_REGION_START("bmx::EvolveFluid");
    BL_PROFILE("bmx::EvolveFluid");

    amrex::Print() << "\n ============   NEW TIME STEP   ============ \n";

    // Extrapolate boundary values 
    for (int lev = 0; lev <= finest_level; lev++)
    {
        m_leveldata[lev]->D_k->FillBoundary(geom[lev].periodicity());
        m_leveldata[lev]->X_k->FillBoundary(geom[lev].periodicity());
    }

    bmx_set_chem_species_bcs(time, get_X_k(), get_D_k());

    // Create temporary multifabs to hold the old-time lap_X and vel_RHS
    //    so we don't have to re-compute them in the corrector
    Vector< MultiFab* > lap_X(finest_level+1);
    Vector< MultiFab* > chem_species_RHS(finest_level+1);

    const int nchem_species = FLUID::nchem_species;

    for (int lev = 0; lev <= finest_level; lev++)
    {
       lap_X[lev]       = new MultiFab(grids[lev], dmap[lev], nchem_species, 0, MFInfo());
       chem_species_RHS[lev] = new MultiFab(grids[lev], dmap[lev], nchem_species, 0, MFInfo());
    }

    dt = fixed_dt;

    // Set new and old time to correctly use in fillpatching
    for (int lev = 0; lev <= finest_level; lev++)
    {
        t_old[lev] = time;
        t_new[lev] = time+dt;
    }

    amrex::Print() << "\n   Step " << nstep+1 << ": from old_time " \
                   << time << " to new time " << time+dt
                   << " with dt = " << dt << "\n" << std::endl;

    // Copy current "new" into "old"
    for (int lev = 0; lev <= finest_level; lev++)
    {
      MultiFab& X_k = *m_leveldata[lev]->X_k;
      MultiFab& X_ko = *m_leveldata[lev]->X_ko;
      MultiFab::Copy(X_ko, X_k, 0, 0, X_k.nComp(), X_ko.nGrow());

      MultiFab& vf_n = *m_leveldata[lev]->vf_n;
      MultiFab& vf_o = *m_leveldata[lev]->vf_o;
      MultiFab::Copy(vf_o, vf_n, 0, 0, 1, vf_n.nGrow());
    }

    // Interpolate chem_species to particle locations
    bmx_calc_txfr_particle(time, dt);

    // Deposit sources/sink from individual particles to grid
    bmx_calc_txfr_fluid(time, dt);
    print_mesh(1);
    //check_mesh_values();

    // Calculate the fraction of each grid cell not occupied by biological cells -- this
    //   1) defines vf_n using the current particle locations
    //   2) updates X_k on the grid to allow for the change in vf
    bmx_calc_volume_fraction();
    print_mesh(2);

    // Average down from fine to coarse to ensure consistency
    for (int lev = finest_level; lev > 0; lev--)
        average_down(*m_leveldata[lev]->vf_n,*m_leveldata[lev-1]->vf_n, 0, 1, refRatio(lev-1));

    print_mesh(3);
    //
    // Time integration step
    //
    Real new_time = time+dt;

    if (m_diff_type == DiffusionType::Implicit) amrex::Print() << "Doing fully implicit diffusion..." << std::endl;
    if (m_diff_type == DiffusionType::Explicit) amrex::Print() << "Doing fully explicit diffusion..." << std::endl;
    if (m_diff_type == DiffusionType::Crank_Nicolson) amrex::Print() << "Doing Crank-Nicolson diffusion..." << std::endl;

    // Local flag for explicit diffusion
    bool l_explicit_diff = (m_diff_type == DiffusionType::Explicit);

    fillpatch_Xk(get_X_k_old(), new_time);
    fillpatch_Dk(get_D_k(), new_time);
    print_mesh(4);

    // *************************************************************************************
    // Compute explicit diffusive terms
    // *************************************************************************************

    // fillpatch_vf(get_vf_new(),new_time);

    if (m_diff_type != DiffusionType::Implicit)
        diffusion_op->ComputeLapX(lap_X, get_X_k_old(), get_D_k_const(), get_vf_new_const());
    else
       for (int lev = 0; lev <= finest_level; lev++)
           lap_X[lev]->setVal(0.);
    print_mesh(5);
     

    // *************************************************************************************
    // Compute right hand side terms on the old status
    // *************************************************************************************

    // *************************************************************************
    // Update 
    // *************************************************************************
    Real l_dt = dt;

    // theta = 1 for fully explicit, theta = 1/2 for Crank-Nicolson
    Real theta;

    if (m_diff_type == DiffusionType::Explicit) theta = 1.0;
    if (m_diff_type == DiffusionType::Implicit) theta = 0.0;
    if (m_diff_type == DiffusionType::Crank_Nicolson) theta = 0.5;

    for (int lev = 0; lev <= finest_level; lev++)
    {
        Real grid_vol = (geom[lev].CellSize(0))*(geom[lev].CellSize(1))*(geom[lev].CellSize(2));
        auto& ld = *m_leveldata[lev];

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(*ld.X_k,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            Box const& bx = mfi.tilebox();
            Array4<Real      > const& X_k_n     = ld.X_k->array(mfi);
            Array4<Real const> const& X_RHS_arr = ld.X_rhs->const_array(mfi);
            Array4<Real const> const& lap_X_arr = lap_X[lev]->const_array(mfi);
            Array4<Real const> const&    vf_arr = ld.vf_n->const_array(mfi);

            ParallelFor(bx, nchem_species, [=]
              AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
                X_k_n(i,j,k,n) += theta * l_dt * lap_X_arr(i,j,k,n) / vf_arr(i,j,k) 
                                               + X_RHS_arr(i,j,k,n) / (vf_arr(i,j,k) * grid_vol);
                if (X_k_n(i,j,k) < 0.0) amrex::Print() << " LOW VAL OLD/NEW/LAP/RHS " << IntVect(i,j,k) << " " << n << " " <<
                    X_k_n(i,j,k,n) - (theta * l_dt * lap_X_arr(i,j,k,n) / vf_arr(i,j,k) 
                                              + X_RHS_arr(i,j,k,n) / (vf_arr(i,j,k) * grid_vol)) << " " << 
                      X_k_n(i,j,k,n) << " " << l_dt * lap_X_arr(i,j,k,n) / vf_arr(i,j,k)  << " " <<
                                               + X_RHS_arr(i,j,k,n) / (vf_arr(i,j,k) * grid_vol) << std::endl;

            });
        } // mfi
    } // lev

    // Average down from fine to coarse to ensure consistency
    for (int lev = finest_level; lev > 0; lev--)
        average_down(*m_leveldata[lev]->X_k,*m_leveldata[lev-1]->X_k, 0, nchem_species, refRatio(lev-1));

    for (int lev = finest_level; lev >= 0; lev--)
    {
        auto& ld = *m_leveldata[lev];

        amrex::Real x0max = ld.X_k->max(0);
        amrex::Real x1max = ld.X_k->max(1);
        amrex::Real x2max = ld.X_k->max(2);
        amrex::Real x0min = ld.X_k->min(0);
        amrex::Real x1min = ld.X_k->min(1);
        amrex::Real x2min = ld.X_k->min(2);
        amrex::Print() << "Max/min of species 0 at level " << lev << " " << x0max << " " << x0min << std::endl;
        amrex::Print() << "Max/min of species 1 at level " << lev << " " << x1max << " " << x1min << std::endl;
        amrex::Print() << "Max/min of species 2 at level " << lev << " " << x2max << " " << x2min << std::endl;

        Real eps = 1.e-8;
        if (x0max > 1.0+eps || x1max > 1.0+eps || x2max > 1.0+eps) amrex::Abort("Species greater than 1");
        if (x0min < 0.0-eps || x1min < 0.0-eps || x2min < 0.0-eps) amrex::Abort("Species    less than 0");
    } // lev

    print_mesh(6);


    // *************************************************************************************
    // If doing implicit diffusion...
    // *************************************************************************************
    if (not l_explicit_diff) 
    {
        bmx_set_chem_species_bcs(time, get_X_k(), get_D_k());

        Real omt = 1. - theta;
        diffusion_op->diffuse_chem_species(get_X_k(), get_D_k_const(), get_vf_new_const(), omt, l_dt);

        for (int lev = 0; lev <= finest_level; lev++)
        {
            Real grid_vol = (geom[lev].CellSize(0))*(geom[lev].CellSize(1))*(geom[lev].CellSize(2));
            auto& ld = *m_leveldata[lev];

            for (MFIter mfi(*ld.X_k,TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                Box const& bx = mfi.tilebox();
                Array4<Real      > const& X_k_o     = ld.X_ko->array(mfi);
                Array4<Real      > const& X_k_n     = ld.X_k->array(mfi);
                Array4<Real const> const& X_RHS_arr = ld.X_rhs->const_array(mfi);

                ParallelFor(bx, nchem_species, [=]
                  AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    if (X_k_n(i,j,k,n) < -1.e-12) 
                    {
                        std::cout << "Implicitly created negative fluid X_k at (i,j,k) " 
                                  << IntVect(i,j,k) << " at level  " << lev << " in component " << n << " " << X_k_n(i,j,k) << std::endl;
                        std::cout << "OLD WAS " << X_k_o(i,j,k,n) << std::endl;
                        std::cout << "RHS WAS " << X_RHS_arr(i,j,k,n) << std::endl;
                        amrex::Abort();
                    } else {
                        X_k_n(i,j,k,n) = std::max(X_k_n(i,j,k,n),0.0);
                    }
                });
            }
        }
    }
    print_mesh(7);

    for (int lev = 0; lev <= finest_level; lev++)
    {
       delete lap_X[lev];
       delete chem_species_RHS[lev];
    }

    BL_PROFILE_REGION_STOP("bmx::EvolveFluid");
}
