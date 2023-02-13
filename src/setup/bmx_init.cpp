//
//     Copyright (c) 2013 Battelle Memorial Institute
//     Licensed under modified BSD License. A copy of this license can be found
//     in the LICENSE file in the top level directory of this distribution.
//
#include <AMReX_ParmParse.H>

#include <bmx.H>
#include <bmx_init_fluid.H>
#include <bmx_bc_parms.H>
#include <bmx_dem_parms.H>
#include <bmx_fluid_parms.H>
#include <bmx_calc_fluid_coeffs.H>
#include <bmx_chem_species_parms.H>

using BMXParIter = BMXParticleContainer::BMXParIter;
using PairIndex = BMXParticleContainer::PairIndex;

void
bmx::InitParams ()
{
  if (ooo_debug) amrex::Print() << "InitParams" << std::endl;

  // Read and process chem_species, fluid and DEM particle model options.
  SPECIES::Initialize();
  FLUID::Initialize();

  BL_ASSERT(FLUID::nchem_species <= SPECIES::NMAX);

  DEM::Initialize();

  // Important! Resize the bc vector for the fluid chem_species mass fractions
  // We have to do it here because the size has to match the number of fluid
  // chem_species
  // NOTE: once we will have a class for BCs this won't be needed anymore
  m_bc_X_k.resize(FLUID::nchem_species, Gpu::DeviceVector<Real>(50, 0));
  m_bc_X_k_ptr.resize(FLUID::nchem_species, nullptr);
  {
      Vector<Real*> tmp(FLUID::nchem_species);
      for (int i = 0; i < FLUID::nchem_species; ++i) {
          tmp[i] = m_bc_X_k[i].data();
      }
      Gpu::copyAsync(Gpu::hostToDevice, tmp.begin(), tmp.end(), m_bc_X_k_ptr.begin());
      Gpu::synchronize();
  }
  bcs_X.resize(2*FLUID::nchem_species);
  bcs_D.resize(2*FLUID::nchem_species);

  BC::Initialize(geom[0]);

  {
    ParmParse pp("bmx");

    // Options to control time stepping
    pp.query("cfl", m_cfl);

    fixed_dt = -1.;
    pp.query("fixed_dt", fixed_dt);
    pp.query("dt_min", dt_min);
    pp.query("dt_max", dt_max);

    // Flag to set verbosity
    m_verbose = 0;
    pp.query("verbose", m_verbose);

    pp.query("ooo_debug", ooo_debug);

    // Initialize random number generator
    int seed = 77389;
    pp.query("seed", seed);
    amrex::ResetRandomSeed(seed+ParallelDescriptor::MyProc()+1);

    // The default type is "AsciiFile" but we can over-write that in the inputs
    // file with "Random"
    pp.query("particle_init_type", particle_init_type);

    Array<int,3> sorting_bin{0, 0, 0};
    pp.query("particle_sorting_bin", sorting_bin);
    particle_sorting_bin = IntVect(sorting_bin);

    // Set the bmx class flag equal to the FLUID parameter
    advect_fluid_chem_species = FLUID::solve_chem_species;

    // Read in number of substeps in chemistry integration
    m_nloop = 4;
    pp.query("substeps",m_nloop);
    amrex::Print() << "SUBSTEPS: " <<m_nloop<<std::endl;

    // We can still turn it off explicitly even if we passed chem_species inputs
    pp.query("advect_fluid_chem_species", advect_fluid_chem_species);

    if (advect_fluid_chem_species)
      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(FLUID::solve_chem_species,
          "Advect fluid chem_species flag is on but no fluid chem_species were provided");

    // The default type is "KnapSack"; alternative is "SFC"
    pp.query("load_balance_type", load_balance_type);
    pp.query("knapsack_weight_type", knapsack_weight_type);
    pp.query("load_balance_fluid", load_balance_fluid);

    // The default for diffusion_type is 2, i.e. the default m_diff_type is DiffusionType::Implicit
    int diffusion_type = 1;
    pp.query("diffusion_type", diffusion_type);
    if (diffusion_type == 0) {
      m_diff_type = DiffusionType::Explicit;
    } else if (diffusion_type == 1) {
      m_diff_type = DiffusionType::Crank_Nicolson;
    } else if (diffusion_type == 2) {
      m_diff_type = DiffusionType::Implicit;
    } else if (diffusion_type == 3) {
      m_diff_type = DiffusionType::Explicit;
    } else {
      amrex::Abort("We currently require diffusion_type be one of the following:\
                   \n   0: explicit,\
                   \n   1: Crank-Nicolson\
                   \n   2: implicit\
                   \n   3: explicit predictor");
    }

    AMREX_ALWAYS_ASSERT(load_balance_type.compare("KnapSack") == 0  or
                        load_balance_type.compare("SFC") == 0);

    AMREX_ALWAYS_ASSERT(knapsack_weight_type.compare("RunTimeCosts") == 0 or
                        knapsack_weight_type.compare("NumParticles") == 0);

    ParmParse amr_pp("amr");
    amr_pp.query("dual_grid", dual_grid);

    if (load_balance_type.compare("KnapSack") == 0)
      pp.query("knapsack_nmax", knapsack_nmax);
  }

  if (DEM::solve)
  {
    ParmParse pp("particles");

    pp.query("max_grid_size_x", particle_max_grid_size_x);
    pp.query("max_grid_size_y", particle_max_grid_size_y);
    pp.query("max_grid_size_z", particle_max_grid_size_z);

    // Keep particles that are initially touching the wall. Used by DEM tests.
    pp.query("removeOutOfRange", removeOutOfRange);
  }

  if ((DEM::solve) and (not FLUID::solve))
  {
    if (fixed_dt <= 0.0)
      amrex::Abort("If running particle-only must specify a positive fixed_dt"
          " in the inputs file");
  }

  if ((DEM::solve) and FLUID::solve)
  {
    ParmParse pp("bmx");

    std::string cnc_deposition_scheme = "one_to_one";
    pp.query("cnc_deposition_scheme", cnc_deposition_scheme);

    if (cnc_deposition_scheme.compare("one_to_one") == 0) {
      m_cnc_deposition_scheme = DepositionScheme::one_to_one;
    }
    else if (cnc_deposition_scheme.compare("trilinear") == 0) {
      m_cnc_deposition_scheme = DepositionScheme::trilinear;
    }
    else if (cnc_deposition_scheme.compare("trilinear-dpvm-square") == 0) {
      m_cnc_deposition_scheme = DepositionScheme::square_dpvm;
    }
    else if (cnc_deposition_scheme.compare("true-dpvm") == 0) {
      m_cnc_deposition_scheme = DepositionScheme::true_dpvm;
    }
    else if (cnc_deposition_scheme.compare("centroid") == 0) {
      m_cnc_deposition_scheme = DepositionScheme::centroid;
    }
    else {
      amrex::Abort("Don't know this deposition_scheme for concentrations!");
    }

    std::string vf_deposition_scheme = "one_to_one";
    pp.query("vf_deposition_scheme", vf_deposition_scheme);

    if (vf_deposition_scheme.compare("one_to_one") == 0) {
      m_vf_deposition_scheme = DepositionScheme::one_to_one;
    }
    else if (vf_deposition_scheme.compare("trilinear") == 0) {
      m_vf_deposition_scheme = DepositionScheme::trilinear;
    }
    else if (vf_deposition_scheme.compare("trilinear-dpvm-square") == 0) {
      m_vf_deposition_scheme = DepositionScheme::square_dpvm;
    }
    else if (vf_deposition_scheme.compare("true-dpvm") == 0) {
      m_vf_deposition_scheme = DepositionScheme::true_dpvm;
    }
    else if (vf_deposition_scheme.compare("centroid") == 0) {
      m_vf_deposition_scheme = DepositionScheme::centroid;
    }
    else {
      amrex::Abort("Don't know this deposition_scheme for volume fraction!");
    }

    m_deposition_scale_factor = 1.;
    pp.query("deposition_scale_factor", m_deposition_scale_factor);
  }

  {
    ParmParse amr_pp("amr");

    amr_pp.query("restart_from_cold_flow", restart_from_cold_flow);

    amr_pp.query("plot_int", plot_int);
    amr_pp.query("plot_per_exact", plot_per_exact);
    amr_pp.query("plot_per_approx", plot_per_approx);

    if ((plot_int       > 0 and plot_per_exact  > 0) or
        (plot_int       > 0 and plot_per_approx > 0) or
        (plot_per_exact > 0 and plot_per_approx > 0) )
      amrex::Abort("Must choose only one of plot_int or plot_per_exact or plot_per_approx");
  }
}

void bmx::Init (Real time)
{
    if (ooo_debug) amrex::Print() << "Init" << std::endl;
    InitIOChkData();
    InitIOPltData();

    /****************************************************************************
     *                                                                          *
     * Generate levels using ErrorEst tagging.                                  *
     *                                                                          *
     ***************************************************************************/

    // This tells the AmrMesh class not to iterate when creating the initial
    // grid hierarchy
    SetIterateToFalse();

    // This tells the Cluster routine to use the new chopping routine which
    // rejects cuts if they don't improve the efficiency
    SetUseNewChop();

    /****************************************************************************
     *                                                                          *
     * BMX-specific grid creation                                              *
     *                                                                          *
     ***************************************************************************/

    // Define coarse level BoxArray and DistributionMap
    // This is an AmrCore member function which recursively makes new levels
    // with MakeNewLevelFromScratch.
    InitFromScratch(time);

    for (int lev = 1; lev <= finestLevel(); lev++)
    {
       if (m_verbose > 0)
            std::cout << "Setting refined region at level " << lev
                      << " to " << grids[lev] << std::endl;

       MakeNewLevelFromScratch(lev, time, grids[lev], dmap[lev]);
    }

    MakeBCArrays();

    /****************************************************************************
     *                                                                          *
     * Create particle container using bmx::ParGDB                             *
     *                                                                          *
     ***************************************************************************/

    if (DEM::solve) {
      pc = new BMXParticleContainer(this);
    }

    /****************************************************************************
     *                                                                          *
     * BMX-Specific Initialization                                             *
     *                                                                          *
     ***************************************************************************/

    // ******************************************************
    // We only do these at level 0
    // ******************************************************

    for (int lev = 0; lev <= finestLevel(); lev++)
        bmx_set_bc_type(lev);
}


BoxArray bmx::MakeBaseGrids () const
{
    if (ooo_debug) amrex::Print() << "MakeBaseGrids" << std::endl;
    BoxArray ba(geom[0].Domain());

    ba.maxSize(max_grid_size[0]);

    // We only call ChopGrids if dividing up the grid using max_grid_size didn't
    //    create enough grids to have at least one grid per processor.
    // This option is controlled by "refine_grid_layout" which defaults to true.
    if ( refine_grid_layout &&
         ba.size() < ParallelDescriptor::NProcs() )
           ChopGrids(geom[0].Domain(), ba, ParallelDescriptor::NProcs());

    if (ba == grids[0]) {
        ba = grids[0];  // to avoid duplicates
    }
    amrex::Print() << "In MakeBaseGrids: BA HAS " << ba.size() << " GRIDS " << std::endl;
    return ba;
}


void bmx::ChopGrids (const Box& domain, BoxArray& ba, int target_size) const
{
    if (ooo_debug) amrex::Print() << "ChopGrids" << std::endl;
    if ( ParallelDescriptor::IOProcessor() )
       amrex::Warning("Using max_grid_size only did not make enough grids for the number of processors");

    // Here we hard-wire the maximum number of times we divide the boxes.
    int n = 10;

    // Here we hard-wire the minimum size in any one direction the boxes can be
    int min_grid_size = 4;

    IntVect chunk(domain.length(0),domain.length(1),domain.length(2));

    int j = -1;
    for (int cnt = 1; cnt <= n; ++cnt)
    {
        if (chunk[0] >= chunk[1] && chunk[0] >= chunk[2])
        {
            j = 0;
        }
        else if (chunk[1] >= chunk[0] && chunk[1] >= chunk[2])
        {
            j = 1;
        }
        else if (chunk[2] >= chunk[0] && chunk[2] >= chunk[1])
        {
            j = 2;
        }
        chunk[j] /= 2;

        if (chunk[j] >= min_grid_size)
        {
            ba.maxSize(chunk);
        }
        else
        {
            // chunk[j] was the biggest chunk -- if this is too small then we're done
            if ( ParallelDescriptor::IOProcessor() )
               amrex::Warning("ChopGrids was unable to make enough grids for the number of processors");
            return;
        }

        // Test if we now have enough grids
        if (ba.size() >= target_size) return;
    }
}


void bmx::MakeNewLevelFromScratch (int lev, Real /*time*/,
                                    const BoxArray& new_grids,
                                    const DistributionMapping& new_dmap)
{
    if (ooo_debug) amrex::Print() << "MakeNewLevelFromScratch" << std::endl;
    if (m_verbose > 0)
    {
        std::cout << "MAKING NEW LEVEL " << lev << std::endl;
        std::cout << "WITH BOX ARRAY   " << new_grids << std::endl;
    }

    SetBoxArray(lev, new_grids);
    SetDistributionMap(lev, new_dmap);
    amrex::Print() << "SETTING NEW GRIDS IN MAKE NEW LEVEL " << new_grids << std::endl;
    amrex::Print() << "SETTING NEW DMAP IN MAKE NEW LEVEL " << new_dmap << std::endl;

    // This is being done by bmx::make_eb_geometry,
    // otherwise it would be done here
    if (lev == 0) MakeBCArrays();
}


void bmx::ReMakeNewLevelFromScratch (int lev,
                                     const BoxArray & new_grids,
                                     const DistributionMapping & new_dmap)
{
    if (ooo_debug) amrex::Print() << "ReMakeNewLevelFromScratch" << std::endl;
    SetBoxArray(lev, new_grids);
    SetDistributionMap(lev, new_dmap);

    if (lev == 0) MakeBCArrays();

    // We need to re-fill these arrays for the larger domain (after replication).
    bmx_set_bc_type(lev);
}

void bmx::InitLevelData (Real /*time*/)
{
    if (ooo_debug) amrex::Print() << "InitLevelData" << std::endl;

    // Allocate the fluid data
    if (FLUID::solve)
       for (int lev = 0; lev <= finestLevel(); lev++)
          AllocateArrays(lev);

    ParmParse pp("particles");
    // Allocate the particle data
    if (DEM::solve)
    {
      Real strt_init_part = ParallelDescriptor::second();

      pc->AllocData();

      if (particle_init_type == "AsciiFile")
      {
        std::string filename;
        pp.query("input_file",filename);
        amrex::Print() << "Reading particles from "<<filename<<" ..." << std::endl;
        pc->InitParticlesAscii(filename);

      } else { 

         amrex::Abort("Bad particle_init_type");
      }

      pc->Redistribute();

      Real end_init_part = ParallelDescriptor::second() - strt_init_part;
      ParallelDescriptor::ReduceRealMax(end_init_part, ParallelDescriptor::IOProcessorNumber());
      amrex::Print() << "Time spent in initializing particles " << end_init_part << std::endl;
    }

    // Used in load balancing
    if (DEM::solve)
    {
      for (int lev(0); lev < particle_cost.size(); lev++)
        if (particle_cost[lev] != nullptr)
          delete particle_cost[lev];

      particle_cost.clear();
      particle_cost.resize(finestLevel()+1, nullptr);

      for (int lev = 0; lev <= finestLevel(); lev++)
      {
        particle_cost[lev] = new MultiFab(pc->ParticleBoxArray(lev),
                                          pc->ParticleDistributionMap(lev), 1, 0);
        particle_cost[lev]->setVal(0.0);
      }
    }

    // Used in load balancing
    if (FLUID::solve)
    {
      for (int lev(0); lev < fluid_cost.size(); lev++)
        if (fluid_cost[lev] != nullptr)
          delete fluid_cost[lev];

      fluid_cost.clear();
      fluid_cost.resize(finestLevel()+1, nullptr);

      for (int lev = 0; lev <= finestLevel(); lev++)
      {
        fluid_cost[lev] = new MultiFab(grids[lev], dmap[lev], 1, 0);
        fluid_cost[lev]->setVal(0.0);
      }
    }
}

void
bmx::PostInit (Real& dt, Real /*time*/, int restart_flag, Real stop_time)
{
    if (ooo_debug) amrex::Print() << "PostInit" << std::endl;

    if (DEM::solve)
    {
        // Auto generated particles may be out of the domain. This call will
        // remove them. 

        // We need to do this *after* restart (hence putting this here not
        // in Init) because we may want to change the particle_max_grid_size on restart.
        if ( dual_grid && particle_max_grid_size_x > 0
                       && particle_max_grid_size_y > 0
                       && particle_max_grid_size_z > 0)
        {
          IntVect particle_max_grid_size(particle_max_grid_size_x,
                                         particle_max_grid_size_y,
                                         particle_max_grid_size_z);

          for (int lev = 0; lev <= finestLevel(); lev++)
          {
            BoxArray particle_ba(geom[lev].Domain());
            particle_ba.maxSize(particle_max_grid_size);
            DistributionMapping particle_dm(particle_ba, ParallelDescriptor::NProcs());
            pc->Regrid(particle_dm, particle_ba);

            if (particle_cost[lev] != nullptr)
              delete particle_cost[lev];

            particle_cost[lev] = new MultiFab(pc->ParticleBoxArray(lev),
                                              pc->ParticleDistributionMap(lev), 1, 0);
            particle_cost[lev]->setVal(0.0);

          }
        }

        if (!FLUID::solve){
            dt = fixed_dt;
        }
    }

    if (FLUID::solve)
        bmx_init_fluid(restart_flag, dt, stop_time);
}

void
bmx::MakeBCArrays ()
{
    for (int lev = 0; lev < bc_ilo.size(); lev++)
    {
      if (bc_ilo[lev] != nullptr) delete bc_ilo[lev];
      if (bc_ihi[lev] != nullptr) delete bc_ihi[lev];
      if (bc_jlo[lev] != nullptr) delete bc_jlo[lev];
      if (bc_jhi[lev] != nullptr) delete bc_jhi[lev];
      if (bc_klo[lev] != nullptr) delete bc_klo[lev];
      if (bc_khi[lev] != nullptr) delete bc_khi[lev];
    }

    if (ooo_debug) amrex::Print() << "MakeBCArrays with finest level " << finestLevel() << std::endl;
    bc_ilo.clear(); bc_ilo.resize(finestLevel()+1, nullptr);
    bc_ihi.clear(); bc_ihi.resize(finestLevel()+1, nullptr);
    bc_jlo.clear(); bc_jlo.resize(finestLevel()+1, nullptr);
    bc_jhi.clear(); bc_jhi.resize(finestLevel()+1, nullptr);
    bc_klo.clear(); bc_klo.resize(finestLevel()+1, nullptr);
    bc_khi.clear(); bc_khi.resize(finestLevel()+1, nullptr);

    for (int lev = 0; lev <= finestLevel(); lev++)
    {
       // Define and allocate the integer MultiFab that is the outside adjacent
       // cells of the problem domain.
       Box domainx(geom[lev].Domain());
       domainx.grow(1,nghost);
       domainx.grow(2,nghost);
       Box box_ilo = amrex::adjCellLo(domainx,0,1);
       Box box_ihi = amrex::adjCellHi(domainx,0,1);

       Box domainy(geom[lev].Domain());
       domainy.grow(0,nghost);
       domainy.grow(2,nghost);
       Box box_jlo = amrex::adjCellLo(domainy,1,1);
       Box box_jhi = amrex::adjCellHi(domainy,1,1);

       Box domainz(geom[lev].Domain());
       domainz.grow(0,nghost);
       domainz.grow(1,nghost);
       Box box_klo = amrex::adjCellLo(domainz,2,1);
       Box box_khi = amrex::adjCellHi(domainz,2,1);

       // Note that each of these is a single IArrayBox so every process has a copy of them
       bc_ilo[lev] = new IArrayBox(box_ilo,2);
       bc_ihi[lev] = new IArrayBox(box_ihi,2);
       bc_jlo[lev] = new IArrayBox(box_jlo,2);
       bc_jhi[lev] = new IArrayBox(box_jhi,2);
       bc_klo[lev] = new IArrayBox(box_klo,2);
       bc_khi[lev] = new IArrayBox(box_khi,2);
   }
}

void
bmx::bmx_init_fluid (int is_restarting, Real /*dt*/, Real /*stop_time*/)
{
    if (ooo_debug) amrex::Print() << "bmx_init_fluid" << std::endl;

    // Set to bogus value just to make sure everything gets filled later 
    if (!is_restarting)
    {
        for (int lev = 0; lev <= finestLevel(); lev++)
        {
             m_leveldata[lev]->X_k->setVal(1.e234);
             m_leveldata[lev]->D_k->setVal(1.e234);
        }
    }

    for (int lev = 0; lev <= finestLevel(); lev++)
    {
       Box domain(geom[lev].Domain());

       Real dx = geom[lev].CellSize(0);
       Real dy = geom[lev].CellSize(1);
       Real dz = geom[lev].CellSize(2);

       const GpuArray<Real, 3> p_lo = geom[lev].ProbLoArray();
       const GpuArray<Real, 3> p_hi = geom[lev].ProbHiArray();

       LevelData& ld = *m_leveldata[lev];

       // We deliberately don't tile this loop since we will be looping
       //    over bc's on faces and it makes more sense to do this one grid at a time
       MultiFab dummy(grids[lev],dmap[lev],1,0);
       for (MFIter mfi(dummy, false); mfi.isValid(); ++mfi) 
       {
          const Box& bx = mfi.validbox();
          const Box& sbx = dummy[mfi].box();

          Array4<Real> const& X_k_arr = ld.X_k->array(mfi);
          Array4<Real> const& D_k_arr = ld.D_k->array(mfi);

          // Set the initial fluid chem_species mass fractions
          if (advect_fluid_chem_species) 
          {
              if (!is_restarting) 
                  set_ic_chem_species(sbx, domain, dx, dy, dz, p_lo, p_hi, X_k_arr);

             calc_D_k(bx, domain, dx, dy, dz, D_k_arr);
          }
       }

       MultiFab::Copy(*m_leveldata[lev]->X_ko, *m_leveldata[lev]->X_k, 0, 0,
             m_leveldata[lev]->X_k->nComp(), 0);

       m_leveldata[lev]->X_k->FillBoundary(geom[lev].periodicity());
       m_leveldata[lev]->D_k->FillBoundary(geom[lev].periodicity());
    }

    if (is_restarting == 0)
    {
      Real time = 0.0;

      if (advect_fluid_chem_species) {
        bmx_set_chem_species_bcs(time, get_X_k(), get_D_k());
        bmx_set_chem_species_bcs(time, get_X_k_old(), get_D_k());
      }
    }

    // Calculate the initial volume fraction
    if (!is_restarting) 
    {
        bool adjust_X = false;
        bmx_calc_volume_fraction(adjust_X);
        for (int lev = 0; lev <= finestLevel(); lev++)
           MultiFab::Copy(*m_leveldata[lev]->vf_o, *m_leveldata[lev]->vf_n, 0, 0, 1, m_leveldata[lev]->vf_n->nGrow());
    }

    // Average down from fine to coarse to ensure consistency
    for (int lev = finestLevel(); lev > 0; lev--)
    {
        average_down(*m_leveldata[lev]->vf_n,*m_leveldata[lev-1]->vf_n, 0, 1, refRatio(lev-1));
        if (!is_restarting) 
            average_down(*m_leveldata[lev]->vf_o,*m_leveldata[lev-1]->vf_o, 0, 1, refRatio(lev-1));
    }
}
