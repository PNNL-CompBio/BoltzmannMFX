//
//     Copyright (c) 2013 Battelle Memorial Institute
//     Licensed under modified BSD License. A copy of this license can be found
//     in the LICENSE file in the top level directory of this distribution.
//
#include <bmx.H>
#include <bmx_fluid_parms.H>
#include <bmx_dem_parms.H>

void
bmx::Regrid ()
{
  if (ooo_debug) amrex::Print() << "Regrid" << std::endl;
  BL_PROFILE_REGION_START("bmx::Regrid()");

  int base_lev = 0;

  if (load_balance_type == "KnapSack" || load_balance_type == "SFC") // Knapsack and SFC
  {
    amrex::Print() << "Load balancing using " << load_balance_type << std::endl;

    if (DEM::solve)
       AMREX_ALWAYS_ASSERT(particle_cost[0] != nullptr);

    if (FLUID::solve)
       AMREX_ALWAYS_ASSERT(fluid_cost[0] != nullptr);

    if (ParallelDescriptor::NProcs() == 1) return;

    if (dual_grid)  //  Beginning of dual grid regridding
    {
      AMREX_ALWAYS_ASSERT(FLUID::solve);

      if (load_balance_fluid > 0)
      {
        for (int lev = base_lev; lev <= finestLevel(); ++lev)
        {
          DistributionMapping new_fluid_dm;

          if ( load_balance_type == "KnapSack" )
          {
            new_fluid_dm = DistributionMapping::makeKnapSack(*fluid_cost[lev],
                                                             knapsack_nmax);
          }
          else
          {
            new_fluid_dm = DistributionMapping::makeSFC(*fluid_cost[lev], false);
          }

          SetDistributionMap(lev, new_fluid_dm);

          RegridArrays(lev);

          if (fluid_cost[lev] != nullptr)
            delete fluid_cost[lev];

          fluid_cost[lev] = new MultiFab(grids[lev], new_fluid_dm, 1, 0);
          fluid_cost[lev]->setVal(0.0);
        }
      }

      for (int lev = base_lev; lev <= finestLevel(); ++lev)
      {
        DistributionMapping new_particle_dm;

        if ( load_balance_type == "KnapSack" )
        {
          new_particle_dm = DistributionMapping::makeKnapSack(*particle_cost[lev],
                                                              knapsack_nmax);
        }
        else
        {
          new_particle_dm = DistributionMapping::makeSFC(*particle_cost[lev],
                                                         false);
        }

        pc->Regrid(new_particle_dm, pc->ParticleBoxArray(lev), lev);

        if (particle_cost[lev] != nullptr)
          delete particle_cost[lev];

        particle_cost[lev] = new MultiFab(pc->ParticleBoxArray(lev),
                                                       new_particle_dm, 1, 0);
        particle_cost[lev]->setVal(0.0);
      }

    }
    else  // Single-grid regridding
    {
      MultiFab costs(grids[base_lev], dmap[base_lev], 1, 0);
      costs.setVal(0.0);

      Print() << "grids = " << grids[base_lev] << std::endl;
      Print() << "costs ba = " << costs.boxArray() << std::endl;

      if(DEM::solve)
        Print() << "particle_cost ba = "
                << particle_cost[base_lev]->boxArray()
                << std::endl;

      //Print() << "fluid cost ba = " << fluid_cost[base_lev]->boxArray() << std::endl;

      if (DEM::solve) {
        // costs.plus(* particle_cost[base_lev], 0, 1, 0);

        // MultiFab particle_cost_loc(grids[base_lev], dmap[base_lev], 1, 0);
        // particle_cost_loc.copy(* particle_cost[base_lev], 0, 0, 1);
        MultiFab particle_cost_loc = MFUtil::regrid(grids[base_lev], dmap[base_lev],
                                                    *particle_cost[base_lev],
                                                    true);

        costs.plus(particle_cost_loc, 0, 1, 0);
      }
      if (FLUID::solve) {
        // costs.plus(* fluid_cost[base_lev], 0, 1, 0);

        // MultiFab fluid_cost_loc(grids[base_lev], dmap[base_lev], 1, 0);
        // fluid_cost_loc.copy(* fluid_cost[base_lev], 0, 0, 1);
        MultiFab fluid_cost_loc = MFUtil::regrid(grids[base_lev], dmap[base_lev],
                                                 *fluid_cost[base_lev],
                                                 true);

        costs.plus(fluid_cost_loc, 0, 1, 0);
      }

      DistributionMapping newdm = DistributionMapping::makeKnapSack(costs,knapsack_nmax);

      SetDistributionMap(base_lev, newdm);

      if (FLUID::solve)
        RegridArrays(base_lev);

      if (FLUID::solve)
      {
        if (fluid_cost[base_lev] != nullptr)
          delete fluid_cost[base_lev];

        fluid_cost[base_lev] = new MultiFab(grids[base_lev], newdm, 1, 0);
        fluid_cost[base_lev]->setVal(0.0);
      }

      if (DEM::solve)
      {
        if (particle_cost[base_lev] != nullptr)
          delete particle_cost[base_lev];

        particle_cost[base_lev] = new MultiFab(grids[base_lev], newdm, 1, 0);
        particle_cost[base_lev]->setVal(0.0);
      }

      if (DEM::solve)
        pc->Regrid(dmap[base_lev], grids[base_lev], base_lev);

      }
  } else {
      amrex::Abort("load_balance_type must be KnapSack or SFC");
  }

  // This call resets the diffusion solvers
  if (FLUID::solve)
    bmx_setup_solvers();

  BL_PROFILE_REGION_STOP("bmx::Regrid()");
}
