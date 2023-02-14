//
//     Copyright (c) 2013 Battelle Memorial Institute
//     Licensed under modified BSD License. A copy of this license can be found
//     in the LICENSE file in the top level directory of this distribution.
//
#include <AMReX_MultiFab.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_VisMF.H>    // amrex::VisMF::Write(MultiFab)
#include <AMReX_VectorIO.H> // amrex::[read,write]IntData(array_of_ints)
#include <AMReX_AmrCore.H>
#include <AMReX_buildInfo.H>
#include <AMReX_Geometry.H>

#include <bmx.H>
#include <bmx_fluid_parms.H>
#include <bmx_dem_parms.H>

namespace
{
    const std::string level_prefix {"Level_"};
}

void
bmx::Restart (std::string& restart_file, int *nstep, Real *dt, Real *time)
{
    if (ooo_debug) amrex::Print() << "Restart" << std::endl;
    BL_PROFILE("bmx::Restart()");

    amrex::Print() << "  Restarting from checkpoint " << restart_file << std::endl;

    Real prob_lo[BL_SPACEDIM];
    Real prob_hi[BL_SPACEDIM];


    /***************************************************************************
     * Load header: set up problem domain (including BoxArray)                 *
     *              load particle data                                         *
     *              allocate bmx memory (bmx::AllocateArrays)                *
     ***************************************************************************/

    {
      std::string File(restart_file + "/Header");

      VisMF::IO_Buffer io_buffer(VisMF::GetIOBufferSize());

      Vector<char> fileCharPtr;
      ParallelDescriptor::ReadAndBcastFile(File, fileCharPtr);
      std::string fileCharPtrString(fileCharPtr.dataPtr());
      std::istringstream is(fileCharPtrString, std::istringstream::in);

      std::string line, word;

      std::getline(is, line);

      int  nlevs;
      int  int_tmp;
      Real real_tmp;

      is >> nlevs;
      GotoNextLine(is);

      // Time stepping controls
      is >> int_tmp;
      *nstep = int_tmp;
      GotoNextLine(is);

      is >> real_tmp;
      *dt = real_tmp;
      GotoNextLine(is);

      is >> real_tmp;
      *time = real_tmp;
      GotoNextLine(is);

        std::getline(is, line);
        {
            std::istringstream lis(line);
            int i = 0;
            while (lis >> word) {
               prob_lo[i++] = std::stod(word);
            }
        }

        std::getline(is, line);
        {
            std::istringstream lis(line);
            int i = 0;
            while (lis >> word) {
               prob_hi[i++] = std::stod(word);
            }
        }

        for (int lev = 0; lev < nlevs; ++lev) {

            RealBox rb(prob_lo,prob_hi);
            Geom(lev).ProbDomain(rb);
            Geom(lev).ResetDefaultProbDomain(rb);

            BoxArray ba;
            ba.readFrom(is);
            GotoNextLine(is);

            // Particle data is loaded into the BMXParticleContainer's base
            // class using amrex::NeighborParticleContainer::Restart

            if ( (DEM::solve) and lev == 0)
              pc->Restart(restart_file, "particles");

            amrex::Print() << "  Finished reading particle data" << std::endl;

            if (FLUID::solve) AllocateArrays(lev);
        }
    }

    amrex::Print() << "  Finished reading header" << std::endl;

    /***************************************************************************
     * Load fluid data                                                         *
     ***************************************************************************/
    if (FLUID::solve)
    {
       // Load the field data
       for (int lev = 0, nlevs=finestLevel()+1; lev < nlevs; ++lev)
       {
           // Read scalar variables
           ResetIOChkData();

           if (advect_fluid_chem_species)
           {
              for (int i = 0; i < chkChemSpeciesVars.size(); i++ )
              {
                 amrex::Print() << "  Loading " << chkChemSpeciesVarsName[i] << " at level " << lev << std::endl;
                 (chkChemSpeciesVars[i][lev])->setVal(0.);
    
                 MultiFab mf;
                 VisMF::Read(mf,
                         amrex::MultiFabFileFullPrefix(lev,
                                                   restart_file, level_prefix,
                                                   chkChemSpeciesVarsName[i]),
                                                   nullptr,
                                                   ParallelDescriptor::IOProcessorNumber());

                 // Copy from the mf we used to read in to the mf we will use going forward
                 (*(chkChemSpeciesVars[i][lev])).ParallelCopy(mf, 0, 0, FLUID::nchem_species,0,0);
              }
              {
                 amrex::Print() << "  Loading volume fraction " << " at level " << lev << std::endl;
    
                 MultiFab mf;
                 VisMF::Read(mf,
                         amrex::MultiFabFileFullPrefix(lev,
                                                   restart_file, level_prefix, "volfrac"),
                                                   nullptr,
                                                   ParallelDescriptor::IOProcessorNumber());

                 // Copy from the mf we used to read in to the mf we will use going forward
                 MultiFab* vf_n = m_leveldata[lev]->vf_n;
                 vf_n->setVal(0.);
                 vf_n->ParallelCopy(mf, 0, 0, 1, 0, 0);
              }
          }
       }

       amrex::Print() << "  Finished reading fluid data" << std::endl;
    }

    // Make sure that the particle BoxArray is the same as the mesh data -- we can
    //      create a dual grid decomposition in the regrid operation
    if (DEM::solve)
    {
        for (int lev = 0; lev <= finestLevel(); lev++)
        {
          pc->SetParticleBoxArray       (lev, grids[lev]);
          pc->SetParticleDistributionMap(lev,  dmap[lev]);
        }
        pc->Redistribute();

       // We need to do this on restart regardless of whether we replicate
       pc->Redistribute();
    }

    if (FLUID::solve)
    {
        for (int lev = 0; lev <= finestLevel(); lev++)
        {
          // Fill the bc's just in case
          m_leveldata[lev]->X_k->FillBoundary(geom[lev].periodicity());
          m_leveldata[lev]->D_k->FillBoundary(geom[lev].periodicity());
        }
    }

    if (load_balance_type == "KnapSack" or load_balance_type == "SFC")
    {
      if (DEM::solve) {
        for (int lev(0); lev < particle_cost.size(); ++lev)
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
      if (FLUID::solve) {
        for (int lev(0); lev < fluid_cost.size(); ++lev)
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
    amrex::Print() << "  Done with bmx::Restart " << std::endl;
}
