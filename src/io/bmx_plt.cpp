#include <AMReX_PlotFileUtil.H>
#include <AMReX_VisMF.H>    // amrex::VisMF::Write(MultiFab)
#include <AMReX_VectorIO.H> // amrex::[read,write]IntData(array_of_ints)
#include <AMReX_ParmParse.H>

#include <bmx.H>
#include <bmx_fluid_parms.H>
#include <bmx_dem_parms.H>

namespace
{
    const std::string level_prefix {"Level_"};
}

void
bmx::InitIOPltData ()
{
    if (ooo_debug) amrex::Print() << "InitIOPltData" << std::endl;

    // Variables to simplify checkpoint IO
    pltVarCount = 0;

    if (FLUID::solve)
    {
        ParmParse pp("amr");

        pp.query("plt_X",     plt_X_gk   );
        pp.query("plt_D",     plt_D_gk   );
 
        if( plt_X_gk == 1)  pltVarCount += FLUID::nspecies;
        if( plt_D_gk == 1)  pltVarCount += FLUID::nspecies;
    }
}

void
bmx::WritePlotFile (std::string& plot_file, int nstep, Real time )
{
    // If we've already written this plotfile, don't do it again!
    if (nstep == last_plt) return;

    // Now set last_plt to nstep ...
    last_plt = nstep;

    BL_PROFILE("bmx::WritePlotFile()");

    const std::string& plotfilename = amrex::Concatenate(plot_file,nstep);

    amrex::Print() << "  Writing plotfile " << plotfilename <<  " at time " << time << std::endl;

    if (pltVarCount > 0) {

      const int ngrow = 0;

      Vector<std::string> pltFldNames;
      Vector< std::unique_ptr<MultiFab> > mf(nlev);

      // Fluid species mass fractions
      if(FLUID::solve_species and plt_X_gk == 1)
        for(std::string specie: FLUID::species)
          pltFldNames.push_back("X_"+specie);

      // Fluid species mass diffusivities
      if(FLUID::solve_species and plt_D_gk == 1)
        for(std::string specie: FLUID::species)
          pltFldNames.push_back("D_"+specie);

      for (int lev = 0; lev < nlev; ++lev)
      {
        // Multifab to hold all the variables -- there can be only one!!!!
        const int ncomp = pltVarCount;
        mf[lev].reset(new MultiFab(grids[lev], dmap[lev], ncomp, ngrow,  MFInfo()));

        int lc=0;

        if(FLUID::solve_species and plt_X_gk == 1 ) {
          for(int n(0); n < FLUID::nspecies; n++) {
            MultiFab::Copy(*mf[lev], *m_leveldata[lev]->X_gk, n, lc+n, 1, 0);
          }
          lc += FLUID::nspecies;
        }

        // Species mass fraction
        if(FLUID::solve_species and plt_D_gk == 1 ) {
          for(int n(0); n < FLUID::nspecies; n++) {
            MultiFab::Copy(*mf[lev], *m_leveldata[lev]->D_gk, n, lc+n, 1, 0);
          }
          lc += FLUID::nspecies;
        }
      }

      Vector<const MultiFab*> mf2(nlev);
      for (int lev = 0; lev < nlev; ++lev) {
        mf2[lev] = mf[lev].get();
      }

      Vector<int> istep;
      istep.resize(nlev,nstep);
      amrex::WriteMultiLevelPlotfile(plotfilename, nlev, mf2, pltFldNames,
                                     Geom(), time, istep, refRatio());

      // no fluid
    } else {

      // Some post-processing tools (such as yt) might still need some basic
      // MultiFab header information to function. We provide this here by
      // creating an "empty" plotfile header (which essentially only contains
      // the BoxArray information). Particle data is saved elsewhere.

      Vector< std::unique_ptr<MultiFab> > mf(finest_level+1);
      Vector<std::string>  names;
      // NOTE: leave names vector empty => header should reflect nComp = 0
      //names.insert(names.end(), "placeholder");

      // Create empty MultiFab containing the right BoxArray (NOTE: setting
      // nComp = 1 here to avoid assertion fail in debug build).
      for (int lev = 0; lev <= finest_level; ++lev)
        mf[lev].reset(new MultiFab(grids[lev], dmap[lev], 1, 0));

      Vector<const MultiFab*> mf2(finest_level+1);

      for (int lev = 0; lev <= finest_level; ++lev)
        mf2[lev] = mf[lev].get();

      // Write only the Headers corresponding to the "empty" mf/mf2 MultiFabs
      Vector<int> istep;
      istep.resize(nlev,nstep);
      amrex::WriteMultiLevelPlotfileHeaders(plotfilename, finest_level+1, mf2, names,
                                            Geom(), time, istep, refRatio());

    }

    WriteJobInfo(plotfilename);

    if ( DEM::solve )
    {
        Vector<std::string> real_comp_names;
        Vector<std::string>  int_comp_names;

        real_comp_names.push_back("velx");
        real_comp_names.push_back("vely");
        real_comp_names.push_back("velz");
        real_comp_names.push_back("radius");
        real_comp_names.push_back("volume");
        real_comp_names.push_back("fluid_A");
        real_comp_names.push_back("fluid_B");
        real_comp_names.push_back("consume_A");
        real_comp_names.push_back("consume_B");

        int_comp_names.push_back("phase");
        int_comp_names.push_back("state");

        Vector<int> write_int_comp  = Vector<int>(intData::count,1);
        Vector<int> write_real_comp = Vector<int>(realData::count,1);

        pc->WritePlotFile(plotfilename, "particles", write_real_comp,
                          write_int_comp, real_comp_names, int_comp_names);

    }

}
