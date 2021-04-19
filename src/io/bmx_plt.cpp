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

        pp.query("plt_X",     plt_X_k   );
        pp.query("plt_D",     plt_D_k   );
 
        if( plt_X_k == 1)  pltVarCount += FLUID::nchem_species;
        if( plt_D_k == 1)  pltVarCount += FLUID::nchem_species;
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

      // Fluid chem_species mass fractions
      if(FLUID::solve_chem_species and plt_X_k == 1)
        for(std::string specie: FLUID::chem_species)
          pltFldNames.push_back("X_"+specie);

      // Fluid chem_species mass diffusivities
      if(FLUID::solve_chem_species and plt_D_k == 1)
        for(std::string specie: FLUID::chem_species)
          pltFldNames.push_back("D_"+specie);

      for (int lev = 0; lev < nlev; ++lev)
      {
        // Multifab to hold all the variables -- there can be only one!!!!
        const int ncomp = pltVarCount;
        mf[lev].reset(new MultiFab(grids[lev], dmap[lev], ncomp, ngrow,  MFInfo()));

        int lc=0;

        if(FLUID::solve_chem_species and plt_X_k == 1 ) {
          for(int n(0); n < FLUID::nchem_species; n++) {
            MultiFab::Copy(*mf[lev], *m_leveldata[lev]->X_k, n, lc+n, 1, 0);
          }
          lc += FLUID::nchem_species;
        }

        // ChemSpecies mass fraction
        if(FLUID::solve_chem_species and plt_D_k == 1 ) {
          for(int n(0); n < FLUID::nchem_species; n++) {
            MultiFab::Copy(*mf[lev], *m_leveldata[lev]->D_k, n, lc+n, 1, 0);
          }
          lc += FLUID::nchem_species;
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

#ifdef NEW_CHEM
        real_comp_names.push_back("a axis");
        real_comp_names.push_back("b axis");
        real_comp_names.push_back("c axis");
        real_comp_names.push_back("psi");
        real_comp_names.push_back("theta");
        real_comp_names.push_back("phi");
        real_comp_names.push_back("cell surface area");
        real_comp_names.push_back("cell volume");
        real_comp_names.push_back("velocity x");
        real_comp_names.push_back("velocity y");
        real_comp_names.push_back("velocity z");
        real_comp_names.push_back("angular velocity x");
        real_comp_names.push_back("angular velocity y");
        real_comp_names.push_back("angular velocity z");
        real_comp_names.push_back("force x");
        real_comp_names.push_back("force y");
        real_comp_names.push_back("force z");
        real_comp_names.push_back("torque x");
        real_comp_names.push_back("torque y");
        real_comp_names.push_back("torque z");
        real_comp_names.push_back("cell surface area growth rate");
        real_comp_names.push_back("cell volume growth rate");
        Vector<int> write_real_comp = Vector<int>(MAX_CHEM_REAL_VAR,1);
        int i;
        for (i=realIdx::count-1; i<MAX_CHEM_REAL_VAR; i++) {
          char c[2];
          c[1] = '\0';
          c[0] = static_cast<char>(static_cast<int>('A')+(i-realIdx::count+1)%26);
          real_comp_names.push_back(c);
        }
        int_comp_names.push_back("num_reals");
        int_comp_names.push_back("num_ints");
        int_comp_names.push_back("total_reals");
        int_comp_names.push_back("total_ints");
        int_comp_names.push_back("first_conc_inc");
        for (i=intIdx::count-1; i<MAX_CHEM_INT_VAR; i++) {
          char c[2];
          c[1] = '\0';
          c[0] = static_cast<char>(static_cast<int>('A')+(i-intIdx::count+1)%26);
          int_comp_names.push_back(c);
        }
        Vector<int> write_int_comp  = Vector<int>(MAX_CHEM_INT_VAR,1);
#else
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
#endif

        pc->WritePlotFile(plotfilename, "particles", write_real_comp,
                          write_int_comp, real_comp_names, int_comp_names);

    }

}
