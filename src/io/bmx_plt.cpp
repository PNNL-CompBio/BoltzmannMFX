//
//     Copyright (c) 2013 Battelle Memorial Institute
//     Licensed under modified BSD License. A copy of this license can be found
//     in the LICENSE file in the top level directory of this distribution.
//
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
        pp.query("plt_grad_X",plt_grad_X);
        pp.query("plt_D",     plt_D_k   );
        pp.query("plt_vf",    plt_vf    );
        pp.query("plt_np",    plt_np    );
 
        if( plt_X_k    == 1)  pltVarCount += FLUID::nchem_species;
        if( plt_grad_X == 1)  pltVarCount += AMREX_SPACEDIM * FLUID::nchem_species;
        if( plt_D_k    == 1)  pltVarCount += FLUID::nchem_species;
        if( plt_np     == 1)  pltVarCount += 1;
        if( plt_vf     == 1)  pltVarCount += 1;
    }
}

void
bmx::WriteAsciiVTK (const std::string& vtk_file, int nstep, Real time) const
{
    pc->WriteToAscii(vtk_file, nstep, time);
}

void
bmx::WritePlotFile (std::string& plot_file, int nstep, Real time )
{
    // If we've already written this plotfile, don't do it again!
    if (nstep == last_plt) return;

    // Now set last_plt to nstep ...
    last_plt = nstep;

    std::string vtk_ascii_file = "vtkfile";
    WriteAsciiVTK(vtk_ascii_file,nstep,time);

    BL_PROFILE("bmx::WritePlotFile()");

    const std::string& plotfilename = amrex::Concatenate(plot_file,nstep);

    amrex::Print() << "  Writing plotfile " << plotfilename <<  " at time " << time << std::endl;

    if (pltVarCount > 0) {

      const int ngrow = 0;

      Vector<std::string> pltFldNames;
      Vector< std::unique_ptr<MultiFab> > mf(finestLevel()+1);

      // Fluid chem_species mass fractions
      if(FLUID::solve_chem_species and plt_X_k == 1)
        for(std::string specie: FLUID::chem_species)
          pltFldNames.push_back("X_"+specie);

      // Fluid chem_species mass diffusivities
      if(FLUID::solve_chem_species and plt_D_k == 1)
        for(std::string specie: FLUID::chem_species)
          pltFldNames.push_back("D_"+specie);

      // Volume fraction
      if(plt_vf == 1)
          pltFldNames.push_back("volfrac");

      // Number of particles per grid cell
      if(plt_np == 1)
          pltFldNames.push_back("particle_count");

      // Gradient of species in fluid
      if(FLUID::solve_chem_species and plt_grad_X == 1)
      {
        for(std::string specie: FLUID::chem_species)
          pltFldNames.push_back("gx_"+specie);
        for(std::string specie: FLUID::chem_species)
          pltFldNames.push_back("gy_"+specie);
        for(std::string specie: FLUID::chem_species)
          pltFldNames.push_back("gz_"+specie);
      }


      for (int lev = 0; lev <= finestLevel(); ++lev)
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

        if(plt_vf == 1)
        {
             MultiFab::Copy(*mf[lev], *m_leveldata[lev]->vf_n, 0, lc, 1, 0);
             lc += 1;
        }

        if(plt_np == 1)
        {
             MultiFab temp_dat(grids[lev], dmap[lev], 1, 0);
             temp_dat.setVal(0);
             pc->Increment(temp_dat, lev);
             MultiFab::Copy(*mf[lev], temp_dat, 0, lc, 1, 0);
             lc += 1;
        }

        if(FLUID::solve_chem_species and plt_grad_X == 1)
        {
             MultiFab gx(grids[lev], dmap[lev], FLUID::nchem_species, 0);
             MultiFab gy(grids[lev], dmap[lev], FLUID::nchem_species, 0);
             MultiFab gz(grids[lev], dmap[lev], FLUID::nchem_species, 0);

             compute_grad_X(lev, time, gx, gy, gz);

             MultiFab::Copy(*mf[lev], gx, 0, lc, FLUID::nchem_species, 0);
             lc += FLUID::nchem_species;
             MultiFab::Copy(*mf[lev], gy, 0, lc, FLUID::nchem_species, 0);
             lc += FLUID::nchem_species;
             MultiFab::Copy(*mf[lev], gz, 0, lc, FLUID::nchem_species, 0);
             lc += FLUID::nchem_species;
        }

      } // lev

      Vector<const MultiFab*> mf2(finestLevel()+1);
      for (int lev = 0; lev <= finestLevel(); ++lev) {
        mf2[lev] = mf[lev].get();
      }

      Vector<int> istep;
      istep.resize(finestLevel()+1,nstep);
      amrex::WriteMultiLevelPlotfile(plotfilename, finestLevel()+1, mf2, pltFldNames,
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
      istep.resize(finestLevel()+1,nstep);
      amrex::WriteMultiLevelPlotfileHeaders(plotfilename, finest_level+1, mf2, names,
                                            Geom(), time, istep, refRatio());

    }

    WriteJobInfo(plotfilename);

    if ( DEM::solve )
    {
        Vector<std::string> real_comp_names;
        Vector<std::string>  int_comp_names;

        real_comp_names.push_back("radius");
        real_comp_names.push_back("cylinder_length");
        real_comp_names.push_back("theta");
        real_comp_names.push_back("phi");
        real_comp_names.push_back("cell_surface_area");
        real_comp_names.push_back("cell_volume");
        real_comp_names.push_back("velocity_x");
        real_comp_names.push_back("velocity_y");
        real_comp_names.push_back("velocity_z");
        real_comp_names.push_back("angular_velocity_x");
        real_comp_names.push_back("angular_velocity_y");
        real_comp_names.push_back("angular_velocity_z");
        real_comp_names.push_back("force_x");
        real_comp_names.push_back("force_y");
        real_comp_names.push_back("force_z");
        real_comp_names.push_back("torque_x");
        real_comp_names.push_back("torque_y");
        real_comp_names.push_back("torque_z");
        real_comp_names.push_back("grad_x");
        real_comp_names.push_back("grad_y");
        real_comp_names.push_back("grad_z");
        real_comp_names.push_back("cell_surface_area_growth_rate");
        real_comp_names.push_back("cell_volume_growth_rate");
        real_comp_names.push_back("splitting_fraction");
        real_comp_names.push_back("bond_scaling");
        real_comp_names.push_back("random_vx");
        real_comp_names.push_back("random_vy");
        real_comp_names.push_back("random_vz");
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
        int_comp_names.push_back("cell_type");
        int_comp_names.push_back("number_bonds");
        int_comp_names.push_back("unique_id_1");
        int_comp_names.push_back("unique_id_2");
        int_comp_names.push_back("unique_id_3");
        int_comp_names.push_back("unique_id_4");
        int_comp_names.push_back("home_proc_1");
        int_comp_names.push_back("home_proc_2");
        int_comp_names.push_back("home_proc_3");
        int_comp_names.push_back("home_proc_4");
        int_comp_names.push_back("site_1");
        int_comp_names.push_back("site_2");
        int_comp_names.push_back("site_3");
        int_comp_names.push_back("site_4");
        int_comp_names.push_back("position");
        int_comp_names.push_back("fusion_flag");
        int_comp_names.push_back("fusion_split_flag");
        int_comp_names.push_back("fusion_new_flag");
        int_comp_names.push_back("fusion_id");
        int_comp_names.push_back("fusion_cpu");
        int_comp_names.push_back("fuse_tip");
        int_comp_names.push_back("deleted_site_1");
        int_comp_names.push_back("deleted_site_2");
        int_comp_names.push_back("deleted_id_1");
        int_comp_names.push_back("deleted_id_2");
        int_comp_names.push_back("deleted_cpu_1");
        int_comp_names.push_back("deleted_cpu_2");
        int_comp_names.push_back("fix_site");
        int_comp_names.push_back("unique_id");
        int_comp_names.push_back("home_cpu");
        for (i=intIdx::count-1; i<MAX_CHEM_INT_VAR; i++) {
          char c[2];
          c[1] = '\0';
          c[0] = static_cast<char>(static_cast<int>('A')+(i-intIdx::count+1)%26);
          int_comp_names.push_back(c);
        }
        Vector<int> write_int_comp  = Vector<int>(MAX_CHEM_INT_VAR,1);

        pc->WritePlotFile(plotfilename, "particles", write_real_comp,
                          write_int_comp, real_comp_names, int_comp_names);

    }

}
