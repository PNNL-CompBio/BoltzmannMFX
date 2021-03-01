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

// This function initializes the attributes pltscalarVars, pltscaVarsName,
//                                          chkScalarVars, chkscaVarsName.
// If new variables need to be added to the output/checkpoint, simply add them
// here and the IO routines will automatically take care of them.
void
bmx::InitIOChkData ()
{
    if (ooo_debug) amrex::Print() << "InitIOChkData" << std::endl;

    chkChemSpeciesVarsName = {"X_gk"};
    //chkChemSpeciesVarsName = {"X_gk", "D_gk"};

    ResetIOChkData();
}


void
bmx::ResetIOChkData ()
{
  chkScalarVars.clear();
  chkScalarVars.resize(chkscaVarsName.size(), Vector< MultiFab*>(nlev));

  chkChemSpeciesVars.clear();
  chkChemSpeciesVars.resize(chkChemSpeciesVarsName.size(), Vector< MultiFab*>(nlev));

  for (int lev(0); lev < nlev; ++lev) {
    if (advect_fluid_chem_species) {
      chkChemSpeciesVars[0][lev] = m_leveldata[lev]->X_gk;
      //chkChemSpeciesVars[2][lev] = m_leveldata[lev]->D_gk;
    }
  }
}

void
bmx::WriteCheckHeader (const std::string& name,
                        int nstep,
                        Real dt,
                        Real time) const
{
   bool is_checkpoint = 1;

    if (ParallelDescriptor::IOProcessor())
    {
      std::string HeaderFileName(name + "/Header");
      VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);
      std::ofstream HeaderFile;

      HeaderFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());

      HeaderFile.open(HeaderFileName.c_str(), std::ofstream::out   |
                      std::ofstream::trunc |
                      std::ofstream::binary);

      if ( ! HeaderFile.good() )
          amrex::FileOpenFailed(HeaderFileName);

      HeaderFile.precision(17);

      if (is_checkpoint)
         HeaderFile << "Checkpoint version: 1\n";
      else
         HeaderFile << "HyperCLaw-V1.1\n";

      const int nlevels = finestLevel()+1;
      HeaderFile << nlevels << "\n";

      // Time stepping controls
      HeaderFile << nstep << "\n";
      HeaderFile << dt << "\n";
      HeaderFile << time << "\n";

      Geometry geometry;

      // Geometry
      for (int i = 0; i < BL_SPACEDIM; ++i)
            HeaderFile << geometry.ProbLo(i) << ' ';
      HeaderFile << '\n';

      for (int i = 0; i < BL_SPACEDIM; ++i)
         HeaderFile << geometry.ProbHi(i) << ' ';
      HeaderFile << '\n';

      // BoxArray
      for (int lev = 0; lev < nlevels; ++lev)
      {
          boxArray(lev).writeOn(HeaderFile);
          HeaderFile << '\n';
      }
    }
}

void
bmx::WriteCheckPointFile (std::string& check_file,
                           int nstep,
                           Real dt,
                           Real time)
{
    BL_PROFILE("bmx::WriteCheckPointFile()");

    const std::string& checkpointname = amrex::Concatenate( check_file, nstep );

    if (ParallelDescriptor::IOProcessor()) {
      std::cout << "\n\t Writing checkpoint " << checkpointname << std::endl;
    }

    const int nlevels = finestLevel()+1;
    amrex::PreBuildDirectorHierarchy(checkpointname, level_prefix, nlevels, true);

    WriteCheckHeader(checkpointname, nstep, dt, time);

    WriteJobInfo(checkpointname);
    if (FLUID::solve)
    {
       ResetIOChkData();

       for (int lev = 0; lev < nlevels; ++lev) {

          if (advect_fluid_chem_species) {
             // Write chem_species variables
             for (int i = 0; i < chkChemSpeciesVars.size(); i++) {
                VisMF::Write( *(chkChemSpeciesVars[i][lev]),
                  amrex::MultiFabFileFullPrefix(lev, checkpointname,
                        level_prefix, chkChemSpeciesVarsName[i]));
             }
          }
       }
    }

    if ( DEM::solve )
    {
       pc->Checkpoint(checkpointname, "particles");
    }
}
