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

#include <bmx.H>
#include <bmx_fluid_parms.H>
#include <bmx_dem_parms.H>

namespace
{
    const std::string level_prefix {"Level_"};
}

void
bmx::GotoNextLine (std::istream& is)
{
    constexpr std::streamsize bl_ignore_max { 100000 };
    is.ignore(bl_ignore_max, '\n');
}

void bmx::WriteJobInfo (const std::string& dir) const
{
    if (ParallelDescriptor::IOProcessor())
    {
       // job_info file with details about the run
       std::ofstream jobInfoFile;
       std::string FullPathJobInfoFile = dir;
       std::string PrettyLine = "========================================="
                                "======================================\n";

       FullPathJobInfoFile += "/job_info";
       jobInfoFile.open(FullPathJobInfoFile.c_str(), std::ios::out);

       // job information
       jobInfoFile << PrettyLine;
       jobInfoFile << " BMX Job Information\n";
       jobInfoFile << PrettyLine;

       jobInfoFile << "number of MPI processes: " << ParallelDescriptor::NProcs() << "\n";
#ifdef _OPENMP
       jobInfoFile << "number of threads:       " << omp_get_max_threads() << "\n";
#endif

       jobInfoFile << "\n\n";

        // build information
       jobInfoFile << PrettyLine;
       jobInfoFile << " Build Information\n";
       jobInfoFile << PrettyLine;

       jobInfoFile << "build date:    " << buildInfoGetBuildDate() << "\n";
       jobInfoFile << "build machine: " << buildInfoGetBuildMachine() << "\n";
       jobInfoFile << "build dir:     " << buildInfoGetBuildDir() << "\n";
       jobInfoFile << "AMReX dir:     " << buildInfoGetAMReXDir() << "\n";

       jobInfoFile << "\n";

       jobInfoFile << "COMP:          " << buildInfoGetComp() << "\n";
       jobInfoFile << "COMP version:  " << buildInfoGetCompVersion() << "\n";

       jobInfoFile << "\n";

       jobInfoFile << "C++ compiler:  " << buildInfoGetCXXName() << "\n";
       jobInfoFile << "C++ flags:     " << buildInfoGetCXXFlags() << "\n";

       jobInfoFile << "\n";

       jobInfoFile << "Link flags:    " << buildInfoGetLinkFlags() << "\n";
       jobInfoFile << "Libraries:     " << buildInfoGetLibraries() << "\n";

       const char* githash1 = buildInfoGetGitHash(1);
       const char* githash2 = buildInfoGetGitHash(2);
       if (strlen(githash1) > 0)
       {
           jobInfoFile << "BMX git hash: " << githash1 << "\n";
       }
       if (strlen(githash2) > 0)
       {
           jobInfoFile << "AMReX git hash: " << githash2 << "\n";
       }

       jobInfoFile << "\n\n";

       // grid information
        jobInfoFile << PrettyLine;
        jobInfoFile << " Grid Information\n";
        jobInfoFile << PrettyLine;

        for (int i = 0; i <= finest_level; i++)
        {
            jobInfoFile << " level: " << i << "\n";
            jobInfoFile << "   number of boxes = " << grids[i].size() << "\n";
            jobInfoFile << "   maximum zones   = ";
            for (int n = 0; n < BL_SPACEDIM; n++)
            {
                jobInfoFile << geom[i].Domain().length(n) << " ";
            }
            jobInfoFile << "\n\n";
        }

        jobInfoFile << "\n\n";

        // runtime parameters
        jobInfoFile << PrettyLine;
        jobInfoFile << " Inputs File Parameters\n";
        jobInfoFile << PrettyLine;

        // ParmParse::dumpTable(jobInfoFile, true);

        jobInfoFile.close();
  }
}

void
bmx::WriteParticleAscii ( std::string& par_ascii_file, int nstep ) const
{
    BL_PROFILE("bmx::WriteParticleASCII()");

    if(DEM::solve) {

        const std::string& par_filename = amrex::Concatenate(par_ascii_file,nstep);
        pc->WriteAsciiFile(par_filename);
    }
}

