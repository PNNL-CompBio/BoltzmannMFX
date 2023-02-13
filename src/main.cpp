//
//     Copyright (c) 2013 Battelle Memorial Institute
//     Licensed under modified BSD License. A copy of this license can be found
//     in the LICENSE file in the top level directory of this distribution.
//
#include <fstream>
#include <iomanip>

#include <AMReX_ParmParse.H>
#include <AMReX_Geometry.H>
#include <AMReX_VisMF.H>
#include <AMReX_buildInfo.H>

#include <AMReX_REAL.H>

#include <bmx.H>
#include <bmx_fluid_parms.H>
#ifdef NEW_CHEM
#include <bmx_chem.H>
#include <bmx_cell_interaction.H>
#endif

int  max_step   = -1;
int  regrid_int = -1;
amrex::Real stop_time  = -1.0;

std::string restart_file {""};

int check_int = -1;
int last_chk  = -1;
std::string check_file {"chk"};

std::string plot_file {"plt"};
std::string static_plt_file {"plt_ls"};

bool plotfile_on_restart = false;

int par_ascii_int = -1;
int last_par_ascii  = -1;
std::string par_ascii_file {"par"};

void set_ptr_to_bmx (bmx& bmx);

void writeBuildInfo ();

void ReadParameters ()
{
  {
     // The "amr" in the constructor means that "amr" is automatically assumed
     // to be a prefix for the variables in this particular instance, so only
     // variables prefixed with "amr" are read by this particular ParmParse
     // object. Note that check_file and check_int etc. appear in the input file
     // as amr.check_file and amr.check_int
     ParmParse pp("amr");

     pp.query("check_file", check_file);
     pp.query("check_int", check_int);

     pp.query("plot_file", plot_file);

     pp.query("plotfile_on_restart", plotfile_on_restart);

     pp.query("par_ascii_file", par_ascii_file);
     pp.query("par_ascii_int", par_ascii_int);

     pp.query("restart", restart_file);

     pp.query("regrid_int",regrid_int);

     if ( regrid_int == 0 )
       amrex::Abort("regrid_int must be > 0 or < 0");
  }

  {
     ParmParse pp("bmx");

     pp.query("stop_time", stop_time);
     pp.query("max_step", max_step);
  }

}

/**
 * Currently just guessing about these parameters.
 * @brief routine to write out information at current time step. This is
 * probably complicated by the fact that the code maybe using different time
 * increments in different parts of the code (grid and particle dynamics) and
 * that grids at different levels are also using different time increments.
 * @param[in] nstep integer value of current time step
 * @param[in] time real value of current time
 * @param[in] dt time step increment
 * @param[in] bmx simulation object
 */
void writeNow (int nstep, Real time, Real dt, bmx& bmx)
{
    int plot_test = 0;
    if (bmx::plot_per_approx > 0.0)
    {
        // Check to see if we've crossed a bmx::plot_per_approx interval by comparing
        // the number of intervals that have elapsed for both the current
        // time and the time at the beginning of this timestep.

        int num_per_old = static_cast<int>( (time-dt) / bmx::plot_per_approx );
        int num_per_new = static_cast<int>( (time   ) / bmx::plot_per_approx );

        // Before using these, however, we must test for the case where we're
        // within machine epsilon of the next interval. In that case, increment
        // the counter, because we have indeed reached the next bmx::plot_per_approx interval
        // at this point.

        const Real eps = std::numeric_limits<Real>::epsilon() * 10.0 * amrex::Math::abs(time);
        const Real next_plot_time = (num_per_old + 1) * bmx::plot_per_approx;

        if ((num_per_new == num_per_old) && amrex::Math::abs(time - next_plot_time) <= eps)
        {
            num_per_new += 1;
        }

        // Similarly, we have to account for the case where the old time is within
        // machine epsilon of the beginning of this interval, so that we don't double
        // count that time threshold -- we already plotted at that time on the last timestep.

        if ((num_per_new != num_per_old) && amrex::Math::abs((time - dt) - next_plot_time) <= eps)
            num_per_old += 1;

        if (num_per_old != num_per_new)
            plot_test = 1;

    }
    else if ( bmx::plot_per_exact  > 0 && (amrex::Math::abs(remainder(time, bmx::plot_per_exact)) < 1.e-12) )
    {
        plot_test = 1;
    }

    if ( (plot_test == 1) || ( ( bmx::plot_int > 0) && ( nstep %  bmx::plot_int == 0 ) ) )
    {
        bmx.WritePlotFile( plot_file, nstep, time );
    }


    if ( ( check_int > 0) && ( nstep %  check_int == 0 ) )
    {
        bmx.WriteCheckPointFile( check_file, nstep, dt, time );
        last_chk = nstep;
    }

    if ( ( par_ascii_int > 0) && ( nstep %  par_ascii_int == 0 ) )
    {
        bmx.WriteParticleAscii( par_ascii_file, nstep );
        last_par_ascii = nstep;
    }
}

int main (int argc, char* argv[])
{
    // check to see if it contains --describe. If so, write out information on
    // the build.
    if (argc >= 2) {
        for (auto i = 1; i < argc; i++) {
            if (std::string(argv[i]) == "--describe") {
                writeBuildInfo();
                return 0;
            }
        }
    }

    // Issue an error if AMR input file is not given
    if ( argc < 2 ) {
       std::cerr << "AMReX input file missing" << std::endl << std::endl;
       std::cerr << "Usage:  " << argv[0] << " inputs [--describe]" << std::endl;
       return -1;
    }

    // AMReX will now read the inputs file and the command line arguments, but the
    //        command line arguments are in bmx-format so it will just ignore them.
    //        The inputs file contents are available any time a ParmParser
    //        object is created. The contents can be filtered by specifying a
    //        prefix in the ParmParser constructor.
    amrex::Initialize(argc,argv,true,MPI_COMM_WORLD);
    { // This start bracket and the end bracket before Finalize are essential so
      // that the bmx object is deleted before Finalize
    BL_PROFILE_VAR("main()", pmain)
    BL_PROFILE_REGION_START("bmx::main()");

    // Write out the BMX git hash (the AMReX git hash is already written (if
    // --define specified?))
    const char* githash_bmx = buildInfoGetGitHash(1);
    amrex::Print() << "BMX git hash: " << githash_bmx<< "\n";

    // Setting format to NATIVE rather than default of NATIVE_32
    FArrayBox::setFormat(FABio::FAB_NATIVE);

    Real strt_time = ParallelDescriptor::second();

    Real time=0.0L;
    int nstep = 0;  // Current time step

    Real dt = -1.;

    // Default constructor. Note inheritance: bmx : AmrCore : AmrMesh
    //                                                             |
    //  => Geometry is constructed here: (constructs Geometry) ----+
    bmx bmx;

    // Parameters have already been read in when initialize was called.
    // These are now available to bmx
    ReadParameters();

    // Set global static pointer to bmx object. Used by fill-patch utility
    set_ptr_to_bmx(bmx);

    // Initialize internals from ParamParse database
    bmx.InitParams();

    // Initialize memory for data-array internals
    bmx.ResizeArrays();

    // Initialize derived internals
    bmx.Init(time);

    // Read in chemistry parameters. Currently hardwiring these so file
    // name is set to arbitrary string
    SPECIES::Initialize();
    amrex::Print() << "Volume threshold for cell division:  " << SPECIES::max_vol << std::endl;
    amrex::Print() << "Length threshold for segment division:  " << SPECIES::max_len << std::endl;
    amrex::Print() << "Maximum segment radius:  " << SPECIES::max_rad << std::endl;
    BMXChemistry *bmxchem = BMXChemistry::instance();
    bmxchem->setParams("NullFile");
    BMXCellInteraction *interaction = BMXCellInteraction::instance();
    interaction->setParams("NullFile");

    // Either init from scratch or from the checkpoint file
    int restart_flag = 0;
    if (restart_file.empty())
    {
        bmx.InitLevelData(time);
    }
    else
    {
        restart_flag = 1;
        bmx.Restart(restart_file, &nstep, &dt, &time);
    }

    if (FLUID::solve) {
        bmx.bmx_init_solvers();
    }

    // This checks if we want to regrid
    if (regrid_int > -1 && nstep%regrid_int == 0)
    {
        amrex::Print() << "Regridding at step " << nstep << std::endl;
        bmx.Regrid();
    }

    bmx.PostInit(dt, time, restart_flag, stop_time);

    Real end_init = ParallelDescriptor::second() - strt_time;
    ParallelDescriptor::ReduceRealMax(end_init, ParallelDescriptor::IOProcessorNumber());

    if (ParallelDescriptor::IOProcessor())
       std::cout << "Time spent in init      " << end_init << std::endl;

    int finish  = 0;

    // Initialize prev_dt here; it will be re-defined by call to evolve_fluid but
    // only if FLUID::solve = T
    Real prev_dt = dt;

    // Write checkpoint and plotfiles with the initial data
    if ( (restart_file.empty() || plotfile_on_restart) &&
         (bmx::plot_int > 0 || bmx::plot_per_exact > 0 || bmx::plot_per_approx > 0) )
    {
       bmx.WritePlotFile(plot_file, nstep, time);
    }

    // We automatically write checkpoint files with the initial data
    //    if check_int > 0
    if ( restart_file.empty() && check_int > 0 )
    {
       bmx.WriteCheckPointFile(check_file, nstep, dt, time);
       last_chk = nstep;
    }

    // We automatically write ASCII files with the particle data
    //    if par_ascii_int > 0
    if ( par_ascii_int > 0 )
    {
       bmx.WriteParticleAscii(par_ascii_file, nstep);
       last_par_ascii = nstep;
    }

    bool do_not_evolve = ( (max_step == 0) ||
                         ( (stop_time >= 0.) && (time >  stop_time) ) ||
                         ( (stop_time <= 0.) && (max_step <= 0) ) );

    if (restart_file.empty())
    {
        amrex::Print() << " " << std::endl;
        bool unused_inputs = ParmParse::QueryUnusedInputs();
        if (unused_inputs)
           amrex::Print() << "We should think about aborting here due to unused inputs" << std::endl;
    }

    { // Start profiling solve here

        BL_PROFILE("bmx_solve");
        BL_PROFILE_REGION("bmx_solve");

        bmx.ComputeAndPrintSums();

        if ( !do_not_evolve)
        {
            while (finish == 0)
            {
                Real strt_step = ParallelDescriptor::second();

                if (regrid_int > -1 && nstep%regrid_int == 0)
                {
                   amrex::Print() << "Regridding at step " << nstep << std::endl;
                   bmx.Regrid();
                }

                // This is probably the key routine to understand if we want to
                // customize code for new models The routine is defined in
                // timestepping/bmx_evolve.cpp
                bmx.Evolve(nstep, dt, prev_dt, time, stop_time);

                Real end_step = ParallelDescriptor::second() - strt_step;
                ParallelDescriptor::ReduceRealMax(end_step, ParallelDescriptor::IOProcessorNumber());
                if (ParallelDescriptor::IOProcessor())
                    std::cout << "   Time per step        " << end_step << std::endl;

                time += prev_dt;
                nstep++;

                writeNow(nstep, time, prev_dt, bmx);

                // Mechanism to terminate BMX normally.
                do_not_evolve =  ( ( (stop_time >= 0.) && (time+0.1*dt >= stop_time) ) ||
                                   ( max_step >= 0 && nstep >= max_step ) );
                if ( do_not_evolve ) finish = 1;
            }
        }
    }

    // Dump plotfile at the final time
    if ( check_int > 0 && nstep != last_chk)
        bmx.WriteCheckPointFile(check_file, nstep, dt, time);
    if ( bmx::plot_int > 0)
        bmx.WritePlotFile(plot_file, nstep, time);
    if ( par_ascii_int > 0  && nstep != last_par_ascii)
        bmx.WriteParticleAscii(par_ascii_file, nstep);

    Real end_time = ParallelDescriptor::second() - strt_time;
    ParallelDescriptor::ReduceRealMax(end_time, ParallelDescriptor::IOProcessorNumber());

    if (ParallelDescriptor::IOProcessor())
    {
        std::cout << "Time spent in main      " << end_time << std::endl;
        std::cout << "Time spent in main-init " << end_time-end_init << std::endl;
    }

    amrex::Print() << " " << std::endl;
    //bool unused_inputs = ParmParse::QueryUnusedInputs(); UNUSED VARIABLE

    BL_PROFILE_REGION_STOP("bmx::main()");
    BL_PROFILE_VAR_STOP(pmain);

    } // This end bracket and the start bracket after Initialize are essential so
      // that the bmx object is deleted before Finalize

    amrex::Finalize();
    return 0;
}
