# Biological Modeling and interface eXchange (BMX)

BMX is an adaptive mesh refinement based code designed to simulate microbial and
fungal systems at the level multiple cells. The code can simulate individual
bacteria or yeast, represented as spheres and fungal networks, represented as
connected rods. The code can evaluate reactions inside the individual cells or
hyphae, transport between the cell interior and exterior medium, diffusion of
material in the exterior medium and, in the case of hyphae, direct transfer
between different fungal segments. The code also supports mechanical
interactions between different cells and fungal elements as well as low Reynolds
number motion of the cells in the suppporting fluid.

The current BMX code has implemented a static refinement scheme that is
specified by the user in the input file and remains fixed for the duration of
the calculation. The underlying AMReX framework, which implements most of the
grid-based calculation and particle motion, has extensive support for adaptive
AMR but this has not been needed to date. BMX requires two files if
starting a calculation from scratch, an input file with a complete listing of
parameters used in the simulation and a file that provides an initial
configuration for the calculation. These two files are described in more detail
below.

## Input File

A typical input file for a BMX calculation is given below

```
# Example input file

# Random number seed
bmx.seed = 43343

#_______________________________________________________________________
# Solver settings

amrex.fpe_trap_invalid = 1

#! Fluid solver
#!-----------------------------------------------------------//

# Simulation time step
bmx.fixed_dt = 0.25

# Maximum number of time steps
bmx.max_step = 400000

dem.solve = bacteria
# Initial configuration file format
bmx.particle_init_type = "AsciiFile"

# Flag to control level of verbosity. Useful for debugging
bmx.verbose = 0

# This deposition scheme flags choose whether material from a cell is
# deposited only into the grid cell containing the biological cell or
# if it is deposited using a weighting scheme into both the grid cell
# containing the biological cell and neighboring grid cells. At this
# point in time the one_to_one option appears to be the most stable.
#_______________________________________________________________________
# Deposition scheme
bmx.cnc_deposition_scheme = one_to_one # trilinear # one_to_one
bmx.vf_deposition_scheme = one_to_one # trilinear # one_to_one

#_______________________________________________________________________
# Geometry / grids / tiles

# Maximum level in grid refinement hierarchy
amr.max_level = 2

geometry.coord_sys   = 0                       # 0: Cartesian
geometry.is_periodic = 1       1       0       # Is periodic in each direction?
geometry.prob_lo     = 0.      0.      0.      # lo corner of physical domain
geometry.prob_hi     = 0.1280  0.1280  0.0640  # hi corner of physical domain

# Number of grid cells in each direction at the coarsest level
amr.n_cell =  16  16  8

# Define region at highest level of refinement
bmx.tag_region    = true
bmx.tag_region_lo = 0.     0.     0.0400
bmx.tag_region_hi = 0.1280 0.1280 0.0520

# Width of buffer region around the tagged cells
amr.n_error_buf = 2

# Just require grid to be coarsenable by 2
amr.blocking_factor = 2

#! Grids
#!-----------------------------------------------------------//
# Maximum allowable size of each fluid subdomain in the problem domain;

# Fluid
amr.max_grid_size_x = 64
amr.max_grid_size_y = 64
amr.max_grid_size_z = 64

# Particles (not with KDTree)
#particles.max_grid_size_x =   32
#particles.max_grid_size_y =   32
#particles.max_grid_size_z =   32

#! Tiles
#!-----------------------------------------------------------//

# Fluid: Maximum tile size within each grid
fabarray.mfiter_tile_size = 1024 1024 1024

# Initial configuration file for particles
particles.input_file = "fungi_init_cfg.dat"
# Particles: Maximum particle tile size
particles.tile_size = 1024 1024 1024

#_______________________________________________________________________
# Particle load balancing
#bmx.load_balance_type = "KnapSack"

amr.dual_grid          = 0
amr.regrid_int         = -1

#! KnapSack settings
#!-----------------------------------------------------------//
#default is "RunTimeCosts"; options include RunTimeCosts or NumParticles
#this option is only relevant if load_balance_type = KnapSack

#bmx.knapsack_weight_type = "NumParticles"
#bmx.knapsack_weight_type = "RunTimeCosts"

#_______________________________________________________________________
# IO / Checkpointing: the variable with _int at the end are the number
#                     of time steps between writing out the files and
#                     _file variables are the root name of the file
#                     or directory that are being exported. If the
#                     time step interval is -1, then the files are not
#                     written.

amr.par_ascii_int  = -1
amr.par_ascii_file ="vis"

amr.plot_int       = 400
amr.plot_file      ="plt"

amr.check_int      = -1
amr.check_file     ="chk"

#! Restart from checkpoint
#! If this parameter is specified, then the simulation starts from a
#! checkpoint file instead of the configuration specified in
#! particles.input_file.
#!-----------------------------------------------------------//
#amr.restart   ="chk00000"

#_______________________________________________________________________
# Fluid model settings
#
fluid.solve = fluid
fluid.chem_species      = A B C
fluid.chem_species_diff = 6.0e-10 0.0 6.0e-10
fluid.init_conc_species = 2.0e-5 0.0 2.0e-6
fluid.surface_location = 0.04800
#_______________________________________________________________________
# Reaction model parameters
#
chem_species.max_vol = 6.872e-10 # 9.5e-13 # 1.0e-12
chem_species.seg_split_length = 30.0e-4 # 9.5e-13 # 1.0e-12
chem_species.max_seg_length = 35.0e-4 # 9.5e-13 # 1.0e-12
chem_species.max_seg_radius = 2.5e-4 # 9.5e-13 # 1.0e-12
chem_species.k1 = 2.0e-5
chem_species.kr1 = 2.0e-5
chem_species.k2 = 4.0
chem_species.kr2 = 0.0006
chem_species.k3 = 2.0e-5
chem_species.kr3 = 2.0e-5
chem_species.kg = 40.0 # 400.0
chem_species.kv = 400000.0

chem_species.mass_transfer_A = 2.0e-5
chem_species.mass_transfer_B = 2.0e-5
chem_species.mass_transfer_C = 2.0e-5

chem_species.branching_probability = 0.0000
chem_species.splitting_probability = 0.2
chem_species.fusion_probability = 0.01
chem_species.max_fusion_separation = 6.0e-4

bmx.substeps = 4

#_______________________________________________________________________
# Cell force parameters
#
cell_force.stiffness = 1.0e6
cell_force.wall_stiffness = 1.0e5

cell_force.boundary_width = 1.0e-4
cell_force.wall_boundary_width = 1.0e-04

cell_force.neighbor_width = 1.0e-2 #set range of neighbor list

cell_force.bond_strength = 1.0
cell_force.bond_cutoff = 5.0e-4

cell_force.viscous_drag = 200.0

cell_force.gravity = 0.01

#cell_force.repulsion_stiffness = 1.0e7 # 5.0e6
#cell_force.adhesion_stiffness = 1.0e6
#cell_force.wall_repulsion_stiffness = 1.0e7 # 5.0e6
#cell_force.wall_adhesion_stiffness = 1.0e6

#Diffusion_type: 0 = explicit, 1 = C-N, 2 = implicit
bmx.diffusion_type = 2

#Verbosity of diffusion solver
diffusion.verbose = 2

#What goes into the plotfile?
amr.plt_X = 1
amr.plt_D = 1
amr.plt_grad_X = 1

#_______________________________________________________________________
# DEM model settings
#

```

# Notes

The material below represents notes on developing, building and running BMX.
They are still fairly crude but may be of some use to others interested in
working with the code.

## AMReX submodule

To pull in AMReX as a submodule, run the following commands in the BoltzmannMFX
directory

- remove the existing subprojects/amrex directory (if there is one)
```
rm -rf subprojects/amrex
```

- initialize and update the subprojects directory
```
git submodule init
```
```
git submodule update
```

The following single command will also probably work
```
git submodule update --init
```
## Configure and Run BMX

To configure BoltzmannMFX using CMake, create a build directory at the top
level, cd into it and run

```
cmake ..
```

(you can add the option  <code>-DCMAKE_BUILD_TYPE=Debug</code> for additional
debugging information.

After configuration is complete, type <code>make</code> to build the executable.

The test case can be run by cd'ing into the <code>exec</code> directory from the
top level directory and typing

```
../build/bmx inputs
```

Depending on where the build directory is located, the location of the
<code>bmx</code> executable may be different.

## Notes for FALLACY project

The following configuration script works on constance and can be used to build
the BMX code for running on CPUs only.

```
rm -rf CMake*

cmake -DCMAKE_BUILD_TYPE=Release \
      -DAMReX_ASSERTIONS=YES \
      -DCMAKE_VERBOSE_MAKEFILE=ON \
      -DAMReX_TINY_PROFILE=ON \
      ..
```
This script works with the GNU compilers and Open MPI.

For running with GPUs, the following script works on newell.

```
rm -rf CMake*

cmake -DCMAKE_BUILD_TYPE=Release \
      -DAMReX_ASSERTIONS=YES \
      -DBMX_GPU_BACKEND=CUDA \
      -DAMReX_GPU_BACKEND=CUDA \
      -DCMAKE_CUDA_FLAGS=--std=c++14 \
      -DCMAKE_VERBOSE_MAKEFILE=ON \
      -DAMReX_CUDA_ARCH=Volta \
      -DNVCC_ARCH_FLAGS="-gencode=arch=compute_70,code=sm_70" \
      -DAMReX_TINY_PROFILE=ON \
      ..
```
This script works with the following environment (on newell)
```
module purge
module load gcc/8.3.0
module load cmake/3.19.6
module load openmpi-gpu/4.1.0
module load cuda/10.2

setenv OMPI_MCA_pml "ucx"
setenv OMPI_MCA_btl "^vader,tcp,openib,uct"
setenv UCX_NET_DEVICES mlx5_1:1,mlx5_3:1
setenv OMP_NUM_THREADS 1
```

Two test problems are currently exported to the build directory. They are located
in BUILD/exec/nlev_large_test and BUILD/exec/nlev_real_test. The nlev_real_test
problem is the smaller of the two. To run these problems, simple cd into
these directories and launch with the desired number of processors (and GPUs, if
applicable). The input deck for both tests is input_real_nlev, so typing
```
mpirun -n 4 ../../bmx input_real_nlev
```
should get the code to run. Both tests are set to run for 100000 steps, which is
quite long. You can shorten the tests to run for a few seconds by setting the
number of steps, which is the `bmx.max_step` parameter in the input file, to
something like 10.

The GPU code can be launched with a script such as
```
#!/bin/csh
#SBATCH -t 02:30:00
#SBATCH -A dmc_biology
#SBATCH -p newell_test
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --gres=gpu:4
#SBATCH -o ./test.out
#SBATCH -e ./test.err

mpirun -n 4 ../../bmx input_real_nlev
```
