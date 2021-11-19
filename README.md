Just collecting some notes at the moment.

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
## Configure and Run BoltzmannMFX

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
