Some notes on updating the version of AMReX used by the CMake build system. The
variable `BMX_HOME` refers to the top-level directory in a clone of the
BoltzmannMFX repository. The following actions will update the version of AMReX
to the latest Github version

```
cd $BMX_HOME
```
```
git submodule foreach git pull git@github.com:AMReX-codes/amrex.git
```
```
cd subprojects
```
```
git add amrex
```
```
git commit
```
```
git push
```
