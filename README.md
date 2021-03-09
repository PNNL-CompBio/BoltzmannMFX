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
