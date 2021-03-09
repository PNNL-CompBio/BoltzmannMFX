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

