//
//     Copyright (c) 2013 Battelle Memorial Institute
//     Licensed under modified BSD License. A copy of this license can be found
//     in the LICENSE file in the top level directory of this distribution.
//
#ifndef _BMX_LEVELDATA_H_
#define _BMX_LEVELDATA_H_

#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_FabFactory.H>
#include <AMReX_MultiFab.H>

class LevelData
{
  public:
    // Default Constructor
    LevelData () = default;

    // Constructor
    LevelData (amrex::BoxArray const& ba,
               amrex::DistributionMapping const& dm,
               const int nghost);

    // Destructor
    ~LevelData ();

    // Gas chem_species (new and old)
    amrex::MultiFab* X_k;
    amrex::MultiFab* X_ko;

    // Gas chem_species source/sink from particles
    amrex::MultiFab* X_rhs;

    // ChemSpecies diffusion coefficients
    amrex::MultiFab* D_k;

    // Volume fraction of agar/biofilm in grid cell (= 1. - volume fraction of biological cell)
    amrex::MultiFab* vf_o;
    amrex::MultiFab* vf_n;
};

#endif
