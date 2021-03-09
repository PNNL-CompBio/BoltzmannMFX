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
};

#endif
