//
//     Copyright (c) 2013 Battelle Memorial Institute
//     Licensed under modified BSD License. A copy of this license can be found
//     in the LICENSE file in the top level directory of this distribution.
//
#include <bmx_leveldata.H>
#include <bmx_fluid_parms.H>
#include <bmx_chem_species_parms.H>

using namespace amrex;

LevelData::LevelData (BoxArray const& ba,
                      DistributionMapping const& dmap,
                      const int nghost)
  : X_k(nullptr)
  , X_ko(nullptr)
  , X_rhs(nullptr)
  , D_k(nullptr)
  , vf_o(nullptr)
  , vf_n(nullptr)
{
    X_k   = new MultiFab(ba, dmap, FLUID::nchem_species , nghost, MFInfo());
    X_ko  = new MultiFab(ba, dmap, FLUID::nchem_species , nghost, MFInfo());
    X_rhs = new MultiFab(ba, dmap, FLUID::nchem_species, nghost, MFInfo());
    D_k   = new MultiFab(ba, dmap, FLUID::nchem_species , nghost, MFInfo());
    vf_o  = new MultiFab(ba, dmap,                    1 , nghost, MFInfo());
    vf_n  = new MultiFab(ba, dmap,                    1 , nghost, MFInfo());
}

LevelData::~LevelData ()
{
    delete X_k;
    delete X_ko;
    delete X_rhs;
    delete D_k;
    delete vf_o;
    delete vf_n;
}
