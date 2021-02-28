#include <bmx_leveldata.H>
#include <bmx_fluid_parms.H>
#include <bmx_species_parms.H>

using namespace amrex;

LevelData::LevelData (BoxArray const& ba,
                      DistributionMapping const& dmap,
                      const int nghost)
  : X_gk(nullptr)
  , X_gko(nullptr)
  , X_rhs(nullptr)
  , D_gk(nullptr)
{
    amrex::Print() << "MAKING ARRAYS " << FLUID::nspecies << std::endl;
    X_gk  = new MultiFab(ba, dmap, FLUID::nspecies, nghost, MFInfo());
    X_gko = new MultiFab(ba, dmap, FLUID::nspecies, nghost, MFInfo());
    X_rhs = new MultiFab(ba, dmap, FLUID::nspecies, nghost, MFInfo());
    D_gk  = new MultiFab(ba, dmap, FLUID::nspecies, nghost, MFInfo());
}

LevelData::~LevelData ()
{
    delete X_gk;
    delete X_gko;
    delete X_rhs;
    delete D_gk;
}
