#include <bmx.H>
#include <bmx_fluid_parms.H>
#include <bmx_fillpatch_bc.H>

#include <AMReX_FillPatchUtil.H>

namespace
{
  bmx* bmx_for_fillpatching;
}

// This interface must match the definition of the interface for
//    CpuBndryFuncFab in amrex/Src/Base/AMReX_PhysBCFunct.H
void set_ptr_to_bmx (bmx& bmx_for_fillpatching_in)
{
   bmx_for_fillpatching = &bmx_for_fillpatching_in;
}


// This interface must match the definition of the interface for
//    CpuBndryFuncFab in amrex/Src/Base/AMReX_PhysBCFunct.H
inline
void SpeciesFillBox (Box const& bx,
                     Array4<amrex::Real> const& dest,
                     const int dcomp,
                     const int numcomp,
                     GeometryData const& geom,
                     const Real time_in,
                     const BCRec* bcr,
                     const int bcomp,
                     const int orig_comp)
{
    if (dcomp != 0)
         amrex::Abort("Must have dcomp = 0 in SpeciesFillBox");
    if (numcomp != FLUID::nspecies)
         amrex::Abort("Must have numcomp = nspecies_g in SpeciesFillBox");

    const Box& domain = geom.Domain();

    // This is a bit hack-y but does get us the right level
    int lev = 0;
    for (int ilev = 0; ilev < 10; ilev++)
    {
       const Geometry& lev_geom = bmx_for_fillpatching->GetParGDB()->Geom(ilev);
       if (domain.length()[0] == (lev_geom.Domain()).length()[0])
       {
         lev = ilev;
         break;
       }
    }

    // We only do this to make it not const
    Real time = time_in;

    FArrayBox dest_fab(dest);
    Elixir eli_dest_fab = dest_fab.elixir();

   if( orig_comp == 0 )
   {
      // bmx_for_fillpatching->set_mass_fractions_g_bcs(time, lev, dest_fab, domain);
   } else if(orig_comp == 1)
      bmx_for_fillpatching->set_species_diffusivities_g_bcs(time, lev, dest_fab, domain);
   else
      amrex::Abort("Unknown component in ScalarFillBox!");

}

// Compute a new multifab by copying array from valid region and filling ghost cells
// works for single level and 2-level cases (fill fine grid ghost by interpolating from coarse)
// NOTE: icomp here refers to whether we are filling 0: fluid species
void
bmx::FillPatchSpecies (int lev,
                        Real time,
                        MultiFab& mf,
                        int icomp,
                        int ncomp,
                        const Vector<BCRec>& bcs)
{
    // Hack so that ghost cells are not undefined
    mf.setVal(covered_val);

    // icomp tells us which scalar we are fill-patching
    // But we send "0) into FillPatch since each scalar is stored in its own array

    if (lev == 0)
    {
        Vector<MultiFab*> smf;
        Vector<Real> stime;

        GetDataSpecies(0, time, smf, icomp, stime);

        CpuBndryFuncFab bfunc(SpeciesFillBox);
        PhysBCFunct<CpuBndryFuncFab> physbc(geom[lev], bcs, bfunc);
        amrex::FillPatchSingleLevel(mf, time, smf, stime, 0, 0, ncomp,
                                    geom[lev], physbc, icomp);
    }
    else
    {
        Vector<MultiFab*> cmf, fmf;
        Vector<Real> ctime, ftime;
        GetDataSpecies(lev-1, time, cmf, icomp, ctime);
        GetDataSpecies(lev  , time, fmf, icomp, ftime);

        CpuBndryFuncFab bfunc(SpeciesFillBox);
        PhysBCFunct<CpuBndryFuncFab> cphysbc(geom[lev-1],bcs,bfunc);
        PhysBCFunct<CpuBndryFuncFab> fphysbc(geom[lev  ],bcs,bfunc);

        Interpolater* mapper = &cell_cons_interp;

        amrex::FillPatchTwoLevels(mf, time, cmf, ctime, fmf, ftime,
                                  0, 0, ncomp, geom[lev-1], geom[lev],
                                  cphysbc, 0, fphysbc, 0,
                                  refRatio(lev-1), mapper, bcs, icomp);
    }
}


// Utility to copy in data from phi_old and/or phi_new into another multifab
void
bmx::GetDataSpecies (int lev,
                      Real time,
                      Vector<MultiFab*>& data,
                      int icomp,
                      Vector<Real>& datatime)
{
    data.clear();
    datatime.clear();

    const Real teps = (t_new[lev] - t_old[lev]) * 1.e-3;

    if (time > t_new[lev] - teps && time < t_new[lev] + teps)
    {
        if (icomp == 0) {
           data.push_back(m_leveldata[lev]->X_gk);
        }
        datatime.push_back(t_new[lev]);
    }
    else if (time > t_old[lev] - teps && time < t_old[lev] + teps)
    {
        if (icomp == 0) {
           data.push_back(m_leveldata[lev]->X_gko);
        }
        datatime.push_back(t_old[lev]);
    }
    else
    {
        if (icomp == 0) {
           data.push_back(m_leveldata[lev]->X_gko);
           data.push_back(m_leveldata[lev]->X_gk);
        }
        datatime.push_back(t_old[lev]);
        datatime.push_back(t_new[lev]);
    }
}

void
bmx::fillpatch_all ( Vector< MultiFab* > const& X_gk_in,
                     Real time)
{

  const int l_nspecies = FLUID::nspecies;

  for (int lev = 0; lev < nlev; lev++) {

    int state_comp, num_comp;

    if (advect_fluid_species) {
      MultiFab Sborder_X(grids[lev], dmap[lev], FLUID::nspecies, nghost, MFInfo());
      Sborder_X.setVal(0);
      state_comp = 0;
      num_comp = l_nspecies;
      FillPatchSpecies(lev, time, Sborder_X, state_comp, num_comp, bcs_X);
      MultiFab::Copy(*X_gk_in[lev], Sborder_X, 0, 0, num_comp, X_gk_in[lev]->nGrow());
    }
  }
}
