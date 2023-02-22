//
//     Copyright (c) 2013 Battelle Memorial Institute
//     Licensed under modified BSD License. A copy of this license can be found
//     in the LICENSE file in the top level directory of this distribution.
//
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
void ChemSpeciesFillBox (Box const& /*bx*/,
                         Array4<amrex::Real> const& dest,
                         const int dcomp,
                         const int numcomp,
                         GeometryData const& geom_data,
                         const Real time_in,
                         const BCRec* /*bcr*/,
                         const int /*bcomp*/,
                         const int /*orig_comp*/)
{
    if (dcomp != 0)
         amrex::Abort("Must have dcomp = 0 in ChemSpeciesFillBox");
    if (numcomp != FLUID::nchem_species)
         amrex::Abort("Must have numcomp = nchem_species in ChemSpeciesFillBox");

    const Box& domain = geom_data.Domain();

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

    bmx_for_fillpatching->set_neumann_bcs(time, lev, dest_fab, geom_data);

}

// This interface must match the definition of the interface for
//    CpuBndryFuncFab in amrex/Src/Base/AMReX_PhysBCFunct.H
inline
void ChemCoeffsFillBox (Box const& /*bx*/,
                        Array4<amrex::Real> const& dest,
                        const int dcomp,
                        const int numcomp,
                        GeometryData const& geom_data,
                        const Real time_in,
                        const BCRec* /*bcr*/,
                        const int /*bcomp*/,
                        const int /*orig_comp*/)
{
    if (dcomp != 0)
         amrex::Abort("Must have dcomp = 0 in ChemCoeffsFillBox");
    if (numcomp != FLUID::nchem_species)
         amrex::Abort("Must have numcomp = nchem_species in ChemCoeffsFillBox");

    const Box& domain = geom_data.Domain();

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

    bmx_for_fillpatching->set_neumann_bcs(time, lev, dest_fab, geom_data);
}

// This interface must match the definition of the interface for
//    CpuBndryFuncFab in amrex/Src/Base/AMReX_PhysBCFunct.H
inline
void VolFracFillBox (Box const& /*bx*/,
                     Array4<amrex::Real> const& dest,
                     const int dcomp,
                     const int numcomp,
                     GeometryData const& geom_data,
                     const Real time_in,
                     const BCRec* /*bcr*/,
                     const int /*bcomp*/,
                     const int /*orig_comp*/)
{
    if (dcomp != 0)
         amrex::Abort("Must have dcomp = 0 in VolFracFillBox");
    if (numcomp != 1)
         amrex::Abort("Must have numcomp = 1 in VolFracFillBox");

    const Box& domain = geom_data.Domain();

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

    bmx_for_fillpatching->set_neumann_bcs(time, lev, dest_fab, geom_data);
}

// Compute a new multifab by copying array from valid region and filling ghost cells
// works for single level and 2-level cases (fill fine grid ghost by interpolating from coarse)
void
bmx::FillPatchChemSpecies (int lev,
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

        GetDataChemSpecies(0, time, smf, icomp, stime);

        CpuBndryFuncFab bfunc(ChemSpeciesFillBox);
        PhysBCFunct<CpuBndryFuncFab> physbc(geom[lev], bcs, bfunc);
        amrex::FillPatchSingleLevel(mf, time, smf, stime, 0, 0, ncomp,
                                    geom[lev], physbc, icomp);
    }
    else
    {
        Vector<MultiFab*> cmf, fmf;
        Vector<Real> ctime, ftime;
        GetDataChemSpecies(lev-1, time, cmf, icomp, ctime);
        GetDataChemSpecies(lev  , time, fmf, icomp, ftime);

        CpuBndryFuncFab bfunc(ChemSpeciesFillBox);
        PhysBCFunct<CpuBndryFuncFab> cphysbc(geom[lev-1],bcs,bfunc);
        PhysBCFunct<CpuBndryFuncFab> fphysbc(geom[lev  ],bcs,bfunc);

        Interpolater* mapper = &cell_cons_interp;

        amrex::FillPatchTwoLevels(mf, time, cmf, ctime, fmf, ftime,
                                  0, 0, ncomp, geom[lev-1], geom[lev],
                                  cphysbc, 0, fphysbc, 0,
                                  refRatio(lev-1), mapper, bcs, icomp);
    }
}

// Compute a new multifab by copying array from valid region and filling ghost cells
// works for single level and 2-level cases (fill fine grid ghost by interpolating from coarse)
void
bmx::FillPatchChemCoeffs (int lev,
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

        GetDataChemCoeffs(0, time, smf, icomp, stime);

        CpuBndryFuncFab bfunc(ChemCoeffsFillBox);
        PhysBCFunct<CpuBndryFuncFab> physbc(geom[lev], bcs, bfunc);
        amrex::FillPatchSingleLevel(mf, time, smf, stime, 0, 0, ncomp,
                                    geom[lev], physbc, icomp);
    }
    else
    {
        Vector<MultiFab*> cmf, fmf;
        Vector<Real> ctime, ftime;
        GetDataChemCoeffs(lev-1, time, cmf, icomp, ctime);
        GetDataChemCoeffs(lev  , time, fmf, icomp, ftime);

        CpuBndryFuncFab bfunc(ChemCoeffsFillBox);
        PhysBCFunct<CpuBndryFuncFab> cphysbc(geom[lev-1],bcs,bfunc);
        PhysBCFunct<CpuBndryFuncFab> fphysbc(geom[lev  ],bcs,bfunc);

        Interpolater* mapper = &cell_cons_interp;

        amrex::FillPatchTwoLevels(mf, time, cmf, ctime, fmf, ftime,
                                  0, 0, ncomp, geom[lev-1], geom[lev],
                                  cphysbc, 0, fphysbc, 0,
                                  refRatio(lev-1), mapper, bcs, icomp);
    }
}

// Compute a new multifab by copying array from valid region and filling ghost cells
// works for single level and 2-level cases (fill fine grid ghost by interpolating from coarse)
void
bmx::FillPatchVolFrac (int lev,
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

        GetDataVolFrac(0, time, smf, icomp, stime);

        CpuBndryFuncFab bfunc(VolFracFillBox);
        PhysBCFunct<CpuBndryFuncFab> physbc(geom[lev], bcs, bfunc);
        amrex::FillPatchSingleLevel(mf, time, smf, stime, 0, 0, ncomp,
                                    geom[lev], physbc, icomp);
    }
    else
    {
        Vector<MultiFab*> cmf, fmf;
        Vector<Real> ctime, ftime;
        GetDataVolFrac(lev-1, time, cmf, icomp, ctime);
        GetDataVolFrac(lev  , time, fmf, icomp, ftime);

        CpuBndryFuncFab bfunc(VolFracFillBox);
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
bmx::GetDataChemSpecies (int lev,
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
           data.push_back(m_leveldata[lev]->X_k);
        }
        datatime.push_back(t_new[lev]);
    }
    else if (time > t_old[lev] - teps && time < t_old[lev] + teps)
    {
        if (icomp == 0) {
           data.push_back(m_leveldata[lev]->X_ko);
        }
        datatime.push_back(t_old[lev]);
    }
    else
    {
        if (icomp == 0) {
           data.push_back(m_leveldata[lev]->X_ko);
           data.push_back(m_leveldata[lev]->X_k);
        }
        datatime.push_back(t_old[lev]);
        datatime.push_back(t_new[lev]);
    }
}

// Utility to copy in data 
void
bmx::GetDataChemCoeffs (int lev,
                        Real /*time*/,
                        Vector<MultiFab*>& data,
                        int icomp,
                        Vector<Real>& datatime)
{
    data.clear();
    datatime.clear();

    if (icomp == 0) {
       data.push_back(m_leveldata[lev]->D_k);
    }
    datatime.push_back(t_new[lev]);
}


// Utility to copy in data
void
bmx::GetDataVolFrac (int lev,
                     Real /*time*/,
                     Vector<MultiFab*>& data,
                     int icomp,
                     Vector<Real>& datatime)
{
    data.clear();
    datatime.clear();

    if (icomp == 0) {
       data.push_back(m_leveldata[lev]->vf_n);
    }
    datatime.push_back(t_new[lev]);
}

void
bmx::fillpatch_Xk ( Vector< MultiFab* > const& X_k_in,
                    Real time)
{

  const int l_nchem_species = FLUID::nchem_species;

  for (int lev = 0; lev <= finestLevel(); lev++) {

    int state_comp, num_comp;

    if (advect_fluid_chem_species) {
      MultiFab Sborder_Xk(grids[lev], dmap[lev], FLUID::nchem_species, nghost, MFInfo());
      Sborder_Xk.setVal(0);
      state_comp = 0;
      num_comp = l_nchem_species;
      FillPatchChemSpecies(lev, time, Sborder_Xk, state_comp, num_comp, bcs_X);
      MultiFab::Copy(*X_k_in[lev], Sborder_Xk, 0, 0, num_comp, X_k_in[lev]->nGrow());
    }
  }
}

void
bmx::fillpatch_Dk ( Vector< MultiFab* > const& D_k_in,
                    Real time)
{

  const int l_nchem_species = FLUID::nchem_species;

  for (int lev = 0; lev <= finestLevel(); lev++) {

    int state_comp, num_comp;

    if (advect_fluid_chem_species) {
        MultiFab Sborder_Dk(grids[lev], dmap[lev], FLUID::nchem_species, nghost, MFInfo());
        Sborder_Dk.setVal(0);
        state_comp = 0;
        num_comp = l_nchem_species;
        FillPatchChemCoeffs(lev, time, Sborder_Dk, state_comp, num_comp, bcs_D);
        MultiFab::Copy(*D_k_in[lev], Sborder_Dk, 0, 0, num_comp, D_k_in[lev]->nGrow());
    }
  }
}

void
bmx::fillpatch_vf ( Vector< MultiFab* > const& vf_in,
                    Real time)
{
  for (int lev = 0; lev <= finestLevel(); lev++) {

    int state_comp, num_comp;

    MultiFab Sborder_vf(grids[lev], dmap[lev], 1, nghost, MFInfo());
    Sborder_vf.setVal(0);
    state_comp = 0;
    num_comp = 1;
    FillPatchVolFrac(lev, time, Sborder_vf, state_comp, num_comp, bcs_f);
    MultiFab::Copy(*vf_in[lev], Sborder_vf, 0, 0, num_comp, vf_in[lev]->nGrow());
  }
}
