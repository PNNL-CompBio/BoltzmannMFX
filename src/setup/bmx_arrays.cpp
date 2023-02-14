//
//     Copyright (c) 2013 Battelle Memorial Institute
//     Licensed under modified BSD License. A copy of this license can be found
//     in the LICENSE file in the top level directory of this distribution.
//
#include <bmx.H>
#include <bmx_fluid_parms.H>
#include <bmx_chem_species_parms.H>

void
bmx::ResizeArrays ()
{
    int nlevs_max = maxLevel() + 1;

    m_leveldata.resize(nlevs_max);
    for (int lev(0); lev < nlevs_max; ++lev)
      m_leveldata[lev].reset(new LevelData());

    bcoeff.resize(nlevs_max);
}

void
bmx::AllocateArrays (int lev)
{
    if (ooo_debug) amrex::Print() << "AllocateArrays" << std::endl;

    // ********************************************************************************
    // Cell- or node-based arrays
    // ********************************************************************************

    m_leveldata[lev].reset(new LevelData(grids[lev], dmap[lev], nghost));

    // ********************************************************************************
    // X-face-based arrays
    // ********************************************************************************

    // ****************************************************************

    BoxArray ba = grids[lev];
    // Create a BoxArray on x-faces.
    // x-face-based coefficient for MAC and diffusive solves
    if (bcoeff[lev][0] != nullptr)
      delete bcoeff[lev][0];

    bcoeff[lev][0] = new MultiFab(BoxArray(ba).surroundingNodes(0), dmap[lev], 1,
                                  nghost, MFInfo());
    bcoeff[lev][0]->setVal(0.);

    // Create a BoxArray on y-faces.
    // y-face-based coefficient for MAC and diffusive solves
    if (bcoeff[lev][1] != nullptr)
      delete bcoeff[lev][1];

    bcoeff[lev][1] = new MultiFab(BoxArray(ba).surroundingNodes(1), dmap[lev], 1,
                                  nghost, MFInfo());
    bcoeff[lev][1]->setVal(0.);

    // Create a BoxArray on z-faces.
    // z-face-based coefficient for MAC and diffusive solves
    if (bcoeff[lev][2] != nullptr)
      delete bcoeff[lev][2];

    bcoeff[lev][2] = new MultiFab(BoxArray(ba).surroundingNodes(2), dmap[lev], 1,
                                  nghost, MFInfo());
    bcoeff[lev][2]->setVal(0.);
}

void
bmx::RegridArrays (int lev)
{
    if (ooo_debug) amrex::Print() << "RegridArrays" << std::endl;

    // ********************************************************************************
    // Cell-based arrays
    // ********************************************************************************
    //
    // Note: after calling copy() using dst_ngrow, we do not need to call FillBoundary().
    //       However, we want to be sure to only use valid regions of the src, so we use
    //       src_ngrow = 0 and dst_ngrow = all
    //
    //
    int src_ngrow = 0;

    if (advect_fluid_chem_species) {
      // Gas chem_species mass fraction
      MultiFab* X_k_new = new MultiFab(grids[lev], dmap[lev],
                                       m_leveldata[lev]->X_k->nComp(),
                                       m_leveldata[lev]->X_k->nGrow(),
                                       MFInfo());
      X_k_new->setVal(0);
      X_k_new->ParallelCopy(*m_leveldata[lev]->X_k, 0, 0,
          m_leveldata[lev]->X_k->nComp(), src_ngrow,
          m_leveldata[lev]->X_k->nGrow(), geom[lev].periodicity());
      std::swap(m_leveldata[lev]->X_k, X_k_new);
      delete X_k_new;

      // Old gas chem_species mass fraction
      MultiFab* X_ko_new = new MultiFab(grids[lev], dmap[lev],
                                       m_leveldata[lev]->X_k->nComp(),
                                       m_leveldata[lev]->X_k->nGrow(),
                                       MFInfo());
      X_ko_new->setVal(0);
      X_ko_new->ParallelCopy(*m_leveldata[lev]->X_ko, 0, 0,
          m_leveldata[lev]->X_ko->nComp(), src_ngrow,
          m_leveldata[lev]->X_ko->nGrow(), geom[lev].periodicity());
      std::swap(m_leveldata[lev]->X_ko, X_ko_new);
      delete X_ko_new;

      // ChemSpecies diffusion coefficients
      MultiFab* D_k_new = new MultiFab(grids[lev], dmap[lev],
                                       m_leveldata[lev]->D_k->nComp(),
                                       m_leveldata[lev]->D_k->nGrow(),
                                       MFInfo());
      D_k_new->setVal(0);
      D_k_new->ParallelCopy(*m_leveldata[lev]->D_k, 0, 0, m_leveldata[lev]->D_k->nComp(),
                     src_ngrow, m_leveldata[lev]->D_k->nGrow(), geom[lev].periodicity());
      std::swap(m_leveldata[lev]->D_k, D_k_new);
      delete D_k_new;
    }
}
