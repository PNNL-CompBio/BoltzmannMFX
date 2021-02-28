#include <bmx.H>
#include <bmx_fluid_parms.H>
#include <bmx_species_parms.H>

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
    bool need_regrid = true;

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

    if (advect_fluid_species) {
      // Gas species mass fraction
      MultiFab* X_gk_new = new MultiFab(grids[lev], dmap[lev],
                                       m_leveldata[lev]->X_gk->nComp(),
                                       m_leveldata[lev]->X_gk->nGrow(),
                                       MFInfo());
      X_gk_new->setVal(0);
      X_gk_new->ParallelCopy(*m_leveldata[lev]->X_gk, 0, 0,
          m_leveldata[lev]->X_gk->nComp(), src_ngrow,
          m_leveldata[lev]->X_gk->nGrow(), geom[lev].periodicity());
      std::swap(m_leveldata[lev]->X_gk, X_gk_new);
      delete X_gk_new;

      // Old gas species mass fraction
      MultiFab* X_gko_new = new MultiFab(grids[lev], dmap[lev],
                                       m_leveldata[lev]->X_gk->nComp(),
                                       m_leveldata[lev]->X_gk->nGrow(),
                                       MFInfo());
      X_gko_new->setVal(0);
      X_gko_new->ParallelCopy(*m_leveldata[lev]->X_gko, 0, 0,
          m_leveldata[lev]->X_gko->nComp(), src_ngrow,
          m_leveldata[lev]->X_gko->nGrow(), geom[lev].periodicity());
      std::swap(m_leveldata[lev]->X_gko, X_gko_new);
      delete X_gko_new;

      // Species diffusion coefficients
      MultiFab* D_gk_new = new MultiFab(grids[lev], dmap[lev],
                                       m_leveldata[lev]->D_gk->nComp(),
                                       m_leveldata[lev]->D_gk->nGrow(),
                                       MFInfo());
      D_gk_new->setVal(0);
      D_gk_new->ParallelCopy(*m_leveldata[lev]->D_gk, 0, 0, m_leveldata[lev]->D_gk->nComp(),
                     src_ngrow, m_leveldata[lev]->D_gk->nGrow(), geom[lev].periodicity());
      std::swap(m_leveldata[lev]->D_gk, D_gk_new);
      delete D_gk_new;
    }
}
