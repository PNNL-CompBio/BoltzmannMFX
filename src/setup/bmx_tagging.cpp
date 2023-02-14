//
//     Copyright (c) 2013 Battelle Memorial Institute
//     Licensed under modified BSD License. A copy of this license can be found
//     in the LICENSE file in the top level directory of this distribution.
//
#include <AMReX_ParmParse.H>

#include <bmx.H>

void bmx::ErrorEst (int lev, TagBoxArray & tags, Real /*time*/, int /*ngrow*/)
{
    if (ooo_debug) amrex::Print() << "ErrorEst" << std::endl;

    static bool first = true;
    static Vector<Real> rhoerr_v, gradrhoerr_v;

    static bool tag_region;

    if (first) 
    {
        first = false;
        ParmParse pp("bmx");

        tag_region_lo.resize(3);
        tag_region_hi.resize(3);

        tag_region = false;
        pp.query("tag_region", tag_region);

        if (tag_region)
        {
            pp.getarr("tag_region_lo", tag_region_lo);
            pp.getarr("tag_region_hi", tag_region_hi);
        }
    }

    const auto   tagval = TagBox::SET;

    if (tag_region) {

      Real xlo = tag_region_lo[0];
      Real ylo = tag_region_lo[1];
      Real xhi = tag_region_hi[0];
      Real yhi = tag_region_hi[1];
      Real zlo = tag_region_lo[2];
      Real zhi = tag_region_hi[2];

      const Real l_dx = geom[lev].CellSize(0);
      const Real l_dy = geom[lev].CellSize(1);
      const Real l_dz = geom[lev].CellSize(2);

      MultiFab dummy(grids[lev],dmap[lev],1,0);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(dummy,TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {
        Box const& bx = mfi.tilebox();
        auto const& tag = tags.array(mfi);

            amrex::ParallelFor(bx,
            [xlo, xhi, ylo, yhi, zlo, zhi, l_dx, l_dy, l_dz,tagval, tag]
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                 Real x = (i+0.5)*l_dx;
                 Real y = (j+0.5)*l_dy;
                 Real z = (k+0.5)*l_dz;

                 // Tag if we are inside the specified box
                 if (x >= xlo && x <= xhi && y >= ylo && y <= yhi && z >= zlo && z <= zhi)
                 {
                    tag(i,j,k) = tagval;
                 }
            });
      } // MFIter
    } // if (tag_region)
}
