//
//     Copyright (c) 2013 Battelle Memorial Institute
//     Licensed under modified BSD License. A copy of this license can be found
//     in the LICENSE file in the top level directory of this distribution.
//
#include <bmx.H>

Real
bmx::volSum ()
{
    BL_PROFILE("bmx::volSum()");
    AMREX_ALWAYS_ASSERT(finest_level <= 1);

    Real sum_across_levels(0.);
    Real sum_on_level;

    MultiFab* mask;
    if (finest_level > 0)
        mask = build_fine_mask();

    for (int lev = 0; lev <= finest_level; lev++) 
    { 
        const MultiFab& vf = *get_vf_new_const()[lev];
        if (lev == finest_level)
        {
            sum_on_level = amrex::ReduceSum(vf, 0,
                [] AMREX_GPU_HOST_DEVICE (Box const & bx,
                                              Array4<const Real> const & vfrc)
                {
                  Real dm = 0.0;
        
                  amrex::Loop(bx, [vfrc,&dm] (int i, int j, int k) noexcept
                      { dm += vfrc(i,j,k); });
        
                  return dm;
                });
        } else {
            sum_on_level = amrex::ReduceSum(vf, *mask, 0,
                [] AMREX_GPU_HOST_DEVICE (Box const & bx,
                                          Array4<const Real> const & vfrc,
                                          Array4<const Real> const & msk)
                {
                  Real dm = 0.0;
        
                  amrex::Loop(bx, [vfrc,msk,&dm] (int i, int j, int k) noexcept
                      { dm += vfrc(i,j,k) * msk(i,j,k); });
        
                  return dm;
                });
        }
        ParallelDescriptor::ReduceRealSum(sum_on_level);

        Real dx = geom[lev].CellSize(0);
        Real dy = geom[lev].CellSize(1);
        Real dz = geom[lev].CellSize(2);
        Real cell_vol = dx * dy * dz;

        sum_across_levels += sum_on_level * cell_vol;
    }

    return sum_across_levels;
}

Real
bmx::volWgtSum (Vector<const MultiFab*> const& mfs, int comp )
{
    BL_PROFILE("bmx::volWgtSum()");
    AMREX_ALWAYS_ASSERT(finest_level <= 1);
    Real sum_across_levels(0.);
    Real sum_on_level;

    MultiFab* mask;
    if (finest_level > 0)
        mask = build_fine_mask();

    for (int lev = 0; lev <= finest_level; lev++) 
    { 
        if (lev == finest_level)
        {
            const MultiFab& vf = *(get_vf_new_const()[lev]);
            sum_on_level = amrex::ReduceSum(*mfs[lev],vf,0,
                [comp] AMREX_GPU_HOST_DEVICE (Box const & bx,
                                              Array4<const Real> const & qty,
                                              Array4<const Real> const & vfrc)
                {
                  Real dm = 0.0;
        
                  amrex::Loop(bx, [qty,vfrc,comp,&dm] (int i, int j, int k) noexcept
                      { dm += qty(i,j,k,comp) * vfrc(i,j,k); });

                  return dm;
                });
        } else {
            const MultiFab& vf = *(get_vf_new_const()[lev]);
            sum_on_level = amrex::ReduceSum(*mfs[lev],vf,*mask,0,
                [comp] AMREX_GPU_HOST_DEVICE (Box const & bx,
                                              Array4<const Real> const & qty,
                                              Array4<const Real> const & vfrc,
                                              Array4<const Real> const & msk)
                {
                  Real dm = 0.0;
        
                  amrex::Loop(bx, [qty,vfrc,msk,comp,&dm] (int i, int j, int k) noexcept
                      { dm += qty(i,j,k,comp) * vfrc(i,j,k) * msk(i,j,k); });

                  return dm;
                });
        }
        ParallelDescriptor::ReduceRealSum(sum_on_level);

        Real dx = geom[lev].CellSize(0);
        Real dy = geom[lev].CellSize(1);
        Real dz = geom[lev].CellSize(2);
        Real cell_vol = dx * dy * dz;

        sum_across_levels += sum_on_level * cell_vol;
    } 

    return sum_across_levels;
}

MultiFab*
bmx::build_fine_mask()
{
    int crse_level = 0;

    if (fine_mask != 0) return fine_mask;

    BoxArray baf = grids[crse_level+1];
    baf.coarsen(refRatio(crse_level));

    const BoxArray& bac = grids[crse_level];
    fine_mask = new MultiFab(bac,dmap[crse_level], 1,0);
    fine_mask->setVal(1.0);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*fine_mask,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        FArrayBox& fab = (*fine_mask)[mfi];

        std::vector< std::pair<int,Box> > isects = baf.intersections(fab.box());

        for (int ii = 0; ii < isects.size(); ii++)
        {
            fab.setVal<RunOn::Device>(0.0,isects[ii].second,0);
        }
    }
    return fine_mask;
}
