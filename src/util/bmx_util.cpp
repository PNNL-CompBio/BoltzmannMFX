#include <bmx.H>

Real
bmx::volSum (int lev, const MultiFab& vf, bool local) const
{
    BL_PROFILE("bmx::volWgtSum()");

    Real sum = amrex::ReduceSum(vf, 0,
        [] AMREX_GPU_HOST_DEVICE (Box const & bx,
                                      Array4<const Real> const & vfrc)
        {
          Real dm = 0.0;

          amrex::Loop(bx, [vfrc,&dm] (int i, int j, int k) noexcept
              { dm += vfrc(i,j,k); });

          return dm;
        });

    if (!local)
        ParallelDescriptor::ReduceRealSum(sum);

    return sum;
}

Real
bmx::volWgtSum (int lev, const MultiFab& vf, const MultiFab& mf, int comp, bool local) const
{
    BL_PROFILE("bmx::volWgtSum()");

    Real sum = amrex::ReduceSum(mf, vf, 0,
        [comp] AMREX_GPU_HOST_DEVICE (Box const & bx,
                                      Array4<const Real> const & qty,
                                      Array4<const Real> const & vfrc)
        {
          Real dm = 0.0;

          amrex::Loop(bx, [qty,vfrc,comp,&dm] (int i, int j, int k) noexcept
              { dm += qty(i,j,k,comp) * vfrc(i,j,k); });

          return dm;
        });

    if (!local)
        ParallelDescriptor::ReduceRealSum(sum);

    return sum;
}

