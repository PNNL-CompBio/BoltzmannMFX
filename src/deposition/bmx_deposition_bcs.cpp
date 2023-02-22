//
//     Copyright (c) 2013 Battelle Memorial Institute
//     Licensed under modified BSD License. A copy of this license can be found
//     in the LICENSE file in the top level directory of this distribution.
//
#include <bmx.H>

void bmx::bmx_deposition_bcs (int lev, amrex::MultiFab& filled_mf)
{
  BL_PROFILE("bmx::bmx_deposition_bcs_scalar()");

  Box domain(Geom(lev).Domain());

  const int minf = bc_list.get_minf();

  for (MFIter mfi(filled_mf, false); mfi.isValid(); ++mfi) {

    const Box& sbx = filled_mf[mfi].box();

    // Extract the lower and upper boundaries of Box sbx and Domain
    const IntVect sbx_lo(sbx.loVect()), sbx_hi(sbx.hiVect());
    const amrex::Dim3& dom_lo = amrex::lbound(domain);
    const amrex::Dim3& dom_hi = amrex::ubound(domain);

    // Create a 2D Box collapsing sbx on x-direction
    IntVect sbx_yz_hi(sbx.hiVect());
    sbx_yz_hi[0] = sbx_lo[0];
    const Box sbx_yz(sbx_lo, sbx_yz_hi);

    // Create a 2D Box collapsing sbx on y-direction
    IntVect sbx_xz_hi(sbx.hiVect());
    sbx_xz_hi[1] = sbx_lo[1];
    const Box sbx_xz(sbx_lo, sbx_xz_hi);

    // Create a 2D Box collapsing sbx on z-direction
    IntVect sbx_xy_hi(sbx.hiVect());
    sbx_xy_hi[2] = sbx_lo[2];
    const Box sbx_xy(sbx_lo, sbx_xy_hi);

    Array4<int> const& bc_ilo_type = bc_ilo[lev]->array();
    Array4<int> const& bc_ihi_type = bc_ihi[lev]->array();
    Array4<int> const& bc_jlo_type = bc_jlo[lev]->array();
    Array4<int> const& bc_jhi_type = bc_jhi[lev]->array();
    Array4<int> const& bc_klo_type = bc_klo[lev]->array();
    Array4<int> const& bc_khi_type = bc_khi[lev]->array();

    const int ncomp = filled_mf.nComp();

    Array4<Real> const& vol = filled_mf.array(mfi);

    if(sbx_lo[0] < dom_lo.x or sbx_hi[0] > dom_hi.x or
       sbx_lo[1] < dom_lo.y or sbx_hi[1] > dom_hi.y or
       sbx_lo[2] < dom_lo.z or sbx_hi[2] > dom_hi.z) {

      const int sbx_yz_npoints = sbx_yz.numPts();
      const auto sbx_yz_lo = amrex::lbound(sbx_yz);
      const auto sbx_yz_len = amrex::length(sbx_yz);

      const int sbx_xz_npoints = sbx_xz.numPts();
      const auto sbx_xz_lo = amrex::lbound(sbx_xz);
      const auto sbx_xz_len = amrex::length(sbx_xz);

      const int sbx_xy_npoints = sbx_xy.numPts();
      const auto sbx_xy_lo = amrex::lbound(sbx_xy);
      const auto sbx_xy_len = amrex::length(sbx_xy);

      const int npoints = amrex::max(sbx_yz_npoints,sbx_xz_npoints,sbx_xy_npoints);

      ParallelFor(npoints, [sbx_xy_npoints,sbx_yz_npoints,sbx_xz_npoints,
          sbx_xy_len,sbx_yz_len,sbx_xz_len,sbx_xy_lo,sbx_yz_lo,sbx_xz_lo,
          dom_lo,dom_hi,bc_ilo_type,bc_ihi_type,bc_jlo_type,bc_jhi_type,
          bc_klo_type,bc_khi_type,vol,ncomp,minf,sbx_hi,sbx_lo]
        AMREX_GPU_DEVICE (int idx) noexcept
      {
        if(idx < sbx_yz_npoints)
        {
          int k = idx / (sbx_yz_len.y);
          int j = idx - k*(sbx_yz_len.y);

          j += sbx_yz_lo.y;
          k += sbx_yz_lo.z;

          if(sbx_lo[0] < dom_lo.x) {
            int i = dom_lo.x;

            // amrex::Print() << "DOING I J K " << IntVect(i,j,k) << std::endl;
            const int bct_ilo = bc_ilo_type(i-1,j,k,0);
            if(bct_ilo == minf)
            {
              for(int n(0); n < ncomp; n++) {
                vol(i,  j,k,n) += vol(i-1,j,k,n);
                vol(i-1,j,k,n) = 0;
              }
            }
          }
          if(sbx_hi[0] > dom_hi.x) {
            int i = dom_hi.x;

            const int bct_ihi = bc_ihi_type(i+1,j,k,0);
            if(bct_ihi == minf)
            {
              for(int n(0); n < ncomp; n++) {
                vol(i,  j,k,n) += vol(i+1,j,k,n);
                vol(i+1,j,k,n) = 0;
              }
            }
          }
        }

        if(idx < sbx_xz_npoints)
        {
          int k = idx / (sbx_xz_len.x);
          int i = idx - k*(sbx_xz_len.x);

          i += sbx_xz_lo.x;
          k += sbx_xz_lo.z;

          if(sbx_lo[1] < dom_lo.y) {
            int j = dom_lo.y;

            const int bct_jlo = bc_jlo_type(i,j-1,k,0);
            if(bct_jlo == minf)
            {
              for(int n(0); n < ncomp; n++) {
                vol(i,j  ,k,n) += vol(i,j-1,k,n);
                vol(i,j-1,k,n) = 0;
              }
            }
          }
          if(sbx_hi[1] > dom_hi.y) {
            int j = dom_hi.y;

            const int bct_jhi = bc_jhi_type(i,j+1,k,0);
            if( bct_jhi == minf)
            {
              for(int n(0); n < ncomp; n++) {
                vol(i,j  ,k,n) += vol(i,j+1,k,n);
                vol(i,j+1,k,n) = 0;
              }
            }
          }
        }

        if(idx < sbx_xy_npoints)
        {
          int j = idx / (sbx_xy_len.x);
          int i = idx - j*(sbx_xy_len.x);

          i += sbx_xy_lo.x;
          j += sbx_xy_lo.y;

          if(sbx_lo[2] < dom_lo.z) {
            int k = dom_lo.z;

            const int bct_klo = bc_klo_type(i,j,k-1,0);
            if(bct_klo == minf)
            {
              for(int n(0); n < ncomp; n++) {
                vol(i,j,k  ,n) += vol(i,j,k-1,n);
                vol(i,j,k-1,n) = 0;
              }
            }
          }
          if(sbx_hi[2] > dom_hi.z) {
            int k = dom_hi.z;

            const int bct_khi = bc_khi_type(i,j,k+1,0);
            if(bct_khi == minf)
            {
              for(int n(0); n < ncomp; n++) {
                vol(i,j,k  ,n) += vol(i,j,k+1,n);
                vol(i,j,k+1,n) = 0;
              }
            }
          }
        }
      });
    }
  }
}
