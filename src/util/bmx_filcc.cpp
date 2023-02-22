//
//     Copyright (c) 2013 Battelle Memorial Institute
//     Licensed under modified BSD License. A copy of this license can be found
//     in the LICENSE file in the top level directory of this distribution.
//
#include <AMReX_FilCC_3D_C.H>

#include <bmx_filcc.H>

namespace bmx_aux {

using namespace amrex;

void filcc (Real* data,
            const int* lo,
            const int* hi,
            const int* dom_lo,
            const int* dom_hi,
            const Real* /*dx*/,
            const Real* /*grd_lo*/,
            const Real* /*time*/,
            const int* bc)
{
  Box domain(IntVect(dom_lo[0], dom_lo[1], dom_lo[2]),
             IntVect(dom_hi[0], dom_hi[1], dom_hi[2]));

  Dim3 data_begin, data_end;
  data_begin.x = lo[0];
  data_begin.y = lo[1];
  data_begin.z = lo[2];
  data_end.x = hi[0];
  data_end.y = hi[1];
  data_end.z = hi[2];

  Array4<Real> q(data, data_begin, data_end, 1);

  Box bx(IntVect(lo[0], lo[1], lo[2]), IntVect(hi[0], hi[1], hi[2]));

  const int bc_lo[3] = {bc[0], bc[1], bc[2]};
  const int bc_hi[3] = {bc[3], bc[4], bc[5]};
  amrex::BCRec bcr(bc_lo, bc_hi);

  const int dcomp = 0;
  const int numcomp = 1;
  const int bcomp = 0;

  amrex::ParallelFor(bx,
    [q,bcr,domain]
    AMREX_GPU_DEVICE (int i, int j, int k) noexcept
  {
    IntVect iv(i,j,k);
    FilccCell()(iv, q, dcomp, numcomp, domain, &bcr, bcomp);
  });
}

} // end namespace bmx_aux
