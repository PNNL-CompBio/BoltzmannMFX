//
//     Copyright (c) 2013 Battelle Memorial Institute
//     Licensed under modified BSD License. A copy of this license can be found
//     in the LICENSE file in the top level directory of this distribution.
//
#ifndef BMX_FILLPATCH_BC_H_
#define BMX_FILLPATCH_BC_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_PhysBCFunct.H>


struct BMXForFill
{
  int probtype;

  AMREX_GPU_HOST
  constexpr BMXForFill (int a_probtype) : probtype(a_probtype) {}

  AMREX_GPU_DEVICE
  void operator() (const amrex::IntVect& /*iv*/, amrex::Array4<amrex::Real> const& /*vel*/,
                   const int /*dcomp*/, const int /*numcomp*/,
                   amrex::GeometryData const& /*geom*/, const amrex::Real /*time*/,
                   const amrex::BCRec* /*bcr*/, const int /*bcomp*/,
                   const int /*orig_comp*/) const
  {}
};

#endif
