//
//     Copyright (c) 2013 Battelle Memorial Institute
//     Licensed under modified BSD License. A copy of this license can be found
//     in the LICENSE file in the top level directory of this distribution.
//
#ifndef _BMX_BCLIST_H_
#define _BMX_BCLIST_H_

// Boundary conditions list auxiliary structure
//
class BCList {
  public:
    BCList (const int undefined = 0,
            const int ig = 9,
            const int minf = 20,
            const int nsw = 100)
      : und_(undefined)
      , ig_(ig)
      , minf_(minf)
      , nsw_(nsw)
    {}

    int get_undefined () const {return und_;};
    int get_ig () const {return ig_;};

    AMREX_GPU_HOST_DEVICE AMREX_INLINE int get_minf () const {return minf_;};
    AMREX_GPU_HOST_DEVICE AMREX_INLINE int get_nsw  () const {return nsw_ ;};

  private:
    const int und_;
    const int ig_;
    const int minf_;
    const int nsw_;
};

#endif
