#ifndef _BMX_BCLIST_H_
#define _BMX_BCLIST_H_

// Boundary conditions list auxiliary structure
//
class BCList {
  public:
    BCList (const int undefined = 0,
            const int ig = 9,
            const int pinf = 10,
            const int pout = 11,
            const int minf = 20,
            const int nsw = 100,
            const int eb  = 111)
      : und_(undefined)
      , ig_(ig)
      , pinf_(pinf)
      , pout_(pout)
      , minf_(minf)
      , nsw_(nsw)
      , eb_(eb)
    {}

    int get_undefined () const {return und_;};
    int get_ig () const {return ig_;};

    AMREX_GPU_HOST_DEVICE AMREX_INLINE int get_pinf () const {return pinf_;};
    AMREX_GPU_HOST_DEVICE AMREX_INLINE int get_pout () const {return pout_;};
    AMREX_GPU_HOST_DEVICE AMREX_INLINE int get_minf () const {return minf_;};
    AMREX_GPU_HOST_DEVICE AMREX_INLINE int get_nsw  () const {return nsw_ ;};
    AMREX_GPU_HOST_DEVICE AMREX_INLINE int get_eb   () const {return eb_  ;};

  private:
    const int und_;
    const int ig_;
    const int pinf_;
    const int pout_;
    const int minf_;
    const int nsw_;
    const int eb_;
};

#endif
