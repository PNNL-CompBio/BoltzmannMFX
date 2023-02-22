//
//     Copyright (c) 2013 Battelle Memorial Institute
//     Licensed under modified BSD License. A copy of this license can be found
//     in the LICENSE file in the top level directory of this distribution.
//
#ifndef _BMX_ALGORITHM_H_
#define _BMX_ALGORITHM_H_

#include <AMReX_GpuUtility.H>

namespace aux {

template <class T>
class is_equal {
  public:
    AMREX_GPU_HOST_DEVICE
    is_equal (const T& value):
      m_value(value)
    {}

    AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE
    bool operator () (const T& other)
    { return other == m_value; }

  private:
    const T m_value;
};


template<class InputIterator, class UnaryPredicate>
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
bool any_of (InputIterator first, InputIterator last, UnaryPredicate pred)
{
  while (first != last) {
    if (pred(*first))
      return true;
    
    ++first;
  }

  return false; 
}

}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
int BMXfind (const int* array, const int N, const int elem)
{
  for (int n(0); n < N; n++)
    if (array[n] == elem)
      return n;

  return -1; // -1 is the InvalidIdx value
}

#endif
