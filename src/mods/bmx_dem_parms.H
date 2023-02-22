//
//     Copyright (c) 2013 Battelle Memorial Institute
//     Licensed under modified BSD License. A copy of this license can be found
//     in the LICENSE file in the top level directory of this distribution.
//
#ifndef BMX_DEM_PARMS_H_
#define BMX_DEM_PARMS_H_

#include <AMReX_REAL.H>
#include <AMReX_Gpu.H>
#include <AMReX_Arena.H>
#include <AMReX_Vector.H>
#include <AMReX_Array.H>

namespace DEM
{
    extern int NPHASE;

    extern int solve;

    extern amrex::Real dtsolid;

    template <typename ARRAY_TYPE>
    struct DeviceArray
    {
        typedef ARRAY_TYPE array_type;

        template <typename... Is>
        AMREX_GPU_DEVICE AMREX_FORCE_INLINE
        amrex::Real& operator() (Is... index) const noexcept {
            return (*a)(index...);
        }

        void alloc () {
#ifdef AMREX_USE_GPU
            auto p = (ARRAY_TYPE*)amrex::The_Device_Arena()->alloc(sizeof(ARRAY_TYPE));
            amrex::single_task(amrex::Gpu::gpuStream(), [p] AMREX_GPU_DEVICE
            {
                new (p) ARRAY_TYPE;
            });
            amrex::Gpu::synchronize();
            a = p;
#else
            a = new ARRAY_TYPE;
#endif
        }

        void free () {
#ifdef AMREX_USE_GPU
            amrex::single_task(amrex::Gpu::gpuStream(), [p=a] AMREX_GPU_DEVICE
            {
                p->~ARRAY_TYPE();
            });
            amrex::Gpu::synchronize();
            amrex::The_Device_Arena()->free(a);
#else
            delete a;
#endif
        }

        ARRAY_TYPE* arrayPtr () { return a; }

    private:
        ARRAY_TYPE * a = nullptr;
    };

    extern amrex::Real small_number;
    extern amrex::Real large_number;
    extern amrex::Real eps;

    extern amrex::Real neighborhood;

    void Initialize ();
    void Finalize ();

}

#endif
