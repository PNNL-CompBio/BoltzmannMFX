#ifndef __BMX_INIT_FLUID_HPP_
#define __BMX_INIT_FLUID_HPP_

#include <bmx.H>

#include <AMReX_REAL.H>
#include <AMReX_BLFort.H>
#include <AMReX_SPACE.H>
#include <AMReX_GpuContainers.H>

void init_fluid (const amrex::Box& sbx,
                 const amrex::Box& bx,
                 const amrex::Box& domain,
                 const amrex::MFIter& mfi,
                 LevelData& ld,
                 const amrex::Real dx,
                 const amrex::Real dy,
                 const amrex::Real dz,
                 const amrex::Real xlength,
                 const amrex::Real ylength,
                 const amrex::Real zlength,
                 const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& plo,
                 const int advect_fluid_chem_species);

void init_fluid_parameters (const Box& bx,
                            const MFIter& mfi,
                            LevelData& ld,
                            const int advect_fluid_chem_species);

void set_ic (const Box& sbx,
             const Box& domain,
             const Real dx,
             const Real dy,
             const Real dz,
             FArrayBox& vel_g);

#endif
