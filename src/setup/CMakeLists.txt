#
#     Copyright (c) 2013 Battelle Memorial Institute
#     Licensed under modified BSD License. A copy of this license can be found
#     in the LICENSE file in the top level directory of this distribution.
#
target_include_directories(bmxcore
   PUBLIC
   $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
   )

target_sources( bmxcore
   PRIVATE 
   bmx_calc_cell.H
   bmx_calc_cell.cpp
   bmx_calc_fluid_coeffs.H
   bmx_calc_fluid_coeffs.cpp
   bmx_init_fluid.H
   bmx_init_fluid.cpp
   bmx_init.cpp
   bmx_arrays.cpp
   bmx_setup_solvers.cpp
   bmx_set_bc_type.cpp
   bmx_tagging.cpp
   )
