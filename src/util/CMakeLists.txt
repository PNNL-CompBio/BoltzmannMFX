#
#     Copyright (c) 2013 Battelle Memorial Institute
#     Licensed under modified BSD License. A copy of this license can be found
#     in the LICENSE file in the top level directory of this distribution.
#
target_include_directories(bmxcore
   PUBLIC
   $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

target_sources( bmxcore
   PRIVATE
   bmx_algorithm.H
   bmx_filcc.H
   bmx_filcc.cpp
   bmx_regrid.cpp
   bmx_fillpatch_bc.H
   bmx_fillpatch.cpp
   bmx_set_chem_species_bcs.cpp
   bmx_util.cpp
   )
