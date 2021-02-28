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
   bmx_set_species_bcs.cpp
   bmx_mf_helpers.H
   bmx_mf_helpers.cpp
   )
