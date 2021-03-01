target_include_directories(bmxcore
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

target_sources( bmxcore
   PRIVATE
   bmx_bc_list.H
   bmx_fluid_parms.cpp
   bmx_fluid_parms.H
   bmx_chem_species_parms.cpp
   bmx_chem_species_parms.H
   bmx_dem_parms.cpp
   bmx_dem_parms.H
   bmx_bc_parms.cpp
   bmx_bc_parms.H
   )
