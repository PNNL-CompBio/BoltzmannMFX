target_include_directories( bmxcore
   PUBLIC
   $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

target_sources( bmxcore
   PRIVATE
   bmx_chem.H
   bmx_chem_K.H
   bmx_cell_interaction.H
   bmx_cell_interaction_K.H
   bmx_chem.cpp
   bmx_cell_interaction.cpp
   )

