#
#     Copyright (c) 2013 Battelle Memorial Institute
#     Licensed under modified BSD License. A copy of this license can be found
#     in the LICENSE file in the top level directory of this distribution.
#
target_include_directories( bmxcore
   PUBLIC
   $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

target_sources( bmxcore
   PRIVATE
   bmx_des_K.H
   bmx_pc.cpp
   bmx_pc.H
   bmx_pc_init.cpp
   bmx_pc_interaction.cpp
   bmx_pc_write.cpp
   bmx_calc_volume_fraction.cpp
   bmx_calc_txfr.cpp
   bmx_split_particles.cpp
   bmx_pc_xchng.cpp
   bmx_pc_fusion.cpp
   )

