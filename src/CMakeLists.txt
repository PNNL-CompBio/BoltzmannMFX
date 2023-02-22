#
#     Copyright (c) 2013 Battelle Memorial Institute
#     Licensed under modified BSD License. A copy of this license can be found
#     in the LICENSE file in the top level directory of this distribution.
#
#

add_library(bmxcore OBJECT)

#
# Add sources and header from current dir and all sub dirs
#

# Uncategorized
target_include_directories(bmxcore
   PUBLIC
   $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

# IMPORTANT: DO NOT ADD main.cpp HERE!
target_sources(bmxcore
   PRIVATE
   bmx.cpp
   #  main.cpp <---- DO NOT ADD THIS HERE
   bmx.H
   )

# I/0
add_subdirectory(io)

# DES
add_subdirectory(interpolation)
add_subdirectory(deposition)
add_subdirectory(des)

# Diffusion
add_subdirectory(diffusion)

# LevelData
add_subdirectory(leveldata)

# Setup
add_subdirectory(setup)

# Utility functions
add_subdirectory(util)

# Global module files
add_subdirectory(mods)

# Time evolution functions
add_subdirectory(timestepping)

# Cell chemistry functions
add_subdirectory(chemistry)

#
# Link to imported AMReX
#
target_link_libraries( bmxcore PUBLIC AMReX::amrex )

#
# Require C++14 or newer and no compiler extensions.
#
set_target_properties(bmxcore PROPERTIES CXX_EXTENSIONS OFF)
target_compile_features(bmxcore PUBLIC cxx_std_14)

#
# If BMX_CUDA, set CUDA-specific properties and
# make C++ files be compiled as CUDA sources
#
if (BMX_CUDA)
   setup_target_for_cuda_compilation(bmxcore)

   #
   # Add manually curand
   # CMake >= 3.17 provides the module FindCUDAToolkit to do this natively.
   #
   find_library(CURAND curand PATHS ${CMAKE_CUDA_IMPLICIT_LINK_DIRECTORIES})
   target_link_libraries(bmxcore PUBLIC ${CURAND})
endif ()

#
# Set the default flags
#
if ( NOT CMAKE_CXX_FLAGS )
   target_link_libraries(bmxcore PUBLIC AMReX::Flags_CXX)
endif ()

if (BMX_FPE)
   target_link_libraries(bmxcore PUBLIC AMReX::Flags_FPE)
endif ()

#
# Add AMReX "build info" capabilities
#
include(AMReXBuildInfo)
generate_buildinfo(bmxcore ${CMAKE_SOURCE_DIR} REQUIRED)
