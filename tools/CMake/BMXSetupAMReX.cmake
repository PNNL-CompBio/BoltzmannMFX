set(AMREX_MINIMUM_VERSION 20.12.1 CACHE INTERNAL "Minimum required AMReX version")

# We first check if we can find an AMReX installation.
# If so, we proceed with STANDALONE mode
# If not, we proceed with SUPERBUILD MODE
find_package( AMReX ${AMREX_MINIMUM_VERSION} CONFIG QUIET )

#
# ~~~~~~~~~~~~~~~~~~~~~~~ STANDALONE MODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
if (AMReX_FOUND)

   # We are in this branch because an AMReX installation has been found.
   # In this scenario, we don't know how AMReX has been built. We only know
   # that some AMReX components, namely 3D DOUBLE PARTICLES PDOUBLE PFTOOLS
   # EB LSOLVERS  MUST be present in the installation.
   set(AMREX_REQUIRED_COMPONENTS 3D DOUBLE PARTICLES PDOUBLE EB LSOLVERS PFTOOLS)

   if (BMX_MPI)
      list(APPEND AMREX_REQUIRED_COMPONENTS MPI)
      if (BMX_MPI_THREAD_MULTIPLE)
         list(APPEND AMREX_REQUIRED_COMPONENTS MPI_THREAD_MULTIPLE)
      endif ()
   endif ()
   if (BMX_OMP)
      list(APPEND AMREX_REQUIRED_COMPONENTS OMP)
   endif ()
   if (NOT BMX_GPU_BACKEND STREQUAL "NONE")
      list(APPEND AMREX_REQUIRED_COMPONENTS ${BMX_GPU_BACKEND})
   endif ()
   if (BMX_HYPRE)
      list(APPEND AMREX_REQUIRED_COMPONENTS HYPRE)
   endif ()


   # We now check again for the AMReX package.
   # This time we mark AMReX + its required components as REQUIRED.
   # If the AMReX installation does not contain all the required components,
   # the configuration step stops with an error message.
   find_package(AMReX ${AMREX_MINIMUM_VERSION} CONFIG
      REQUIRED ${AMREX_REQUIRED_COMPONENTS}
      )

   message(STATUS "AMReX found: configuration file located at ${AMReX_DIR}")

   # IS it worth checking this???
   if ( NOT ( "${CMAKE_BUILD_TYPE}" STREQUAL "${AMReX_BUILD_TYPE}" ) )
      message (WARNING "BMX build type (${CMAKE_BUILD_TYPE}) type does not match AMReX build type (${AMReX_BUILD_TYPE})")
   endif ()

   # We load this here so we have the CUDA helper functions
   # available everywhere we need it
   if (BMX_CUDA)
      include(AMReXTargetHelpers)
   endif ()

   #
   # Define full path to fextract
   #
   get_target_property(_config AMReX::fextract IMPORTED_CONFIGURATIONS)
   get_target_property(FEXTRACT AMReX::fextract IMPORTED_LOCATION_${_config})
   set(FEXTRACT ${FEXTRACT} CACHE INTERNAL "")
   unset(_config)

else ()
#
# ~~~~~~~~~~~~~~~~~~~~~~~ SUPERBUILD MODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#

   message(STATUS "No AMReX installation found: cloning AMReX repo")

   if (NOT EXISTS  "${PROJECT_SOURCE_DIR}/.git")
      message(FATAL_ERROR
         "${PROJECT_SOURCE_DIR} is not a Git repo: missing .git")
   endif ()

   set(AMREX_SRC_DIR "${PROJECT_SOURCE_DIR}/subprojects/amrex"
     CACHE INTERNAL "Path to AMReX source (submodule)")

   if (NOT EXISTS "${AMREX_SRC_DIR}/.git")
      message(STATUS "Initializing git submodule for AMReX")

      find_package(Git REQUIRED)

      execute_process(
         COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive subprojects/amrex
         WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
         RESULT_VARIABLE GIT_SUBMOD_RESULT
         )

      if ( NOT GIT_SUBMOD_RESULT EQUAL "0")
         message(FATAL_ERROR "git submodule update --init failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
      endif()

      unset(GIT_SUBMOD_RESULT)

   endif ()


   set(AMReX_SPACEDIM              3)
   set(AMReX_PRECISION             DOUBLE)
   set(AMReX_FORTRAN               OFF)
   set(AMReX_MPI                   ${BMX_MPI})
   set(AMReX_MPI_THREAD_MULTIPLE   ${BMX_MPI_THREAD_MULTIPLE})
   set(AMReX_OMP                   ${BMX_OMP})
   set(AMReX_GPU_BACKEND           ${BMX_GPU_BACKEND} CACHE STRING "" FORCE)
   set(AMReX_PARTICLES             ON)
   set(AMReX_PARTICLES_PRECISION   DOUBLE)
   set(AMReX_EB                    ON)
   set(AMReX_LINEAR_SOLVERS        ON)
   set(AMReX_HYPRE                 ${BMX_HYPRE})
   set(AMReX_BUILD_TUTORIALS       OFF)
   set(AMReX_PLOTFILE_TOOLS        ON)

   list(APPEND CMAKE_MODULE_PATH ${AMREX_SRC_DIR}/Tools/CMake)

   # If CUDA is required, enable the language BEFORE adding the AMReX directory
   # Since AMReX_SetupCUDA has an include guard, it will be included only once here.
   # The reason for enabling CUDA before adding the AMReX subdirectory is that
   # the top-most directory needs to setup the CUDA language before a CUDA-enabled target
   # from a sub-project is included via add_subdirectory.
   # IMPORTANT: if you don't do this, AMReX will perform this step in a sub-scope and therefore
   # it will not setup CUDA here!
   if(BMX_CUDA)
      include(AMReX_SetupCUDA)
   endif ()

   # Append -w to AMREX flags (but not BMX flags)
   #
   # In order to enable build flags for BMX compiler warnings but not not
   # enable AMReX compiler warnings, we save the CMAKE_CXX_FLAGS value set for
   # BMX, append -w to disable all warnings (at this point where AMREX subdir
   # is added), then restore the previous CMAKE_CXX_FLAGS value for later when
   # CMake reaches the point in the build scripts where the BMX build targets
   # are defined.
   set(CMAKE_CXX_FLAGS_SAVE "${CMAKE_CXX_FLAGS}")
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -w")
   add_subdirectory(${AMREX_SRC_DIR})

   # Restore CMAKE_CXX_FLAGS (without -w) for BMX
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_SAVE}")

   #
   # Define full path to fextract
   #
   get_target_property(FEXTRACT AMReX::fextract BINARY_DIR)
   set(FEXTRACT ${FEXTRACT}/fextract CACHE INTERNAL "")

endif ()
