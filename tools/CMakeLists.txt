add_executable( fjoin_par EXCLUDE_FROM_ALL C++/fjoin_par.cpp )
target_link_libraries( fjoin_par AMReX::amrex )

if (BMX_CUDA)
   setup_target_for_cuda_compilation(fjoin_par)
endif ()

# Define full path to fjoin_par
set(FJOIN_PAR ${CMAKE_CURRENT_BINARY_DIR}/fjoin_par CACHE INTERNAL "")
