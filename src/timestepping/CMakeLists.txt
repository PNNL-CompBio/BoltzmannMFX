target_include_directories( bmxcore
   PUBLIC
   $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

target_sources( bmxcore
   PRIVATE 
   bmx_evolve.cpp
   bmx_evolve_fluid.cpp
   )
