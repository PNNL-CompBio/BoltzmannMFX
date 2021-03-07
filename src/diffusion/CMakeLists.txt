target_include_directories(bmxcore
   PUBLIC
   $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

target_sources( bmxcore
   PRIVATE
   bmx_define_coeffs_on_faces.cpp
   bmx_diffusion_op.cpp
   bmx_diffusion_op.H
   bmx_diffuse_chem_species.cpp
   )
