#
# Macro to create tags
#
function (add_tags_targets )

   # Source directories where to look
   set(_src_dirs ${PROJECT_SOURCE_DIR}/src)

   # Add AMReX source if in superbuild mode -- don't use the alias AMReX::amrex here
   # otherwise it will always be true even if for standalone mode 
   if (TARGET amrex)
      set(_src_dirs "${_src_dirs}"  "${PROJECT_SOURCE_DIR}/subprojects/amrex/Src")
   endif ()


  # Add rule to generate TAGS -- on macOS ctags-exuberant is just ctags
  set(_tags_exe FALSE)

  find_program(_ctags_exu "ctags-exuberant")
  if (_ctags_exu)
     set(_tags_exe ctags-exuberant)
  else ()
     find_program(_ctags "ctags")
     if(_ctags)
        set(_tags_exe ctags)
     endif()
  endif ()

  if(_tags_exe)   
     add_custom_target ( tags
        COMMAND ${_tags_exe} -R    --fortran-kinds=+i  ${_src_dirs}
        COMMAND ${_tags_exe} -R -e --fortran-kinds=+i  ${_src_dirs}
        COMMENT "Generating tags files"
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/ )
     # Some file systems are case-insensitive so the above will fail:
     #  the first command would overwrite the second!
     add_custom_target ( ctags
        COMMAND ${_tags_exe} -R    --fortran-kinds=+i ${_src_dirs}
        COMMENT "Generating only ctags file"
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/ )
     add_custom_target( etags
        COMMAND ${_tags_exe} -R -e --fortran-kinds=+i ${_src_dirs}
        COMMENT "Generating only etags file"
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/ )
  endif ()


endfunction ()
