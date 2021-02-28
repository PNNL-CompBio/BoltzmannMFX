#
# Find git info
#
macro( get_git_info ) # EXTRA ARGS: branch commit

  # Find branch
  execute_process(
    COMMAND git branch
    COMMAND grep \*
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    OUTPUT_VARIABLE out
    ERROR_VARIABLE  err
    )

  if(err)
    message(WARNING "Failing to retrieve BMX  Git branch")
  else()
    string( REPLACE "*" "" out ${out} )
    string( STRIP ${out} out )
    message( STATUS "BMX  branch: ${out}" )
    if( ${ARGC} GREATER 0 ) # branch
      set( ${ARGV0} ${out} )
    endif()
  endif()

  # Find commit
  execute_process(
    COMMAND git rev-parse HEAD
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    OUTPUT_VARIABLE out
    ERROR_VARIABLE  err
    )

  if(err)
    message(WARNING "Failing to retrieve BMX  Git commit")
  else()
    string(STRIP ${out} out)
    message( STATUS "BMX  commit: ${out}" )
    if( ${ARGC} EQUAL 2 ) # commit
      set( ${ARGV1} ${out} )
    endif()
  endif()

  unset(out)

endmacro()

#
# Print Configuration Summary
#
function( print_bmx_configuration_summary bmx_libname )

  if(NOT TARGET ${bmx_libname})
    message(AUTHOR_WARNING "Target ${bmx_libname} is not defined.")
    return()
  endif()

  if(NOT TARGET AMReX::amrex)
    message(AUTHOR_WARNING "Target ${bmx_libname} is not defined.")
    return()
  endif()

  string(TOUPPER "${CMAKE_BUILD_TYPE}" CMAKE_BUILD_TYPE_UPPERCASE)


  # Get AMReX cmake functions
  include(AMReXGenexHelpers)
  include(AMReXTargetHelpers)

  get_target_properties_flattened(${bmx_libname} _includes _defines _flags _link_line)

  set(_lang CXX)
  set(_prop _includes _defines _flags _link_line )


  # Loop over all combinations of language and property and extract
  # what you need
  foreach( _p IN LISTS _prop )
     foreach( _l IN LISTS _lang )

        string(TOLOWER ${_l} _ll) # Lower case language name

        # _${_ll}${_p} is a variable named as _lang_property,
        # both lower case.
        set(_${_ll}${_p} "${${_p}}")
        eval_genex( _${_ll}${_p} ${_l} ${CMAKE_${_l}_COMPILER_ID}
           COMP_VERSION ${CMAKE_${_l}_COMPILER_VERSION}
           CONFIG       ${CMAKE_BUILD_TYPE}
           INTERFACE    BUILD)

        if (_${_ll}${_p})

           list(REMOVE_DUPLICATES _${_ll}${_p})

           if ("${_p}" STREQUAL "_defines")
              string(REPLACE ";" " -D" _${_ll}${_p} "-D${_${_ll}${_p}}")
           elseif ("${_p}" STREQUAL "_includes")
              string(REPLACE ";" " -I" _${_ll}${_p} "-I${_${_ll}${_p}}")
           else()
              string(REPLACE ";" " " _${_ll}${_p} "${_${_ll}${_p}}")
           endif ()

        endif ()

     endforeach()
  endforeach ()


  #
  # Config summary
  #
  message( STATUS "BMX configuration summary: " )
  message( STATUS "  C++ defines           = ${_cxx_defines}" )
  message( STATUS "  C++ compiler          = ${CMAKE_CXX_COMPILER}" )
  message( STATUS "  C++ flags             = ${CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE_UPPERCASE}} ${CMAKE_CXX_FLAGS} ${_cxx_flags}" )
  if (BMX_CUDA)
     message( STATUS "  CUDA flags            = ${CMAKE_CUDA_FLAGS_${CMAKE_BUILD_TYPE_UPPERCASE}} ${CMAKE_CUDA_FLAGS}" )
  endif ()
  message( STATUS "   BMX includes         = ${_cxx_includes}" )
  message( STATUS "   BMX link line        = ${_cxx_link_line}" )

endfunction()
