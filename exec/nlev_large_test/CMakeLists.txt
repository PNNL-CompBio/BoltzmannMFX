#add_custom_command(
#  OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/input_real_nlev"
#  COMMAND ${CMAKE_COMMAND} -E copy
#  ${PROJECT_SOURCE_DIR}/exec/nlev_large_test/input_real_nlev
#  ${CMAKE_CURRENT_BINARY_DIR}
#
#  DEPENDS
#  ${PROJECT_SOURCE_DIR}/exec/nlev_large_test/input_real_nlev
#)
#
#add_custom_command(
#  OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/abc_real_input.dat"
#  COMMAND ${CMAKE_COMMAND} -E copy
#  ${PROJECT_SOURCE_DIR}/exec/nlev_large_test/abc_real_input.dat
#  ${CMAKE_CURRENT_BINARY_DIR}
#
#  DEPENDS
#  ${PROJECT_SOURCE_DIR}/exec/nlev_large_test/abc_real_input.dat
#)
#
#add_custom_target(nlev_large.input
#  DEPENDS
#  ${PROJECT_SOURCE_DIR}/exec/nlev_large_test/input_real_nlev
#  ${PROJECT_SOURCE_DIR}/exec/nlev_large_test/abc_real_input.dat
#)
#
#add_dependencies(bmx nlev_large.input)

configure_file(
  ${PROJECT_SOURCE_DIR}/exec/nlev_large_test/input_real_nlev
  ${CMAKE_CURRENT_BINARY_DIR}/input_real_nlev
)

configure_file(
  ${PROJECT_SOURCE_DIR}/exec/nlev_large_test/abc_real_input.dat
  ${CMAKE_CURRENT_BINARY_DIR}/abc_real_input.dat
)
