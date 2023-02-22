#add_custom_command(
#  OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/input_fungi"
#  COMMAND ${CMAKE_COMMAND} -E copy
#  ${PROJECT_SOURCE_DIR}/exec/fungi/input_fungi
#  ${CMAKE_CURRENT_BINARY_DIR}
#
#  DEPENDS
#  ${PROJECT_SOURCE_DIR}/exec/fungi/input_fungi
#)
#
#add_custom_command(
#  OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/fungi_init_cfg.dat"
#  COMMAND ${CMAKE_COMMAND} -E copy
#  ${PROJECT_SOURCE_DIR}/exec/fungi/fungi_init_cfg.dat
#  ${CMAKE_CURRENT_BINARY_DIR}
#
#  DEPENDS
#  ${PROJECT_SOURCE_DIR}/exec/fungi/fungi_init_cfg.dat
#)
#
#add_custom_target(fungi.input
#  DEPENDS
#  ${PROJECT_SOURCE_DIR}/exec/fungi/input_fungi
#  ${PROJECT_SOURCE_DIR}/exec/fungi/fungi_init_cfg.dat
#)
#
#add_dependencies(bmx fungi.input)

configure_file(
  ${PROJECT_SOURCE_DIR}/exec/fungi/input_fungi
  ${CMAKE_CURRENT_BINARY_DIR}/input_fungi
)

configure_file(
  ${PROJECT_SOURCE_DIR}/exec/fungi/fungi_init_cfg.dat
  ${CMAKE_CURRENT_BINARY_DIR}/fungi_init_cfg.dat
)
