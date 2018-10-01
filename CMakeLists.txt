cmake_minimum_required(VERSION 3.0)

include(Common)

project("orbitTools")

set(PROJECT_LIB_NAME orbitTools)

discover_variable(QD_ROOT PATH "QD math library")
discover_variable(ENABLE_QD_QD_INTEGRATION STRING "enable Qd library qd_real replacement for the double")
discover_variable(ENABLE_QD_DD_INTEGRATION STRING "enable Qd library dd_real replacement for the double")


file(GLOB_RECURSE all_headers ${CMAKE_CURRENT_LIST_DIR}/*.h*)
file(GLOB_RECURSE all_sources ${CMAKE_CURRENT_LIST_DIR}/*.c*)

add_library(${PROJECT_LIB_NAME} STATIC ${all_headers};${all_sources})

initialize_library_target_defaults(${PROJECT_LIB_NAME} "32bit")

source_groups_from_dir_list("Header Files" FILES ${CMAKE_CURRENT_LIST_DIR} *.h*)
source_groups_from_dir_list("Source Files" FILES ${CMAKE_CURRENT_LIST_DIR} *.c*)

target_include_directories(${PROJECT_LIB_NAME}
  PUBLIC
    ${CMAKE_CURRENT_LIST_DIR}/core
)

# integration with QD

add_target_subdirectory(qd ${QD_ROOT} ${CMAKE_BUILD_ROOT}/_3dparty/math/qd)

add_dependencies(${PROJECT_LIB_NAME} qd)

target_link_libraries(${PROJECT_LIB_NAME}
  PUBLIC
    qd
)

if (ENABLE_QD_QD_INTEGRATION)
  add_target_compile_definitions(${PROJECT_LIB_NAME} *
    PUBLIC
      ENABLE_QD_QD_INTEGRATION
  )
endif()
if (ENABLE_QD_DD_INTEGRATION)
  add_target_compile_definitions(${PROJECT_LIB_NAME} *
    PUBLIC
      ENABLE_QD_DD_INTEGRATION
  )
endif()
