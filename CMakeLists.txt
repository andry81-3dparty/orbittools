cmake_minimum_required(VERSION 3.0)

include(tacklelib/Project)

project("orbitTools")

set(PROJECT_LIB_NAME orbitTools)

tkl_check_var(OPTIONAL PATH QD_ROOT)
tkl_check_var(OPTIONAL . ENABLE_QD_INTEGRATION)
tkl_check_var(OPTIONAL . ENABLE_QD_QD_INTEGRATION)
tkl_check_var(OPTIONAL . ENABLE_QD_DD_INTEGRATION)


file(GLOB_RECURSE all_headers ${CMAKE_CURRENT_LIST_DIR}/*.h*)
file(GLOB_RECURSE all_sources ${CMAKE_CURRENT_LIST_DIR}/*.c*)

add_library(${PROJECT_LIB_NAME} STATIC ${all_headers};${all_sources})

tkl_initialize_library_target_defaults(${PROJECT_LIB_NAME} "32bit")

tkl_source_groups_from_dir_list("Header Files" FILES ${CMAKE_CURRENT_LIST_DIR} *.h*)
tkl_source_groups_from_dir_list("Source Files" FILES ${CMAKE_CURRENT_LIST_DIR} *.c*)

target_include_directories(${PROJECT_LIB_NAME}
  PUBLIC
    ${CMAKE_CURRENT_LIST_DIR}/core
)

# integration with QD

tkl_add_target_subdirectory(QD_ROOT qd ${CMAKE_BUILD_ROOT}/_3dparty/math/qd)

add_dependencies(${PROJECT_LIB_NAME} qd)

target_link_libraries(${PROJECT_LIB_NAME}
  PUBLIC
    qd
)

if (DEFINED ENABLE_QD_INTEGRATION)
  tkl_add_target_compile_definitions(${PROJECT_LIB_NAME} *
    PUBLIC
      ENABLE_QD_INTEGRATION=${ENABLE_QD_INTEGRATION}
  )
endif()
if (DEFINED ENABLE_QD_QD_INTEGRATION)
  tkl_add_target_compile_definitions(${PROJECT_LIB_NAME} *
    PUBLIC
      ENABLE_QD_QD_INTEGRATION=${ENABLE_QD_QD_INTEGRATION}
  )
endif()
if (DEFINED ENABLE_QD_DD_INTEGRATION)
  tkl_add_target_compile_definitions(${PROJECT_LIB_NAME} *
    PUBLIC
      ENABLE_QD_DD_INTEGRATION=${ENABLE_QD_DD_INTEGRATION}
  )
endif()
