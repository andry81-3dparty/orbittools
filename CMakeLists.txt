cmake_minimum_required(VERSION 3.0)

include(tacklelib/Project)

project("orbitTools")

set(PROJECT_LIB_NAME orbitTools)

tkl_check_var(REQUIRED PATH   ORBITTOOLS_BRANCH_NAME) # ex: `1_4_0_20170422`

set(PROJECT_EXTERNALS_DIR "_externals")
set(ORBITTOOLS_EXTERNAL_SOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}/${PROJECT_EXTERNALS_DIR}/${ORBITTOOLS_BRANCH_NAME}")

tkl_check_var(OPTIONAL PATH QD_ROOT)
tkl_check_var(OPTIONAL . ENABLE_QD_INTEGRATION)
tkl_check_var(OPTIONAL . ENABLE_QD_QD_INTEGRATION)
tkl_check_var(OPTIONAL . ENABLE_QD_DD_INTEGRATION)


file(GLOB_RECURSE all_headers ${ORBITTOOLS_EXTERNAL_SOURCE_DIR}/*.h*)
file(GLOB_RECURSE all_sources ${ORBITTOOLS_EXTERNAL_SOURCE_DIR}/*.c*)

add_library(${PROJECT_LIB_NAME} STATIC ${all_headers};${all_sources})

tkl_initialize_library_target_defaults(${PROJECT_LIB_NAME} "anybit")

tkl_source_groups_from_dir_list("Header Files" FILES ${ORBITTOOLS_EXTERNAL_SOURCE_DIR} *.h*)
tkl_source_groups_from_dir_list("Source Files" FILES ${ORBITTOOLS_EXTERNAL_SOURCE_DIR} *.c*)

target_include_directories(${PROJECT_LIB_NAME}
  PUBLIC
    ${ORBITTOOLS_EXTERNAL_SOURCE_DIR}/core
)

if (DEFINED QD_ROOT)
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
endif()