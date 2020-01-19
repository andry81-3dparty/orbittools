cmake_minimum_required(VERSION 3.0)

include(tacklelib/Project)

project("orbitTools")

set(PROJECT_LIB_NAME orbitTools)

tkl_check_var(OPTIONAL PATH QD_ROOT)


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
