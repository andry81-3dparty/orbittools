cmake_minimum_required(VERSION 3.0)

include(Common)

project("orbitTools")

set(PROJECT_LIB_NAME orbitTools)

discover_variable(QD_ROOT PATH "QD math library")
discover_variable(ORBIT_TOOLS_ENABLE_QD_QD_INTEGRATION STRING "orbitTools library qd_real replacement for the double")
discover_variable(ORBIT_TOOLS_ENABLE_QD_DD_INTEGRATION STRING "orbitTools library dd_real replacement for the double")


file(GLOB_RECURSE all_headers ${CMAKE_CURRENT_LIST_DIR}/*.h*)
file(GLOB_RECURSE all_sources ${CMAKE_CURRENT_LIST_DIR}/*.cpp)

add_library(${PROJECT_LIB_NAME} STATIC ${all_headers};${all_sources})

source_groups_from_dir_list("Header Files" FILES ${CMAKE_CURRENT_LIST_DIR} *.h*)
source_groups_from_dir_list("Source Files" FILES ${CMAKE_CURRENT_LIST_DIR} *.cpp)

target_include_directories(${PROJECT_LIB_NAME}
  PUBLIC
    ${CMAKE_CURRENT_LIST_DIR}/core
)

if(GCC)
  add_target_compile_definitions(${PROJECT_LIB_NAME} DEBUG
    PUBLIC
      _DEBUG
  )

  add_target_compile_properties(${PROJECT_LIB_NAME} DEBUG
    -g
  )

  add_target_compile_properties(${PROJECT_LIB_NAME} *
    -m32        # compile 32 bit target on 64 bit host
  )

  add_target_link_properties(${PROJECT_LIB_NAME} *
    -m32        # link 32 bit target on 64 bit host
  )
endif()

# integration with QD

add_target_subdirectory(qd ${QD_ROOT} ${CMAKE_BUILD_ROOT}/_3dparty/math/qd)

add_dependencies(${PROJECT_LIB_NAME} qd)

target_link_libraries(${PROJECT_LIB_NAME}
  PUBLIC
    qd
)

if (ORBIT_TOOLS_ENABLE_QD_QD_INTEGRATION)
  add_target_compile_definitions(${PROJECT_LIB_NAME} *
    PUBLIC
      ORBIT_TOOLS_ENABLE_QD_QD_INTEGRATION
  )
endif()
if (ORBIT_TOOLS_ENABLE_QD_DD_INTEGRATION)
  add_target_compile_definitions(${PROJECT_LIB_NAME} *
    PUBLIC
      ORBIT_TOOLS_ENABLE_QD_DD_INTEGRATION
  )
endif()
