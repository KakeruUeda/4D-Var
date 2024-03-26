cmake_minimum_required(VERSION 3.16)

set(CMAKE_CXX_COMPILER "/opt/petsc/bin/mpic++")
set(base_compile_options ${base_compile_options} "-O2")
string(REPLACE " " " " compile_options ${base_compile_options})

project(DataAssimilation C CXX)
set(CMAKE_CXX_STANDARD 11)

set(CMAKE_BUILD_TYPE Release)
set(CMAKE_MODULE_PATH  ${PROJECT_SOURCE_DIR}/cmake)

message(" ")
message("${PROJECT_NAME}")
message( STATUS "Destination PATH      : "      ${CMAKE_INSTALL_PREFIX})
message( STATUS "CMAKE_CXX_COMPILER    : "      ${CMAKE_CXX_COMPILER})
message( STATUS "COMPILE_OPTIONS       : "      ${compile_options})
message(" ")

include_directories(
  /opt/homebrew/Cellar/metis/5.1.0/include
  /opt/homebrew/include/eigen3
  /opt/TextParser/include
  /opt/petsc/include
)

link_directories(
  /opt/homebrew/Cellar/metis/5.1.0/lib
  /opt/TextParser/lib
  /opt/petsc/lib
  /opt/homebrew/opt/libomp/lib
)

add_subdirectory(src)
add_subdirectory(test)
