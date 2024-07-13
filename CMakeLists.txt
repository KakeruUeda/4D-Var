cmake_minimum_required(VERSION 3.16)

set(CMAKE_CXX_COMPILER mpic++)
set(base_compile_options ${base_compile_options} "-O2 -std=c++17")
string(REPLACE " " " " compile_options ${base_compile_options})

project(DataAssimilation C CXX)
set(CMAKE_CXX_STANDARD 11)

set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_MODULE_PATH  ${PROJECT_SOURCE_DIR}/cmake)

message(" ")
message("${PROJECT_NAME}")
message( STATUS "Destination PATH      : "      ${CMAKE_INSTALL_PREFIX})
message( STATUS "CMAKE_CXX_COMPILER    : "      ${CMAKE_CXX_COMPILER})
message( STATUS "COMPILE_OPTIONS       : "      ${compile_options})
message(" ")

# TextParser
find_package(TP REQUIRED)
include_directories(${TP_INCLUDE_DIR})

# Eigen
find_package(Eigen REQUIRED)
include_directories(${EIGEN_INCLUDE_DIR})

include_directories(
  /usr/local/petsc/include
)

link_directories(
  /usr/local/petsc/lib
  /opt/homebrew/opt/libomp/lib
)

add_subdirectory(src)
