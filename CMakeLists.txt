
cmake_minimum_required(VERSION 2.8.0 FATAL_ERROR)

set(CMAKE_VERBOSE_MAKEFILE 0)

project(BANNT)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmakeModules")
message(STATUS "Using cmake module path '${CMAKE_MODULE_PATH}'")
include(CommonMacros)

enforce_out_of_source_build()

#C++ Flags
set(CMAKE_CXX_FLAGS "-g -Wall")

# redirect output files
set(LIBRARY_OUTPUT_PATH "${CMAKE_BINARY_DIR}/lib")
message(STATUS "Using library output path '${LIBRARY_OUTPUT_PATH}'")
set(EXECUTABLE_OUTPUT_PATH "${CMAKE_BINARY_DIR}/bin")
message(STATUS "Using executable output path '${EXECUTABLE_OUTPUT_PATH}'")

# setup ROOT includes and libraries
find_package(ROOT 5.26 REQUIRED)

# find yaml
find_package(YamlCpp 0.3.0 REQUIRED)

set(BANNT_LIB "bannt")
set(BANNT_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/src)
add_subdirectory(src)
add_subdirectory(tools)

