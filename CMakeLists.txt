# check if cmake has the required version
message(STATUS ">>> This is CMake version ${CMAKE_VERSION}.")
cmake_minimum_required(VERSION 2.8.0 FATAL_ERROR)
cmake_policy(SET CMP0042 NEW) # can be removed once required CMake version is
                              # increased to or beyond version 3.0


# set verbosity
set(CMAKE_VERBOSE_MAKEFILE 0)


# define project
message(STATUS "")
message(STATUS ">>> Setting up project 'ANTOK'.")
project(ANTOK)


# set path, where to look first for cmake modules, before ${CMAKE_ROOT}/Modules/ is checked
message(STATUS "")
message(STATUS ">>> Setting up Cmake modules.")
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmakeModules")
message(STATUS "Using cmake module path '${CMAKE_MODULE_PATH}'.")
# load cmake macros
include(CommonMacros)


enforce_out_of_source_build()


# require C++11 compatible compiler and set compiler flags
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX11)
	set(CMAKE_CXX_FLAGS "-g -Wall -Werror -std=c++11")
elseif(COMPILER_SUPPORTS_CXX0X)
	set(CMAKE_CXX_FLAGS "-g -Wall -Werror")
else()
	message(FATAL_ERROR "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()


# define buid type 'ASAN' that enables compilation with AddressSanitizer
# to enable, add -DCMAKE_BUILD_TYPE=ASAN to cmake call
# see https://github.com/google/sanitizers/wiki/AddressSanitizer
# and http://www.stablecoder.ca/2018/02/01/analyzer-build-types.html
set(CMAKE_CXX_FLAGS_ASAN
	"${CMAKE_CXX_FLAGS} -fsanitize=address -fno-omit-frame-pointer -fno-optimize-sibling-calls -fsanitize-address-use-after-scope"
	CACHE STRING "Flags used by the C++ compiler for 'ASAN' build types that enable the AddressSanitizer."
	FORCE)
set(CMAKE_LINKER_FLAGS_ASAN
	"${CMAKE_LINKER_FLAGS} -fsanitize=address -fno-omit-frame-pointer -fno-optimize-sibling-calls -fsanitize-address-use-after-scope"
	CACHE STRING "Flags used by the C++ linker for 'ASAN' build types that enable the AddressSanitizer."
	FORCE)


# redirect output files
message(STATUS "")
message(STATUS ">>> Setting up output paths.")
set(LIBRARY_OUTPUT_PATH "${CMAKE_BINARY_DIR}/lib")
message(STATUS "Using library output path '${LIBRARY_OUTPUT_PATH}'.")
set(EXECUTABLE_OUTPUT_PATH "${CMAKE_BINARY_DIR}/bin")
message(STATUS "Using executable output path '${EXECUTABLE_OUTPUT_PATH}'.")


# setup ROOT includes and libraries
message(STATUS "")
message(STATUS ">>> Setting up ROOT.")
find_package(ROOT 5.34.04 REQUIRED)
# ROOT version 6.00.X and 6.03.X have to be excluded, as they have no FindFixBin implementation
if(ROOT_MAJOR_VERSION EQUAL 6)
	if(    (ROOT_MINOR_VERSION EQUAL 00)
	    OR (ROOT_MINOR_VERSION EQUAL 02 AND ROOT_VERSION VERSION_LESS 6.02.05)
	    OR (ROOT_MINOR_VERSION EQUAL 03 AND ROOT_VERSION VERSION_LESS 6.03.04) )
		message(FATAL_ERROR "ROOT version ${ROOT_VERSION} cannot be used!")
	endif()
endif()


# setup yaml-cpp
# environment variable $YAML_CPP is expected to point to yaml-cpp location
message(STATUS "")
message(STATUS ">>> Setting up yaml-cpp library.")
find_package(YamlCpp 0.6.0 REQUIRED)


# setup Python
message(STATUS "")
message(STATUS ">>> Setting up Python.")
option(USE_PYTHON "En/Disables Python bindings" OFF)
find_package(Python 2.6)
if(NOT PYTHONLIBS_FOUND)
	set(USE_PYTHON OFF)
	message(STATUS "Cannot find Python libraries. No Python bindings will be generated.")
else()
	set(USE_PYTHON ON)
endif()


# setup Boost
# environment variable $BOOST_ROOT is expected to point to non-standard locations
message(STATUS "")
message(STATUS ">>> Setting up Boost library.")
set(Boost_USE_STATIC_LIBS      OFF)
set(Boost_USE_MULTITHREADED    ON)
set(Boost_USE_STATIC_RUNTIME   OFF)
set(Boost_DETAILED_FAILURE_MSG ON)
# set(Boost_DEBUG 1)
if(     ("$ENV{BOOST_ROOT}"       STREQUAL "")
    AND ("$ENV{BOOSTROOT}"        STREQUAL "")
    AND ("$ENV{BOOST_INCLUDEDIR}" STREQUAL "")
    AND ("$ENV{BOOST_LIBRARYDIR}" STREQUAL ""))
	set(Boost_NO_SYSTEM_PATHS OFF)
else()
	set(Boost_NO_SYSTEM_PATHS ON)
	set(Boost_NO_BOOST_CMAKE  ON)  # required starting from Boost 1.70.0 on
endif()
find_package(Boost 1.67.0 REQUIRED OPTIONAL_COMPONENTS "python${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR}")
if(Boost_FOUND)
	message(STATUS "Using Boost include directory '${Boost_INCLUDE_DIRS}'.")
	message(STATUS "Using Boost library directory '${Boost_LIBRARY_DIRS}'.")
	message(STATUS "Using Boost libraries '${Boost_LIBRARIES}'.")
	if(NOT ${Boost_python${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR}_FOUND})
		set(USE_PYTHON OFF)
		message(STATUS "Cannot find Boost.Python library in "
		  "'${Boost_LIBRARY_DIRS}'. Please consult INSTALL on how to build Boost.Python.")
	endif()
endif()
if(USE_PYTHON)
	message(STATUS "Using Boost.Python libraries '${Boost_PYTHON_LIBS}'.")
	message(STATUS "Enabling compilation of Python components.")
	add_definitions(-DUSE_PYTHON)
	set(PYTHON_LIBRARY_OUTPUT_PATH "${CMAKE_BINARY_DIR}/pyLib")
	file(MAKE_DIRECTORY ${PYTHON_LIBRARY_OUTPUT_PATH})
	file(MAKE_DIRECTORY ${PYTHON_LIBRARY_OUTPUT_PATH}/pytok)
	message(STATUS "Using python bindings output path '${PYTHON_LIBRARY_OUTPUT_PATH}'.")
else()
	message(WARNING "Did not find Python or Boost.Python. No Python bindings will be generated. "
	  "Please consult INSTALL.")
endif()


set(ANTOK_LIB "antok")
if(USE_PYTHON)
	set(PYINTERFACE_LIBS "pytok")
endif()
set(ANTOK_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/src")
set(USER_DIR "${CMAKE_SOURCE_DIR}/user")
aux_source_directory("${USER_DIR}" USER_SOURCES)
add_subdirectory("src")
add_subdirectory("tools")
if(USE_PYTHON)
	add_subdirectory("pytok/bindings")
endif()
