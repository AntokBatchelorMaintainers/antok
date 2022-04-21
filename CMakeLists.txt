# check if cmake has the required version
message(STATUS ">>> This is CMake version ${CMAKE_VERSION}.")
cmake_minimum_required(VERSION 3.11.0 FATAL_ERROR)
cmake_policy(VERSION 3.11)


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


# the c++ standard is given by the C++ sttandard with which ROOT was compiled
set(SUPPORTED_CXX_STANDARDS 11 14 17)
if(NOT ROOT_CPP_STD IN_LIST SUPPORTED_CXX_STANDARDS)
  message(FATAL_ERROR "The C++ standard '${ROOT_CPP_STD}' with which ROOT was compiled must be one of '${SUPPORTED_CXX_STANDARDS}'.")
endif()
set(CMAKE_CXX_STANDARD ${ROOT_CPP_STD})
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# set common compiler flags
# consider flags -pedantic -Wsuggest-attribute=pure -Wsuggest-attribute=noreturn -Wsuggest-attribute=const
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Werror")
# options switchable via command line
# enable option by running cmake .. -D <VAR>=<any value>
if(ARCH_NATIVE)
	# !Note! may cause trouble in heterogeneous environments
	# `gcc -march=native -Q --help=target` shows the compiler switches set by this option on your machine
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
endif()
# set additional flags for specific build types
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g")
# consider adding -DNODEBUG (disables all asserts) and -ffast-math
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELEASE} -g")
# report global build settings
foreach(_LANG "C" "CXX")
	message(STATUS "${_LANG} compiler ID is '${CMAKE_${_LANG}_COMPILER_ID} V${CMAKE_${_LANG}_COMPILER_VERSION}'.")
	message(STATUS "Using ${_LANG} compiler '${CMAKE_${_LANG}_COMPILER}'.")
	message(STATUS "Using general ${_LANG} compiler flags '${CMAKE_${_LANG}_FLAGS}'.")
endforeach()
unset(_LANG)
message(STATUS "Build type is '${CMAKE_BUILD_TYPE}'.")
message(STATUS "Using additional CXX compiler flags '${CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE}}' "
  "for build type ${CMAKE_BUILD_TYPE}.")
message(STATUS "Using additional linker flags '${CMAKE_CXX_LDFLAGS_${CMAKE_BUILD_TYPE}}' "
  "for build type ${CMAKE_BUILD_TYPE}.")
if(DEBUG_OUTPUT)
	include(CMakePrintSystemInformation)
endif()
message(STATUS "Includes defined via environment variables:")
message(STATUS "    CPATH (corresponds to -I...): '$ENV{CPATH}'")
foreach(_LANG "C" "CPLUS")
	message(STATUS "    ${_LANG}_INCLUDE_PATH (corresponds to -isystem ...): '$ENV{${_LANG}_INCLUDE_PATH}'")
endforeach()
unset(_LANG)


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


# setup Python
message(STATUS "")
message(STATUS ">>> Setting up Python.")
option(USE_PYTHON "En/Disables Python bindings" ON)
if(USE_PYTHON)
	find_package(PYTHON "${ROOT_PYTHON_MAJOR_VERSION}.${ROOT_PYTHON_MINOR_VERSION}")
	if(NOT PYTHONLIBS_FOUND)
		set(USE_PYTHON OFF)
		message(STATUS "Cannot find Python libraries. No Python bindings will be generated.")
	endif()
else()
	message(STATUS "'USE_PYTHON' is set to false.")
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
find_package(Boost 1.67.0 REQUIRED OPTIONAL_COMPONENTS "python${Python_VERSION_MAJOR}${Python_VERSION_MINOR}")
if(Boost_FOUND AND USE_PYTHON)
	message(STATUS "Using Boost include directory '${Boost_INCLUDE_DIRS}'.")
	message(STATUS "Using Boost library directory '${Boost_LIBRARY_DIRS}'.")
	message(STATUS "Using Boost libraries '${Boost_LIBRARIES}'.")
	if(NOT ${Boost_python${Python_VERSION_MAJOR}${Python_VERSION_MINOR}_FOUND})
		set(USE_PYTHON OFF)
		message(STATUS "Cannot find Boost.Python library in "
		  "'${Boost_LIBRARY_DIRS}'. Please consult INSTALL on how to build Boost.Python.")
	else()
		set(USE_PYTHON ON)
	endif()
endif()
if(USE_PYTHON)
	message(STATUS "Using Boost.Python libraries '${Boost_PYTHON_LIBS}'.")
	message(STATUS "Enabling compilation of Python components.")
	add_definitions(-DUSE_PYTHON)
	set(Python_LIBRARY_OUTPUT_PATH "${CMAKE_BINARY_DIR}/pyLib")
	file(MAKE_DIRECTORY ${Python_LIBRARY_OUTPUT_PATH})
	file(MAKE_DIRECTORY ${Python_LIBRARY_OUTPUT_PATH}/pytok)
	message(STATUS "Using python bindings output path '${Python_LIBRARY_OUTPUT_PATH}'.")
else()
	message(WARNING "Did not find Python or Boost.Python. No Python bindings will be generated. "
	  "Please consult INSTALL.")
endif()


# setup yaml-cpp
# environment variable $YAML_CPP is expected to point to yaml-cpp location
message(STATUS "")
message(STATUS ">>> Setting up yaml-cpp library.")
find_package(YamlCpp 0.6.0 REQUIRED)


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
