#///////////////////////////////////////////////////////////////////////////
#//
#//    Copyright 2010
#//
#//    This file is part of rootpwa
#//
#//    rootpwa is free software: you can redistribute it and/or modify
#//    it under the terms of the GNU General Public License as published by
#//    the Free Software Foundation, either version 3 of the License, or
#//    (at your option) any later version.
#//
#//    rootpwa is distributed in the hope that it will be useful,
#//    but WITHOUT ANY WARRANTY; without even the implied warranty of
#//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#//    GNU General Public License for more details.
#//
#//    You should have received a copy of the GNU General Public License
#//    along with rootpwa.  If not, see <http://www.gnu.org/licenses/>.
#//
#///////////////////////////////////////////////////////////////////////////
#//-------------------------------------------------------------------------
#//
#// Description:
#//      build file for amplitude library
#//
#//
#// Author List:
#//      Boris Grube          TUM            (original author)
#//
#//
#//-------------------------------------------------------------------------


message(STATUS "")
message(STATUS ">>> Setting up 'pyInterface' directory.")


# set include directories
set(INCLUDE_DIRECTORIES
	${PYTHON_INCLUDE_DIRS}
	SYSTEM ${Boost_INCLUDE_DIRS}
	${ROOT_INCLUDE_DIR}
	${CMAKE_CURRENT_SOURCE_DIR}
	${ANTOK_INCLUDE_DIR}
	${YamlCpp_INCLUDE_DIR}
	)
include_directories(${INCLUDE_DIRECTORIES})


# source files that are compiled into library
set(SOURCES
	basic_calcs_py.cxx
	constants_py.cxx
	initializer_py.cxx
	pytok.cxx
	rootConverters_py.cxx
)

# library
set(THIS_LIB "pytok")
make_shared_library(
	"${THIS_LIB}"
	"${SOURCES}"
	"${PYTHON_LIBRARIES}"
	"${Boost_PYTHON_LIBRARY}"
	"${ROOT_LIBS}"
	"${ANTOK_LIB}"
	# no idea why this is needed explicitly here and not in rootpwa
	"${ROOT_LIBRARY_DIR}/libPyROOT.so"
	)

# link the shared library with the python bindings
add_custom_target(
	libPytokLink ALL
	COMMAND ${CMAKE_COMMAND} -E create_symlink ${LIBRARY_OUTPUT_PATH}/${CMAKE_SHARED_LIBRARY_PREFIX}pytok${CMAKE_SHARED_LIBRARY_SUFFIX} ${PYTHON_LIBRARY_OUTPUT_PATH}/pytok/libPytok.so
)

# produce python byte-code and move it to the build directory
add_custom_target(
	pytokModuleLinks ALL
	COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_SOURCE_DIR}/pytok/helperScripts/python_compiler.py ${CMAKE_SOURCE_DIR}/pytok/module ${PYTHON_LIBRARY_OUTPUT_PATH}/pytok
)

add_dependencies(libPytokLink pytok pytokModuleLinks )
