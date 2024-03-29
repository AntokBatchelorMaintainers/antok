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
#//      cmake module for finding Python 2.x installation
#//      replaces the official FindPYTHON that comes with Cmake which for several reasons is unusable
#//      requires executable of Python interpreter 'python' to be in PATH
#//
#//      following variables are defined:
#//      PYTHONINTERP_FOUND   - Was the Python executable found
#//      PYTHON_EXECUTABLE    - path to the Python interpreter
#//      PYTHON_VERSION       - Python version found e.g. 2.5.2
#//      PYTHON_VERSION_MAJOR - Python major version found e.g. 2
#//      PYTHON_VERSION_MINOR - Python minor version found e.g. 5
#//      PYTHON_VERSION_PATCH - Python patch version found e.g. 2
#//      PYTHONLIBS_FOUND     - have the Python libs been found
#//      PYTHON_LIBRARIES     - path to the Python library
#//      PYTHON_INCLUDE_DIRS  - path to includes
#//      PYTHONLIBS_VERSION   - version of the Python libs found
#//
#//      Example usage:
#//          find_package(PYTHON 2.7 REQUIRED)
#//
#//
#// Author List:
#//      Boris Grube          TUM            (original author)
#//
#//
#//-------------------------------------------------------------------------


set(PYTHON_FOUND        FALSE)
set(PYTHONINTERP_FOUND  FALSE)
set(PYTHONLIBS_FOUND    FALSE)
set(PYTHON_ERROR_REASON "")

set(PYTHON_EXECUTABLE    NOTFOUND)
set(PYTHON_VERSION       NOTFOUND)
set(PYTHON_VERSION_MAJOR NOTFOUND)
set(PYTHON_VERSION_MINOR NOTFOUND)
set(PYTHON_VERSION_PATCH NOTFOUND)
set(PYTHON_LIBRARIES     NOTFOUND)
set(PYTHON_INCLUDE_DIRS  NOTFOUND)
set(PYTHONLIBS_VERSION   NOTFOUND)


# find Python interpreter
set(_PYTHON_EXEC_NAMES "python" "python2")  # define search order for multiple python executable names
foreach(_PYTHON_EXEC_NAME ${_PYTHON_EXEC_NAMES})
	find_program(PYTHON_EXECUTABLE ${_PYTHON_EXEC_NAME})
	if(PYTHON_EXECUTABLE)
		break()
	endif()
endforeach()
if(NOT PYTHON_EXECUTABLE)
	set(PYTHON_ERROR_REASON "${PYTHON_ERROR_REASON} Cannot find executable of Python interpreter with name(s) '${_PYTHON_EXEC_NAMES}' in path. Make sure Python is setup correctly.")
else()
	set(PYTHONINTERP_FOUND TRUE)
endif()
unset(_PYTHON_EXEC_NAMES)
unset(_PYTHON_EXEC_NAME)


# get version number
if(PYTHONINTERP_FOUND)
	execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "import sys; print(sys.version_info[0])"
		OUTPUT_VARIABLE PYTHON_VERSION_MAJOR
		OUTPUT_STRIP_TRAILING_WHITESPACE)
	execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "import sys; print(sys.version_info[1])"
		OUTPUT_VARIABLE PYTHON_VERSION_MINOR
		OUTPUT_STRIP_TRAILING_WHITESPACE)
	execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "import sys; print(sys.version_info[2])"
		OUTPUT_VARIABLE PYTHON_VERSION_PATCH
		OUTPUT_STRIP_TRAILING_WHITESPACE)
	execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "import sys; print(sys.version_info[4] if sys.version_info[3] == 'candidate' else '')"
		OUTPUT_VARIABLE PYTHON_VERSION_RC
		OUTPUT_STRIP_TRAILING_WHITESPACE)
	if(PYTHON_VERSION_RC)
		set(PYTHON_VERSION "${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}.${PYTHON_VERSION_PATCH}.${PYTHON_VERSION_RC}")
	else()
		set(PYTHON_VERSION "${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}.${PYTHON_VERSION_PATCH}")
	endif()
endif()

# find the library
if(PYTHONINTERP_FOUND)
	set(_PYTHON_LIBRARY_DIRS)
	if(PYTHON_VERSION VERSION_LESS 3.2)
		# get name of shared library
		execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "from distutils import sysconfig; print('python' + sysconfig.get_config_var('VERSION'))"
			OUTPUT_VARIABLE _PYTHON_LIBRARY_NAME
			OUTPUT_STRIP_TRAILING_WHITESPACE)
		# get path of shared library
		execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "from distutils import sysconfig; print(sysconfig.get_config_var('LIBPL'))"
			OUTPUT_VARIABLE _PYTHON_LIBRARY_DIR
			OUTPUT_STRIP_TRAILING_WHITESPACE)
		list(APPEND _PYTHON_LIBRARY_DIRS ${_PYTHON_LIBRARY_DIR})
		execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "from distutils import sysconfig; print(sysconfig.get_config_var('LIBDIR'))"
			OUTPUT_VARIABLE _PYTHON_LIBRARY_DIR
			OUTPUT_STRIP_TRAILING_WHITESPACE)
		list(APPEND _PYTHON_LIBRARY_DIRS ${_PYTHON_LIBRARY_DIR})
	else()
		# get name of shared library
		execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "import sys; import sysconfig; print('python' + sysconfig.get_config_var('VERSION') + sys.abiflags)"
			OUTPUT_VARIABLE _PYTHON_LIBRARY_NAME
			OUTPUT_STRIP_TRAILING_WHITESPACE)
		# get path of shared library
		set(_PYTHON_LIBRARY_DIRS)
		execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "import sysconfig; print(sysconfig.get_config_var('LIBPL'))"
			OUTPUT_VARIABLE _PYTHON_LIBRARY_DIR
			OUTPUT_STRIP_TRAILING_WHITESPACE)
		list(APPEND _PYTHON_LIBRARY_DIRS ${_PYTHON_LIBRARY_DIR})
		execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "import sysconfig; print(sysconfig.get_config_var('LIBDIR'))"
			OUTPUT_VARIABLE _PYTHON_LIBRARY_DIR
			OUTPUT_STRIP_TRAILING_WHITESPACE)
		list(APPEND _PYTHON_LIBRARY_DIRS ${_PYTHON_LIBRARY_DIR})
	endif()
	find_library(PYTHON_LIBRARIES
		NAMES ${_PYTHON_LIBRARY_NAME}
		PATHS ${_PYTHON_LIBRARY_DIRS}
		NO_DEFAULT_PATH)
	if(NOT PYTHON_LIBRARIES)
		set(PYTHONLIBS_FOUND FALSE)
		set(PYTHONINTERP_FOUND FALSE)
		set(PYTHON_ERROR_REASON "${PYTHON_ERROR_REASON} Cannot find Python shared library '${_PYTHON_LIBRARY_NAME}' in '${_PYTHON_LIBRARY_DIRS}'. Make sure Python is setup correctly.")
	else()
		set(PYTHONLIBS_FOUND TRUE)
	endif()
	unset(_PYTHON_LIBRARY_NAME)
	unset(_PYTHON_LIBRARY_DIR)
	unset(_PYTHON_LIBRARY_DIRS)
endif()


# find the include directory and get library version
if(PYTHONINTERP_FOUND)
	if(PYTHON_VERSION VERSION_LESS 3.2)
		execute_process(
			COMMAND ${PYTHON_EXECUTABLE} -c "from distutils import sysconfig; print('{0};{1}'.format(sysconfig.get_python_inc(), sysconfig.get_python_inc(plat_specific=True)))"
			OUTPUT_VARIABLE PYTHON_INCLUDE_DIRS
			OUTPUT_STRIP_TRAILING_WHITESPACE)
	else()
		execute_process(
			COMMAND ${PYTHON_EXECUTABLE} -c "import sysconfig; print('{};{}'.format(sysconfig.get_path('include', 'posix_prefix'), sysconfig.get_path('platinclude', 'posix_prefix')))"
			OUTPUT_VARIABLE PYTHON_INCLUDE_DIRS
			OUTPUT_STRIP_TRAILING_WHITESPACE)
	endif()
	list(REMOVE_DUPLICATES PYTHON_INCLUDE_DIRS)
	# filter out non-existing directories
	set(_NONEXISTING_PYTHON_INCLUDE_DIRS "")
	foreach(_PYTHON_INCLUDE_DIR ${PYTHON_INCLUDE_DIRS})
		if(NOT EXISTS "${_PYTHON_INCLUDE_DIR}")
			list(APPEND _NONEXISTING_PYTHON_INCLUDE_DIRS "${_PYTHON_INCLUDE_DIR}")
		endif()
	endforeach()
	if(_NONEXISTING_PYTHON_INCLUDE_DIRS)
		list(REMOVE_ITEM PYTHON_INCLUDE_DIRS ${_NONEXISTING_PYTHON_INCLUDE_DIRS})
	endif()
	if(NOT PYTHON_INCLUDE_DIRS)
		set(PYTHONLIBS_FOUND FALSE)
		set(PYTHON_ERROR_REASON "${PYTHON_ERROR_REASON} None of the include directories '${_NONEXISTING_PYTHON_INCLUDE_DIRS}' exist.")
	endif()
	unset(_NONEXISTING_PYTHON_INCLUDE_DIRS)
	# select directories containing patchlevel.h with same version as python interpreter
	set(_NO_PATCHLEVELH_PYTHON_INCLUDE_DIRS "")
	foreach(_PYTHON_INCLUDE_DIR ${PYTHON_INCLUDE_DIRS})
		# get library version from patchlevel.h
		parse_version_from_file("${_PYTHON_INCLUDE_DIR}/patchlevel.h" "#define[ \t]+PY_VERSION[ \t]+" PYTHONLIBS_VERSION)
		# cut release candidate number, i.e. 4th version number, from PYTHONLIBS_VERSION
		STRING(REGEX REPLACE "^([0-9]+)\\.([0-9]+)\\.([0-9]+)(\\.[0-9]+)*$" "\\1.\\2.\\3" PYTHONLIBS_VERSION ${PYTHONLIBS_VERSION})
		if(NOT PYTHONLIBS_VERSION VERSION_EQUAL PYTHON_VERSION)
			list(APPEND _NO_PATCHLEVELH_PYTHON_INCLUDE_DIRS "${_PYTHON_INCLUDE_DIR}")
		endif()
	endforeach()
	unset(_PYTHON_INCLUDE_DIR)
	if(_NO_PATCHLEVELH_PYTHON_INCLUDE_DIRS)
		list(REMOVE_ITEM PYTHON_INCLUDE_DIRS ${_NO_PATCHLEVELH_PYTHON_INCLUDE_DIRS})
	endif()
	if(NOT PYTHON_INCLUDE_DIRS)
		set(PYTHONLIBS_FOUND FALSE)
		set(PYTHON_ERROR_REASON "${PYTHON_ERROR_REASON} None of the include directories '${_NO_PATCHLEVELH_PYTHON_INCLUDE_DIRS}' contains a patchlevel.h of the same version ${PYTHON_VERSION} as the Python interpreter.")
	endif()
	unset(_NO_PATCHLEVELH_PYTHON_INCLUDE_DIRS)
	if(PYTHON_INCLUDE_DIRS)
		set(PYTHONLIBS_FOUND TRUE)
	endif()
endif()


if(PYTHON_ERROR_REASON AND NOT PYTHON_FIND_QUIETLY)
	message(STATUS "Problems while finding the requested PYTHON installation:${PYTHON_ERROR_REASON}")
endif()
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PYTHON
	FOUND_VAR PYTHON_FOUND
	REQUIRED_VARS PYTHON_EXECUTABLE PYTHONINTERP_FOUND PYTHONLIBS_FOUND PYTHON_VERSION PYTHON_VERSION_MAJOR PYTHON_VERSION_MINOR PYTHON_VERSION_PATCH PYTHON_LIBRARIES PYTHON_INCLUDE_DIRS PYTHONLIBS_VERSION
	VERSION_VAR PYTHON_VERSION)
if(NOT PYTHON_FOUND)
	set(PYTHONINTERP_FOUND FALSE)
	set(PYTHONLIBS_FOUND   FALSE)
endif()
# additional reporting
if(PYTHONLIBS_FOUND AND NOT PYTHON_FIND_QUIETLY)
	message(STATUS "Found Python libraries version ${PYTHONLIBS_VERSION}.")
	message(STATUS "Using Python libraries '${PYTHON_LIBRARIES}'.")
	message(STATUS "Using Python include directory '${PYTHON_INCLUDE_DIRS}'.")
endif()


# hide variables from normal GUI
mark_as_advanced(
	PYTHON_EXECUTABLE
	PYTHON_VERSION
	PYTHON_VERSION_MAJOR
	PYTHON_VERSION_MINOR
	PYTHON_VERSION_PATCH
	PYTHON_LIBRARIES
	PYTHON_INCLUDE_DIRS
	PYTHONLIBS_VERSION
	)


if(NOT PYTHON_FOUND)
	unset(PYTHON_EXECUTABLE)
	unset(PYTHON_VERSION)
	unset(PYTHON_VERSION_MAJOR)
	unset(PYTHON_VERSION_MINOR)
	unset(PYTHON_VERSION_PATCH)
	unset(PYTHON_LIBRARIES)
	unset(PYTHON_INCLUDE_DIRS)
	unset(PYTHONLIBS_VERSION)
endif()
