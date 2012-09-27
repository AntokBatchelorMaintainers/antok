# Locate yaml-cpp
#
# This module defines
#  YAMLCPP_FOUND, if false, do not try to link to yaml-cpp
#  YAMLCPP_LIBRARY, where to find yaml-cpp
#  YAMLCPP_INCLUDE_DIR, where to find yaml.h
#
# By default, the dynamic libraries of yaml-cpp will be found. To find the static ones instead,
# you must set the YAMLCPP_STATIC_LIBRARY variable to TRUE before calling find_package(YAMLCPP ...).
#
# If yaml-cpp is not installed in a standard path, you can use the YAMLCPP_DIR CMake variable
# to tell CMake where yaml-cpp is.

# attempt to find static library first if this is set
if(YAMLCPP_STATIC_LIBRARY)
    set(YAMLCPP_STATIC libyaml-cpp.a)
endif()

set(YAML_CPP $ENV{YAML_CPP})

# find the yaml-cpp include directory
find_path(YAMLCPP_INCLUDE_DIR yaml-cpp/yaml.h
          PATH_SUFFIXES include
          HINTS
          ~/Library/Frameworks/yaml-cpp/include/
          /sw/yaml-cpp/         # Fink
          /opt/local/yaml-cpp/  # DarwinPorts
          /opt/csw/yaml-cpp/    # Blastwave
          /opt/yaml-cpp/
          ${YAMLCPP_DIR}/include/
		  ${YAML_CPP}/include)

# find the yaml-cpp library
find_library(YAMLCPP_LIBRARY
             NAMES ${YAMLCPP_STATIC} yaml-cpp
             PATH_SUFFIXES lib64 lib
             HINTS  ~/Library/Frameworks
                    /Library/Frameworks
                    /sw
                    /opt/local
                    /opt/csw
                    /opt
                    ${YAMLCPP_DIR}/build
                    ${YAML_CPP}/build)

set(YamlCpp_FOUND True)
if (NOT YAMLCPP_INCLUDE_DIR)
	set(YAMLCPP_ERROR_REASON "yaml-cpp include directory not found.")
	set(YamlCpp_FOUND False)
endif()
if(NOT YAMLCPP_LIBRARY)
	set(YAMLCPP_ERROR_REASON "yaml-cpp library not found.")
	set(YamlCpp_FOUND False)
endif()

if(${YamlCpp_FOUND})
	string(REPLACE "/include" "" YAMLCPP ${YAMLCPP_INCLUDE_DIR})
endif()

if(YamlCpp_FOUND)
	message(STATUS "Found yaml in ${YAMLCPP}.")
else()
	if(YamlCpp_FIND_REQUIRED)
		message(FATAL_ERROR "Unable to find yaml-cpp: ${YAMLCPP_ERROR_REASON}")
	else()
		if(NOT YamlCpp_FIND_QUIETLY)
			message(STATUS "Unable to find yaml-cpp.")
		endif()
	endif()
endif()

mark_as_advanced(YAMLCPP_INCLUDE_DIR YAMLCPP_LIBRARY)

