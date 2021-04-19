message(STATUS "")
message(STATUS ">>> Setting up 'pyInterface' directory.")

include_directories(
	${ANTOK_INCLUDE_DIR}
	${CMAKE_CURRENT_SOURCE_DIR}
	SYSTEM
	${Boost_INCLUDE_DIRS}
	${PYTHON_INCLUDE_DIRS}
	${ROOT_INCLUDE_DIR}
	${YamlCpp_INCLUDE_DIR}
	)

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
	"${Boost_PYTHON${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR}_LIBRARY}"
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

add_dependencies(libPytokLink pytok pytokModuleLinks)
