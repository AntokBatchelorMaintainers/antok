
set(INCLUDE_DIRECTORIES
	${CMAKE_CURRENT_SOURCE_DIR}
	${ROOT_INCLUDE_DIR}
	${YAMLCPP_INCLUDE_DIR}
   )
include_directories(${INCLUDE_DIRECTORIES})

set(SOURCES
	basic_calcs.cxx
	constants.cxx
	cutter.cxx
	data.cxx
	event.cxx
	generators_functions.cxx
	initializer.cxx
	object_manager.cxx
	plot.cxx
	plotter.cxx
   )

make_shared_library("antok" "${SOURCES}" "${ROOT_LIBS}" "${YAMLCPP_LIBRARY}")

