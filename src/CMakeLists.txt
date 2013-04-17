
set(INCLUDE_DIRECTORIES
	${CMAKE_CURRENT_SOURCE_DIR}
	${Boost_INCLUDE_DIRS}
	${ROOT_INCLUDE_DIR}
	${YAMLCPP_INCLUDE_DIR}
   )
include_directories(${INCLUDE_DIRECTORIES})

set(SOURCES
	basic_calcs.cxx
	beamfile_generator_helpers.cxx
	beamfile_generator_5dBin.cxx
	constants.cxx
	cutter.cxx
	data.cxx
	event.cxx
	generators_cuts.cxx
	generators_functions.cxx
	generators_plots.cxx
	initializer.cxx
	object_manager.cxx
	plotter.cxx
   )

make_shared_library("antok" "${SOURCES}" "${ROOT_LIBS}" "${YAMLCPP_LIBRARY}")

