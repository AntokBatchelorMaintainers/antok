
set(INCLUDE_DIRECTORIES
	${ROOT_INCLUDE_DIR}
	${CMAKE_CURRENT_SOURCE_DIR}
   )
include_directories(${INCLUDE_DIRECTORIES})

set(SOURCES
	basic_calcs.cxx
	cutter.cxx
	event.cxx
	initializer.cxx
	plot.cxx
	plotter.cxx
   )

make_shared_library("antok" "${SOURCES}" "${ROOT_LIBS}")

