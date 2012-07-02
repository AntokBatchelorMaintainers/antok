
set(INCLUDE_DIRECTORIES
	${ROOT_INCLUDE_DIR}
	${CMAKE_CURRENT_SOURCE_DIR}
   )
include_directories(${INCLUDE_DIRECTORIES})

set(SOURCES
	basic_calcs.cxx
   )

make_shared_library("bannt" "${SOURCES}" "${ROOT_LIBS}")
