
set(INCLUDE_DIRECTORIES
	${CMAKE_CURRENT_SOURCE_DIR}
	SYSTEM ${Boost_INCLUDE_DIRS}
	${ROOT_INCLUDE_DIR}
	${YamlCpp_INCLUDE_DIR}
	${USER_DIR}
   )
include_directories(${INCLUDE_DIRECTORIES})

file(GLOB SOURCES *.cxx ${USER_SOURCES})

make_shared_library("antok" "${SOURCES}" "${ROOT_LIBS}" "${YamlCpp_LIBS}")

