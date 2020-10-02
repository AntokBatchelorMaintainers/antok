message(STATUS "")
message(STATUS ">>> Setting up 'src' directory.")

include_directories(
	${CMAKE_CURRENT_SOURCE_DIR}
	${USER_DIR}
	SYSTEM
	${Boost_INCLUDE_DIRS}
	${ROOT_INCLUDE_DIR}
	${YamlCpp_INCLUDE_DIR}
	)

file(GLOB SOURCES *.cxx)
set(SOURCES ${SOURCES} ${USER_SOURCES})

make_shared_library("antok" "${SOURCES}" "${ROOT_LIBS}" "${YamlCpp_LIBS}")
