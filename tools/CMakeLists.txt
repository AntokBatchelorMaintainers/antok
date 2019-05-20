
set(INCLUDE_DIRECTORIES
	${CMAKE_CURRENT_SOURCE_DIR}
	SYSTEM ${Boost_INCLUDE_DIRS}
	${ANTOK_INCLUDE_DIR}
	${ROOT_INCLUDE_DIR}
	${YamlCpp_INCLUDE_DIR}
	${USER_DIR}
   )
include_directories(${INCLUDE_DIRECTORIES})

make_executable(treereader			treereader.cxx			${ANTOK_LIB})
make_executable(gen_kin				gen_kin.cxx				${ANTOK_LIB})
make_executable(root_to_root		root_to_root.cxx		${ANTOK_LIB})
make_executable(generateBeamfile	generateBeamfile.cxx	${ANTOK_LIB})
make_executable(fitRPDResolutions   fitRPDResolutions.cxx   ${ANTOK_LIB})
file( COPY runParallel.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/runParallel )
add_custom_target( runParallel ALL
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/runParallel.py ${EXECUTABLE_OUTPUT_PATH}/runParallel
    )
file( COPY run.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/run)
add_custom_target( run ALL
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/run.py ${EXECUTABLE_OUTPUT_PATH}/run
    )
