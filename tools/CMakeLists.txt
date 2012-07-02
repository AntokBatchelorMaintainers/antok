
set(INCLUDE_DIRECTORIES
	${ROOT_INCLUDE_DIR}
	${CMAKE_CURRENT_SOURCE_DIR}
	${BANNT_INCLUDE_DIR}
   )
include_directories(${INCLUDE_DIRECTORIES})

make_executable(treereader			treereader.cxx		${BANNT_LIB})
make_executable(gen_kin				gen_kin.cxx			${BANNT_LIB})
make_executable(root_to_txt			root_to_txt.cxx		${BANNT_LIB})
