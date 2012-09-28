
set(INCLUDE_DIRECTORIES
	${ROOT_INCLUDE_DIR}
	${CMAKE_CURRENT_SOURCE_DIR}
	${ANTOK_INCLUDE_DIR}
   )
include_directories(${INCLUDE_DIRECTORIES})

make_executable(treereader			treereader.cxx		${ANTOK_LIB})
make_executable(gen_kin				gen_kin.cxx			${ANTOK_LIB})
make_executable(root_to_root		root_to_root.cxx	${ANTOK_LIB})

