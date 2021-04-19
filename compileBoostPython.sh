#!/bin/bash


# check for environment variable that specifies where to find Boost
if [[ ! "${BOOST_ROOT}" ]]
then
	echo '!!! error: ${BOOST_ROOT} environment variable is not set or empty. cannot compile BOOST libraries.'
	exit 1
fi


# check whether Boost directory exists and define log file location
BOOST_ROOT_DIR=$(readlink --canonicalize ${BOOST_ROOT})
LOG_FILE=${BOOST_ROOT_DIR}/compileBoostLibs.log
if [[ -d "${BOOST_ROOT_DIR}" ]]
then
	echo ">>> compiling Boost libraries in '${BOOST_ROOT_DIR}'"
	cd "${BOOST_ROOT_DIR}"
else
	echo "!!! error: '${BOOST_ROOT_DIR}' does not exist"
	exit 1
fi


# make sure Boost MPI module is included in build
JAM_CONFIG_FILE="${BOOST_ROOT_DIR}/tools/build/v2/user-config.jam"
if [[ ! -e ${JAM_CONFIG_FILE} ]]
then
	echo "!!! error: '${JAM_CONFIG_FILE}' does not exist"
	exit 1
fi
LINE='using mpi ;'
if ! grep --quiet --line-regexp "${LINE}" ${JAM_CONFIG_FILE}
then
	echo "${LINE}" >> ${JAM_CONFIG_FILE}
fi


echo "    view result in ${LOG_FILE}"
./bootstrap.sh --with-libraries=python --prefix=. &> ${LOG_FILE}
./bjam >> ${LOG_FILE} 2>&1


exit 0
