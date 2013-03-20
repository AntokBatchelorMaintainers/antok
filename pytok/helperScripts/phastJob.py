#!/usr/bin/env python

import datetime
import os
import sys

import pytok


if __name__ == "__main__":

	level = "Debug"
	logger = pytok.Logger(level)

	startingTime = datetime.datetime.now()

	logger.info("Job starting")

	if len(sys.argv) != 2:
		logger.fatal("Got " + str(len(sys.argv)) + " arguments instead of 2. Aborting...")
		pytok.exit(5, logger)

	argumentsFileName = sys.argv[1]

	if not os.path.isfile(argumentsFileName):
		logger.fatal("Could not find arguments file '" + argumentsFileName + "'. Aborting...")
		pytok.exit(5, logger)

	logger.info("Reading arguments file '" + argumentsFileName + "'.")
	with open(argumentsFileName, 'r') as argumentsFile:
		arguments = [x[:-1] for x in argumentsFile.readlines()]

	if len(arguments) < 6:
		logger.fatal("Arguments file contained " + str(len(arguments)) + " lines instead of at least 6.")
		pytok.exit(5, logger)

	try:
		workingDirectory = os.environ['TMPDIR']
	except KeyError:
		logger.fatal("$TMPDIR variable not set. Aborting...")
		pytok.exit(5, logger)

	initScript = arguments[0]
	phastCommand = arguments[1]
	phastOptions = arguments[2]
	phastHist = arguments[3]
	phastOut = arguments[4]
	(initScriptPath, initScriptName) = initScript.rsplit('/', 1)

	loggerString = "Found the following items in the job argument file:\n"
	loggerString += "Initialization script	  : " + initScript + "\n"
	loggerString += "Phast command            : " + phastCommand + "\n"
	loggerString += "Phast options            : " + phastOptions + "\n"
	loggerString += "Phast histogram path     : " + phastHist + "\n"
	loggerString += "Phast mDST path          : " + phastHist + "\n"
	logger.info(loggerString[:-1])

	loggerString = "Writing input files paths to temporary filelist.txt file:\n"
	with open(workingDirectory + "/filelist.txt", 'w') as fileList:
		for line in arguments[5:]:
			fileList.write(line + "\n")
			loggerString += "Wrote '" + line + "'.\n"
	logger.debug(loggerString[:-1])

	command = ""
	command += "cd " + initScriptPath + "\n"
	command += ". " + initScriptName + "\n"
	command += "cd " + workingDirectory + "\n"
	command += phastCommand + " " + phastOptions
	if phastHist != "":
		(phastHistDir, phastHistName) = phastHist.rsplit('/', 1)
		command += " -h " + phastHistName
	if phastOut != "":
		(phastOutDir, phastOutName) = phastOut.rsplit('/', 1)
		command += " -o " + phastOutName
	command += " -l filelist.txt\n"
	if phastHist != "":
		(phastHistDir, phastHistName) = phastHist.rsplit('/', 1)
		command += "rfcp " + phastHistName + " " + phastHistDir + "\n"
	if phastOut != "":
		(phastOutDir, phastOutName) = phastOut.rsplit('/', 1)
		command += "rfcp " + phastOutName + " " + phastOutDir + "\n"

	pytok.runCommand("Running phast", command, logger)

	command = "du -hs " + workingDirectory
	pytok.runCommand("Getting local disk usage", command, logger)

	endingTime = datetime.datetime.now()

	logger.info("Job ended.")
	logger.info("Job running time: " + str(endingTime - startingTime))

	pytok.exit(0, logger)
