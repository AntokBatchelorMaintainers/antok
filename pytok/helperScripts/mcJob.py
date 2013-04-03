#!/usr/bin/env python

import ConfigParser
import datetime
import os
import sys

import pytok
import pytok.mc


if __name__ == "__main__":

	if len(sys.argv) != 2:
		print("usage: ./mcJob.py <config_file>")
		sys.exit(1)

	configFile = ConfigParser.ConfigParser()
	configFile.read(sys.argv[1])

	level = "Debug"
	if configFile.has_section("Logging"):
		if configFile.has_option("Logging", "Level"):
			level = configFile.get("Logging", "Level")

	logger = pytok.Logger(level)

	config = pytok.mc.Configuration(logger)
	if not config.read(configFile):
		logger.fatal("Could not read config file. Aborting...")

	try:
		config.wDir = os.environ['TMPDIR']
	except KeyError:
		logger.fatal("$TMPDIR variable not set. Aborting...")
		pytok.exit(5, logger)

	try:
		config.taskID = int(os.environ['SGE_TASK_ID'])
	except KeyError:
		logger.warning("Could not find $SGE_TASK_ID. Using '1'.")
		config.taskID = 1

	try:
		config.compassFiles = os.environ['COMPASS_FILES']
	except KeyError:
		logger.fatal("Could not find $COMPASS_FILES.")
		pytok.exit(10, logger)

	if not config.check():
		logger.fatal("Config file contains errors.")
		pytok.exit(10, logger)

	startingTime = datetime.datetime.now()

	logger.info("Found the following parameters:\n" + str(config))

	environmentString = ""
	for key in os.environ:
		environmentString += '%-20s : %s\n' % (key, os.environ[key])
	logger.debug("Environment as seen by python:\n" + environmentString)

	logger.info("Job starting...")

	# Generate events if necessary
	if config.generateInputYourself:
		logger.info("Generating event file...")
		commandString = config.generatorCommand + " "
		commandString += config.generatorOptions + " "
		commandString += config.generatorOutputFileSwitch + " " + config.eventsDir + "/" + str(config.taskID) + ".gen "
		commandString += config.generatorEventNumberSwitch + " " + str(config.nEvents) + " "
		commandString += config.generatorSeedSwitch + " " + str(config.randomSeed)
		pytok.runCommand("Generating events", commandString, logger)

	# Prepare comgeant directory
	logger.info("Building comgeant preparation command string...")
	commandString = ""
	commandString += "mkdir -v " + config.wDir + "/comgeant\n"
	commandString += "mkdir -v " + config.wDir + "/coral\n"
	commandString += "cd " + config.wDir + "/comgeant\n"
	commandString += "ln -vs " + config.comgeantBuildDir + "/omgbatch\n"
	commandString += "ln -vs " + config.comgeantBuildDir + "/omginter\n"
	commandString += "ln -vs " + config.comgeant + "/data/kinematics/particles/particle_table.ffr fort.16\n"
	commandString += "ln -vs " + config.compassFiles + "/maps/mag_fields/SM1m/SM1M.map.172.data fort.18\n"
	commandString += "ln -vs " + config.compassFiles + "/maps/mag_fields/SM2/FSM.map.5000.data fort.19\n"
	commandString += "ln -vs " + config.comgeant + "/data/geom/geom_general.ffr fort.21\n"
	commandString += "ln -vs " + config.comgeantGeom22 + " fort.22\n"
	commandString += "ln -vs " + config.comgeantGeom23 + " fort.23\n"
	if config.comgeantGeom24 != "":
		commandString += "ln -vs " + config.comgeantGeom24 + " fort.24\n"
	else:
		logger.info("Option 'comgeant_specific++_geom_file' not present or empty, not linking 'fort.24'.")
	commandString += "ln -vs " + config.eventsDir + "/" + str(config.taskID) + ".fort.26 fort.26\n"
	commandString += "cp -v " + config.comgeantOptions + " fort.15\n"
	pytok.runCommand("Comgeant preparation", commandString, logger)

	if config.detectorsDat == "":
		# Produce detectors.dat
		with open(config.wDir + "/comgeant/fort.15", 'a') as fort15:
			fort15.write("TRIG 1")
		commandString = "cd " + config.wDir + "/comgeant\n"
		commandString += "ln -vs " + config.oneEvent + " fort.26\n"
		commandString += "touch fort.32\n"
		commandString += "./omgbatch\n"
		commandString += "mv -v detectors.dat ../coral\n"
		commandString += "cp " + config.comgeantOptions + " fort.15\n"
		commandString += "rm -v fort.32\n"
		commandString += "rm -v fort.26\n"
		commandString += "ln -vs " + config.eventsDir + "/" + str(config.taskID) + ".fort.26 fort.26\n"
		pytok.runCommand("Producing detectors.dat", commandString, logger)
	else:
		# Link detectors.dat
		commandString = "ln -vs " + config.detectorsDat + " " + config.wDir + "/coral/detectors.dat\n"
		pytok.runCommand("Linking detectors.dat '" + config.detectorsDat + "'", commandString, logger)

	# Prepare coral
	with open(config.wDir + "/comgeant/fort.15", 'a') as fort15:
			fort15.write("TRIG " + str(config.nEvents) + "\n")
			fort15.write("RNDMSEQ " + str(config.randomSeed) + "\n")
	commandString = "cd " + config.wDir + "/coral\n"
	commandString += "cp -v " + config.coralOptions + " coralOptions.opt\n"
	commandString += "ln -vs " + config.coralDico + " dicofit.out\n"
	commandString += "ln -vs " + config.coral + "/coral.exe coral.exe\n"
	commandString += "mkfifo eventPipe.fz\n"
	commandString += "ln -vs " + config.wDir + "/coral/eventPipe.fz " + config.wDir + "/comgeant/fort.32\n"
	pytok.runCommand("Preparing coral", commandString, logger)
	with open(config.wDir + "/coral/coralOptions.opt", 'a') as coralOptions:
		coralOptions.write("Monte Carlo file           " + config.wDir + "/coral/eventPipe.fz\n")
		coralOptions.write("detector table             " + config.wDir + "/coral/detectors.dat\n")
		coralOptions.write("TraF Dicofit               " + config.wDir + "/coral/dicofit.out\n")
		coralOptions.write("TraF Graph         [0]      0\n")
		coralOptions.write("mDST file                  " + config.wDir + "/" + str(config.taskID) + config.mDSTFilenameSuffix + "\n")
		coralOptions.write("random number seed         " + str(config.randomSeed))

	commandString = ""
	(initPath, initScriptFile) = config.initScript.rsplit('/', 1)
	commandString += "cd " + initPath + "\n"
	commandString += ". " + initScriptFile + "\n"
	commandString += "cd " + config.wDir + "/comgeant\n"
	commandString += "./omgbatch &>comgeant.out &\n"
	commandString += "cd ../coral\n"
	commandString += "./coral.exe coralOptions.opt &>coral.out\n"
	commandString += "cat ../comgeant/comgeant.out\n"
	commandString += "cat coral.out\n"
	commandString += "cd " + config.wDir + "\n"
	commandString += "rfcp " + str(config.taskID) + config.mDSTFilenameSuffix + " " + config.outputDir + "\n"
	pytok.runCommand("Run Monte Carlo", commandString, logger)

	commandString = "du -hs " + config.wDir
	pytok.runCommand("Getting local disk usage", commandString, logger)

	endingTime = datetime.datetime.now()

	logger.info("Job ended.")
	logger.info("Job running time: " + str(endingTime - startingTime))

	pytok.exit(0, logger)