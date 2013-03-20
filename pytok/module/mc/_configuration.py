
import os.path


class ConfigurationException(Exception):

	text = ""

	def __init__(self, section, option):
		self.text = "[" + str(section) + "] " + str(option)

	def __str__(self):
		return self.text


class Configuration(object):

	_configFile = None

	_logger = None

	comgeant = ""
	comgeantBuildDir = ""
	coral = ""
	initScript = ""
	detectorsDat = ""
	coralDico = ""
	eventsDir = ""
	outputDir = ""
	oneEvent = ""
	coralOptions = ""
	comgeantOptions = ""
	comgeantGeom22 = ""
	comgeantGeom23 = ""
	comgeantGeom24 = ""
	nEvents = -1
	mDSTFilenameSuffix = ""
	_taskID = -1
	randomSeed = -1
	wDir = ""
	compassFiles = ""
	phastJobScript = ""
	phastCommand = ""
	phastOptions = ""
	phastHistDir = ""
	phastOutDir = ""
	phastHistPrefix = ""
	phastOutPrefix = ""
	phastSubmissionCommand = ""
	phastCheckForRunningJobs = ""
	phastmDSTsPerJob = -1
	phastSubmitterIntervall = -1
	generateInputYourself = False
	generatorCommand = ""
	generatorOptions = ""
	generatorOutputFileSwitch = ""
	generatorEventNumberSwitch = ""
	generatorSeedSwitch = ""
	_generatorFilesToWatchString = ""
	generatorFilesToWatch = []
	submissionCommand = ""
	submissionJobNameSwitch = ""
	submissionArraySwitch = ""
	submissionOutputDirSwitch = ""
	submissionJobScript = ""
	bookkeepingDir = ""
	dbFile = ""

	def __init__(self, logger):
		self._logger = logger

	def _getTaskID(self):
		return self._taskID

	def _setTaskID(self, taskID):
		self._taskID = taskID
		self.randomSeed = 1000000 - self._taskID

	taskID = property(_getTaskID, _setTaskID)

	def _checkSection(self, section):
		if not self._configFile.has_section(section):
			self._logger.error("Could not find '" + section + "' section in config file.")
			return False
		return True

	def _getOption(self, section, option, optional = False):
		if not self._configFile.has_option(section, option):
			if not optional:
				self._logger.error("Could not find '" + option + "' option in '" + section + "' section.")
				raise(ConfigurationException(section, option))
			else:
				self._logger.warning("Could not find optional '" + option + "' option in '" + section +
				               "' section.\nSet it to an empty string to avoid this warning.")
				return ""
		return self._configFile.get(section, option)

	def read(self, configFile):

		self._logger.debug("Reading config file.")

		self._configFile = configFile
		if not self._checkSection("Submission") and \
		       self._checkSection("Environment") and \
		       self._checkSection("OptionFiles") and \
		       self._checkSection("Parameters") and \
		       self._checkSection("Generator"):
			return False;
		try:
			self.comgeant = self._getOption("Environment", "comgeant")
			self.comgeantBuildDir = self._getOption("Environment", "comgeant_build_directory")
			self.coral = self._getOption("Environment", "coral")
			self.initScript = self._getOption("Environment", "init_script")
			self.detectorsDat = self._getOption("Option Files", "detectors_dat", True)
			self.coralDico = self._getOption("Option Files", "coral_dico")
			self.eventsDir = self._getOption("Environment", "events_directory")
			self.outputDir = self._getOption("Environment", "output_directory")
			self.oneEvent = self._getOption("Environment", "1_event_file", True)
			self.coralOptions = self._getOption("Option Files", "coral_options")
			self.comgeantOptions = self._getOption("Option Files", "comgeant_options")
			self.comgeantGeom22 = self._getOption("Option Files", "comgeant_specific_geom_file")
			self.comgeantGeom23 = self._getOption("Option Files", "comgeant_specific+_geom_file")
			self.comgeantGeom24 = self._getOption("Option Files", "comgeant_specific++_geom_file", True)
			self.nEvents = int(self._getOption("Parameters", "number_of_events"))
			self.mDSTFilenameSuffix = self._getOption("Parameters", "mdst_filename_suffix")
			self.phastJobScript = self._getOption("Phast", "phast_job_script", True)
			self.phastCommand = self._getOption("Phast", "phast_executable", True)
			self.phastOptions = self._getOption("Phast", "phast_options", True)
			self.phastHistDir = self._getOption("Phast", "hist_dir", True)
			self.phastOutDir = self._getOption("Phast", "output_dir", True)
			self.phastHistPrefix = self._getOption("Phast", "hist_file_prefix", True)
			self.phastOutPrefix = self._getOption("Phast", "out_file_prefix", True)
			self.phastSubmissionCommand = self._getOption("Phast", "submission_command", True)
			self.phastCheckForRunningJobs = self._getOption("Phast", "check_for_running_jobs", True)
			self.phastmDSTsPerJob = self._getOption("Phast", "mDSTs_per_job", True)
			if self.phastmDSTsPerJob == "":
				self.phastmDSTsPerJob = -1
			else:
				self.phastmDSTsPerJob = int(self.phastmDSTsPerJob)

			self.phastSubmitterIntervall = self._getOption("Phast", "loop_interval", True)
			if self.phastSubmitterIntervall == "":
				self.phastSubmitterIntervall = -1
			else:
				self.phastSubmitterIntervall = int(self.phastSubmitterIntervall)

			self.generatorCommand = self._getOption("Generator", "generator_command", True)
			self.generatorOptions = self._getOption("Generator", "generator_options", True).replace('\n', ' ')
			if self.generatorCommand != "":
				self.generatorOutputFileSwitch = self._getOption("Generator", "output_file_switch")
				self.generatorEventNumberSwitch = self._getOption("Generator", "event_number_switch")
				self.generatorSeedSwitch = self._getOption("Generator", "seed_switch")
				self._generatorFilesToWatchString = self._getOption("Generator", "files_to_watch", True)
			self.submissionCommand = self._getOption("Submission", "command", True)
			self.submissionJobNameSwitch = self._getOption("Submission", "job_name_switch", True)
			self.submissionArraySwitch = self._getOption("Submission", "array_switch", True)
			self.submissionOutputDirSwitch = self._getOption("Submission", "output_dir_switch", True)
			self.submissionJobScript = self._getOption("Submission", "job_script", True)
			self.bookkeepingDir = self._getOption("Submission", "bookkeeping_directory", True)
			self.dbFile = self._getOption("Submission", "database_file", True)
		except ConfigurationException:
			return False
		return True

	def check(self, noJob = False):
		self._logger.debug("Checking config file.")
		error = False
		if not os.path.isdir(self.comgeant):
			self._logger.error("Comgeant directory '" + str(self.comgeant) + "' does not exist or is no directory.")
			error = True
		if not os.path.isdir(self.comgeantBuildDir):
			self._logger.error("Comgeant build directory '" + str(self.comgeantBuildDir) + "' does not exist or is no directory.")
			error = True
		if not os.path.isfile(self.comgeantBuildDir + "/omgbatch"):
			self._logger.error("'omgbatch' not found in '" + str(self.comgeantBuildDir) + "'.")
			error = True
		if not os.path.isfile(self.comgeantBuildDir + "/omginter"):
			self._logger.error("'omginter' not found in '" + str(self.comgeantBuildDir) + "'.")
			error = True
		if not os.path.isfile(self.coral + "/coral.exe"):
			self._logger.error("'coral.exe' not found in '" + str(self.coral) + "'.")
			error = True
		if not os.path.isfile(self.initScript):
			self._logger.error("The init_script '" + str(self.initScript) + "' does not appear to exist.")
			error = True
		if not os.path.isdir(self.eventsDir):
			self._logger.error("The directory with the generated events '" + str(self.eventsDir) + "' does not appear to exist.")
			error = True
		if self.oneEvent == "" and self.detectorsDat == "":
			self._logger.error("Either 'detectors_dat' or '1_event_file' have to be set.")
			error = True
		else:
			if self.oneEvent != "" and not os.path.isfile(self.oneEvent):
				self._logger.error("Could not get the file with one generated event '" + str(self.oneEvent) + "'.")
				error = True
			if self.detectorsDat != "" and not os.path.isfile(self.detectorsDat):
				self._logger.error("Could not find detectors.dat '" + self.detectorsDat + "'.")
				error = True
		if not os.path.isfile(self.coralOptions):
			self._logger.error("Could not find coral options file '" + str(self.coralOptions) + "'.")
			error = True
		if not os.path.isfile(self.comgeantGeom22):
			self._logger.error("Could not find comgeant geometry file (22) '" + str(self.comgeantGeom22) + "'.")
			error = True
		if not os.path.isfile(self.comgeantGeom23):
			self._logger.error("Could not find comgeant geometry file (23) '" + str(self.comgeantGeom23) + "'.")
			error = True
		if self.comgeantGeom24 != "" and not os.path.isfile(self.comgeantGeom24):
			self._logger.error("Could not find comgeant geometry file (24) '" + str(self.comgeantGeom24) + "'.")
			error = True
		if not os.path.isfile(self.coralDico):
			self._logger.error("Could not find coral dico file '" + str(self.coralDico) + "'.")
			error = True
		else:
			if self.detectorsDat != "":
				if os.path.getmtime(self.coralDico) < os.path.getmtime(self.detectorsDat):
					self._logger.error("Coral dico file is older than detectors.dat. Regenerate dico file.")
					error = True
			else:
				if os.path.getmtime(self.coralDico) < os.path.getmtime(self.comgeantGeom23):
					self._logger.error("Coral dico file is older than comgeant geometry file. Regenerate dico file.")
					error = True
		if not os.path.isdir(self.outputDir):
			if os.system("rfdir " + self.outputDir):
				self._logger.error("Output directory '" + self.outputDir + "' not found.")
				error = True
		if self._generatorFilesToWatchString != "":
			self.generatorFilesToWatch = self._generatorFilesToWatchString.replace('\n', ' ').replace(',', ' ').split()
			for watchFile in self.generatorFilesToWatch:
				if not os.path.isfile(watchFile):
					self._logger.error("Could not find file '" + str(watchFile) + "', which was supposed to exist for the generator.")
					error = True
		if not noJob:
			if self._taskID < 0:
				self._logger.error("$SGE_TASK_ID '" + self._taskID + "' is invalid.")
				error = True
			if not os.path.isfile(self.eventsDir + "/" + str(self._taskID) + ".fort.26"):
				if self.generatorCommand == "":
					self._logger.error("Input file '" + self.eventsDir + "/" + str(self._taskID) + ".fort.26' not found and no generator specified.")
					error = True
				else:
					if not os.path.isfile(self.generatorCommand):
						self._logger.error("Generator executable '" + self.generatorCommand + "' not found.")
						error = True
					else:
						self.generateInputYourself = True
			if os.path.isfile(self.outputDir + "/" + str(self._taskID) + ".mDST"):
				self._logger.error("Output file '" + self.outputDir + "/" + str(self._taskID) + ".mDST' already exists.")
				error = True
			if not os.path.isdir(self.wDir):
				self._logger.error("Working directory '" + str(self.wDir) + "' is not a directory.")
				error = True
			if self.randomSeed != 1000000 - self._taskID:
				self._logger.error("Internal inconsistency with the random seed and the $SGE_TASK_ID.")
				error = True
		self._logger.debug("Config file checked.")
		return not error

	def checkSubmitterOptions(self):
		try:
			self.submissionCommand = self._getOption("Submission", "command")
			self.submissionJobNameSwitch = self._getOption("Submission", "job_name_switch")
			self.submissionArraySwitch = self._getOption("Submission", "array_switch")
			self.submissionOutputDirSwitch = self._getOption("Submission", "output_dir_switch")
			self.submissionJobScript = self._getOption("Submission", "job_script")
			self.bookkeepingDir = self._getOption("Submission", "bookkeeping_directory")
			self.dbFile = self._getOption("Submission", "database_file")
		except ConfigurationException:
			return False
		error = False
		if not os.path.isdir(self.bookkeepingDir):
			self._logger.error("Could not find bookkeeping directory '" + self.bookkeepingDir + "'.")
			error = True
		return not error

	def checkPhastSubmitterOptions(self):
		try:
			self.phastJobScript = self._getOption("Phast", "phast_job_script")
			self.phastCommand = self._getOption("Phast", "phast_executable")
			self.phastOptions = self._getOption("Phast", "phast_options")
			self.phastCheckForRunningJobs = self._getOption("Phast", "check_for_running_jobs")
			self.phastSubmissionCommand = self._getOption("Phast", "submission_command")
			self.phastmDSTsPerJob = int(self._getOption("Phast", "mDSTs_per_job"))
			self.phastSubmitterIntervall = int(self._getOption("Phast", "loop_interval"))
		except ConfigurationException:
			return False
		error = False
		if not os.path.isfile(self.phastJobScript):
			self._logger.error("Phast job script '" + self.phastJobScript + "' not found.")
			error = True
		if not os.path.isfile(self.phastCommand):
			self._logger.error("Phast executable '" + self.phastCommand + "' not found.")
			error = True
		if self.phastHistDir == "" and self.phastOutDir == "":
			self._logger.error("Phast histogram directory and output directory cannot both be empty.")
			error = True
		else:
			if self.phastHistDir != "":
				if not os.path.isdir(self.phastHistDir):
					if os.system("rfdir " + self.phastHistDir + " &>/dev/null"):
						self._logger.error("Could not find phast histogram directory '" + self.phastHistDir + "'.")
						error = True
			if self.phastOutDir != "":
				if not os.path.isdir(self.phastOutDir):
					if os.system("rfdir " + self.phastOutDir + " &>/dev/null"):
						self._logger.error("Could not find phast output directory '" + self.phastOutDir + "'.")
						error = True
		return not error

	def __str__(self):
		retval  = ""
		retval += "Comgeant                           : " + self.comgeant + "\n"
		retval += "Comgeant executables               : " + self.comgeantBuildDir + "\n"
		retval += "Coral                              : " + self.coral + "\n"
		retval += "Init script                        : " + self.initScript + "\n"
		retval += "detectors.dat                      : " + self.detectorsDat + "\n"
		retval += "Coral dico file                    : " + self.coralDico + "\n"
		retval += "Events directory                   : " + self.eventsDir + "\n"
		retval += "Output directory                   : " + self.outputDir + "\n"
		retval += "1 Event File                       : " + self.oneEvent + "\n"
		retval += "Comgeant options                   : " + self.comgeantOptions + "\n"
		retval += "Comgeant geometry (22)             : " + self.comgeantGeom22 + "\n"
		retval += "Comgeant geometry (23)             : " + self.comgeantGeom23 + "\n"
		retval += "Comgeant geometry (24)             : " + self.comgeantGeom24 + "\n"
		retval += "Coral option file                  : " + self.coralOptions + "\n"
		retval += "Phast job script                   : " + self.phastJobScript + "\n"
		retval += "Phast executable                   : " + self.phastCommand + "\n"
		retval += "Phast options                      : " + self.phastOptions + "\n"
		retval += "Phast histogram directory          : " + self.phastHistDir + "\n"
		retval += "Phast output directory             : " + self.phastOutDir + "\n"
		retval += "Phast histogram file prefix        : " + self.phastHistPrefix + "\n"
		retval += "Phast output file prefix           : " + self.phastOutPrefix + "\n"
		retval += "Phast submission command           : " + self.phastSubmissionCommand + "\n"
		retval += "Phast mDSTs per job                : " + str(self.phastmDSTsPerJob) + "\n"
		retval += "Phast submitter loop interval      : " + str(self.phastSubmitterIntervall) + "\n"
		retval += "Generator command                  : " + self.generatorCommand + "\n"
		retval += "Generator options                  : " + repr(self.generatorOptions) + "\n"
		retval += "Generator event number switch      : " + self.generatorEventNumberSwitch + "\n"
		retval += "Generator seed switch              : " + self.generatorSeedSwitch + "\n"
		retval += "Generator files to watch           : "
		if len(self.generatorFilesToWatch) == 0:
			retval += ("\n")
		else:
			i = 0
			for watchFile in self.generatorFilesToWatch:
				if i != 0:
					retval += "                                     "
				retval += watchFile + "\n"
				i += 1
		retval += "Number of events                   : " + str(self.nEvents) + "\n"
		retval += "mDST filename suffix               : " + self.mDSTFilenameSuffix + "\n"
		retval += "Random seed                        : " + str(self.randomSeed) + "\n"
		retval += "$SGE_TASK_ID                       : " + str(self._taskID) + "\n"
		retval += "$TMPDIR                            : " + self.wDir + "\n"
		retval += "$COMPASS_FILES                     : " + self.compassFiles + "\n"
		retval += "Submission command                 : " + self.submissionCommand + "\n"
		retval += "Submission job name switch         : " + self.submissionJobNameSwitch + "\n"
		retval += "Submission array switch            : " + self.submissionArraySwitch + "\n"
		retval += "Submission output directory switch : " + self.submissionOutputDirSwitch + "\n"
		retval += "Submission job script              : " + self.submissionJobScript + "\n"
		retval += "Bookkeeping directory              : " + self.bookkeepingDir + "\n"
		retval += "Database file                      : " + self.dbFile + "\n"
		return retval[:-1]
