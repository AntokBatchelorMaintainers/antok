
import bz2
import StringIO
import subprocess
import sys


def exit(exitCode, logger):
	if exitCode == 0 and not logger.printSummary():
		exitCode = 255
	sys.exit(exitCode)


def runCommand(name, commandString, logger):

	logger.info(name + " command String: \n" + commandString)
	commandString = "errHandler() { (( errcount++ )); }; trap errHandler ERR\n" + commandString.rstrip('\n') + "\nexit $errcount"
	process = subprocess.Popen(commandString, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
	(output, stderr) = process.communicate()
	errorCode = process.returncode
	if errorCode != 0:
		logger.error(name + " terminated with " + str(errorCode) + " error(s). Output:\n" + output)
		return False
	else:
		logger.info(name + " executed. Output:\n" + output[:-1])
		return True


def promptYesOrNo(question, default = None):
	if default == None:
		prompt = " [y/n] "
	elif default == "yes":
		prompt = " [Y/n] "
	elif default == "no":
		prompt = " [y/N] "
	else:
		raise Exception("Invalid default value '" + default + "'.")
	retvals = {"yes": True, "y": True, "ye": True, "no": False, "n": False}
	while True:
		sys.stdout.write(question + prompt)
		choice = raw_input().lower()
		if default is not None and choice == '':
			return retvals[default]
		elif choice.lower() in retvals:
			return retvals[choice.lower()]
		else:
			print("Please respond with 'yes' or 'no' (or 'y' or 'n').")


def compressConfigFiles(config):
	files = {}
	files["detectorsDat"] = config.detectorsDat
	files["comgeantOpts"] = config.comgeantOptions
	files["comgeantGeom22"] = config.comgeantGeom22
	files["comgeantGeom23"] = config.comgeantGeom23
	files["comgeantGeom24"] = config.comgeantGeom24
	files["coralOpts"] = config.coralOptions
	for inFile in files.keys():
		if files[inFile] == '':
			files[inFile] = buffer('')
		else:
			with open(files[inFile], 'r') as inputFile:
				files[inFile] = buffer(bz2.compress(inputFile.read(), 9))
	memoryFile = StringIO.StringIO()
	config._configFile.write(memoryFile)
	files["scriptConf"] = buffer(bz2.compress(memoryFile.getvalue(), 9))
	generatorConf = ""
	for inputFileName in config.generatorFilesToWatch:
		header = "#          START FILE '" + inputFileName + "'          #"
		separator = len(header) * '#'
		generatorConf += separator + "\n" + header + "\n" + separator + "\n\n\n"
		with open(inputFileName, 'r') as inputFile:
			generatorConf += inputFile.read()
		generatorConf += "\n\n" + separator + "\n"
		generatorConf += "#          END FILE '" + inputFileName + "'            #\n"
		generatorConf += separator + "\n"
	files["generatorConf"] = buffer(bz2.compress(generatorConf, 9))
	return files

def decompressConfigFiles(files):

	for key in files.keys():
		try:
			files[key] = bz2.decompress(str(files[key]))
		except IOError:
			files[key] = ''
	return files
