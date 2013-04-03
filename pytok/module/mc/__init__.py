
import bz2 as _bz2
import StringIO as _StringIO

from _configuration import Configuration
from _database import Database

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
				files[inFile] = buffer(_bz2.compress(inputFile.read(), 9))
	memoryFile = _StringIO.StringIO()
	config._configFile.write(memoryFile)
	files["scriptConf"] = buffer(_bz2.compress(memoryFile.getvalue(), 9))
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
	files["generatorConf"] = buffer(_bz2.compress(generatorConf, 9))
	return files

def decompressConfigFiles(files):

	for key in files.keys():
		try:
			files[key] = _bz2.decompress(str(files[key]))
		except IOError:
			files[key] = ''
	return files
