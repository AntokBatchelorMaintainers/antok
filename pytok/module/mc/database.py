
import datetime
import sqlite3

import mcBatchLib


class Database(object):

	_dbCursor = None
	_logger = None

	def __init__(self, filename, logger, createTable = False):
		self._conn = sqlite3.connect(filename)
		self._dbCursor = self._conn.cursor()
		self._logger = logger
		if createTable:
			self._dbCursor.execute('''CREATE TABLE configFiles (
			                            id INTEGER PRIMARY KEY,
			                            scriptConf BLOB,
			                            generatorConf BLOB,
			                            detectorsDat BLOB,
			                            comgeantOpts BLOB,
			                            comgeantGeom22 BLOB,
			                            comgeantGeom23 BLOB,
			                            comgeantGeom24 BLOB,
			                            coralOpts BLOB
			                         )''')
			self._dbCursor.execute('''CREATE TABLE mcJobs (
			                            id INTEGER PRIMARY KEY,
			                            rangeLow INTEGER,
			                            rangeHigh INTEGER,
			                            description TEXT,
			                            creationDate TIMESTAMP,
			                            configId INTEGER,
			                            FOREIGN KEY(configId) REFERENCES configFiles(id)
			                         )''')
			self._dbCursor.execute('''CREATE TABLE filteredFiles (
			                            id INTEGER PRIMARY KEY,
			                            histogramPath TEXT,
			                            outputFilePath TEXT,
			                            creationDate TIMESTAMP,
			                            status INTEGER
			                         )''')
			self._dbCursor.execute('''CREATE TABLE mDSTs (
			                            id INTEGER PRIMARY KEY,
			                            path TEXT,
			                            filteredFileId INTEGER,
			                            mcJobId INTEGER,
			                            FOREIGN KEY(filteredFileId) REFERENCES filteredFiles(id),
			                            FOREIGN KEY(mcJobId) REFERENCES mcJobs(id)
			                         )''')

	def commit(self):
		self._conn.commit()

	def getFreeSlots(self, nJobs = 1):
		ranges = []
		for row in self._dbCursor.execute("SELECT rangeLow, rangeHigh FROM mcJobs"):
			ranges.append([x for x in row])
		self._logger.debug("Found ranges '" + str(ranges) + "'.")
		if len(ranges) == 0:
			return 1
		ranges.sort(key=lambda x: x[0])
		freeRegions = []
		for i in range(len(ranges) - 1):
			freeRegions.append([ranges[i][1], ranges[i+1][0]])
		self._logger.debug("Found free regions '" + str(freeRegions) + "'.")
		for freeRegion in freeRegions:
			if freeRegion[1] - freeRegion[0] - 1 > nJobs:
				self._logger.debug("Lower job id: " + str(freeRegion[0] + 1))
				return freeRegion[0] + 1
		self._logger.debug("Lower job id: " + str(ranges[-1][1] + 1))
		return ranges[-1][1] + 1

	def insertConfigFiles(self, config):
		currentConfigs = mcBatchLib.compressConfigFiles(config)
		nConfigs = int(self._dbCursor.execute('''SELECT COUNT(*) FROM configFiles WHERE
		                                          scriptConf=:scriptConf AND
		                                          detectorsDat=:detectorsDat AND
		                                          generatorConf=:generatorConf AND
		                                          comgeantOpts=:comgeantOpts AND
		                                          comgeantGeom22=:comgeantGeom22 AND
		                                          comgeantGeom23=:comgeantGeom23 AND
		                                          comgeantGeom24=:comgeantGeom24 AND
		                                          coralOpts=:coralOpts''',
		                                      currentConfigs).fetchone()[0])
		if nConfigs == 0:
			self._logger.warningPrompt("New set of config files detected.")
			self._dbCursor.execute('''INSERT INTO configFiles VALUES
			                           (NULL,
			                            :scriptConf,
			                            :generatorConf,
			                            :detectorsDat,
			                            :comgeantOpts,
			                            :comgeantGeom22,
			                            :comgeantGeom23,
			                            :comgeantGeom24,
			                            :coralOpts)''',
			                       currentConfigs)
		else:
			self._logger.info("Config files found in database.")
		return self._dbCursor.execute('''SELECT id FROM configFiles WHERE
		                                  scriptConf=:scriptConf AND
		                                  generatorConf=:generatorConf AND
		                                  detectorsDat=:detectorsDat AND
		                                  comgeantOpts=:comgeantOpts AND
		                                  comgeantGeom22=:comgeantGeom22 AND
		                                  comgeantGeom23=:comgeantGeom23 AND
		                                  comgeantGeom24=:comgeantGeom24 AND
		                                  coralOpts=:coralOpts''',
		                              currentConfigs).fetchone()[0]

	def insertJob(self, firstJobId, lastJobId, description, configId):
		values = (None, firstJobId, lastJobId, description, datetime.datetime.now(), configId)
		self._dbCursor.execute("INSERT INTO mcJobs VALUES (?, ?, ?, ?, ?, ?)", values)

	def getConfigFilesForJob(self, jobID):
		configID = self._dbCursor.execute("SELECT * FROM mcJobs WHERE rangeLow<=(?) AND rangeHigh>=(?)", (jobID, jobID)).fetchone()
		self._logger.debug("Job " + str(jobID) + " is part of set " +
		                   str(configID[0]) + " (range [" + str(configID[1]) +
		                   ", " + str(configID[2]) + "]), created on " + configID[4] + ".")
		if configID is None:
			self._logger.error("JobID '" + str(jobID) + "' not found in database.")
			return {}
		configFiles = self._dbCursor.execute("SELECT * FROM configFiles WHERE id=(?)", (configID[5],)).fetchone()
		if configFiles is None:
			self._logger.error("Entry in configFiles table with id=" + str(configID[5]) + " does not exist for jobID=" + str(jobID))
			return {}
		files = {"scriptConf": configFiles[1],
		         "generatorConf": configFiles[2],
		         "detectorsDat": configFiles[3],
		         "comgeantOpts": configFiles[4],
		         "comgeantGeom22": configFiles[5],
		         "comgeantGeom23": configFiles[6],
		         "comgeantGeom24": configFiles[7],
		         "coralOpts": configFiles[8]
		        }
		return mcBatchLib.decompressConfigFiles(files)

	def getJobInfo(self, jobID):
		data = self._dbCursor.execute("SELECT * FROM mcJobs WHERE rangeLow<=(?) AND rangeHigh>=(?)", (jobID, jobID)).fetchone()
		if data is None:
			self._logger.error("JobID '" + str(jobID) + "' not found in database.")
			return {}
		retval = {"id": data[0],
		          "rangeLow": data[1],
		          "rangeHigh": data[2],
		          "description": data[3],
		          "creationDate": data[4],
		          "configId": data[5],
		         }
		return retval

	def getJobScriptConf(self, jobID):
		configID = self._dbCursor.execute("SELECT configId FROM mcJobs WHERE rangeLow<=(?) AND rangeHigh>=(?)", (jobID, jobID)).fetchone()
		if configID is None:
			self._logger.error("JobID '" + str(jobID) + "' not found in database.")
			return None
		configFile = self._dbCursor.execute("SELECT scriptConf FROM configFiles WHERE id=(?)", (configID[0],)).fetchone()
		configFile = {"scriptConf": configFile[0]}
		return mcBatchLib.decompressConfigFiles(configFile)["scriptConf"]

	def mDSTPathAlreadyRegistered(self, path):
		mDST = self._dbCursor.execute("SELECT * FROM mDSTs WHERE path=(?)", [path]).fetchone()
		return not mDST is None

	def insertmDST(self, path, mDSTSuffix):
		mDSTFileName = path.rsplit('/', 1)[-1]
		mDSTNumber = None
		try:
			mDSTNumber = int(mDSTFileName.rstrip(mDSTSuffix))
		except ValueError:
			self._logger.warning("Could not extract mDST's id. Setting its JobID to \"NULL\"")
			mDSTNumber = None
		jobId = None
		if mDSTNumber is not None:
			jobInfo = self.getJobInfo(mDSTNumber)
			if jobInfo == {}:
				self._logger.warning("Could not get the JobID for mDST " + str(mDSTNumber) + ". Setting it to \"NULL\"")
			else:
				jobId = jobInfo["id"]
		self._dbCursor.execute("INSERT INTO mDSTs VALUES (NULL, (?), NULL, (?))", (path, jobId))

	def getUnsubmittedmDSTs(self):
		retval = []
		for row in self._dbCursor.execute("SELECT id FROM mDSTs WHERE filteredFileId IS NULL"):
			retval.append(row[0])
		return retval

	def getmDSTPathsForJob(self, jobId):
		retval = []
		for row in self._dbCursor.execute("SELECT path FROM mDSTs WHERE filteredFileId=(?)", (jobId,)):
			retval.append(row[0])
		return retval

	def getNewFilteredFile(self, mDSTs, histPath="", outPath="", histName = "", outName = ""):
		now = datetime.datetime.now()
		self._dbCursor.execute("INSERT INTO filteredFiles VALUES (NULL, NULL, NULL, (?), 0)", [now])
		newJobId = int(self._dbCursor.execute("SELECT id FROM filteredFiles WHERE creationDate=(?)", [now]).fetchone()[0])
		update = [(newJobId, int(x)) for x in mDSTs]
		if histPath != "":
			dot = ""
			if histName != "":
				dot = "."
			histogramPath = histPath + "/" + histName + dot + str(newJobId) + ".root"
		if outPath != "":
			dot = ""
			if outName != "":
				dot = "."
			outputFilePath = outPath + "/" + outName + dot + str(newJobId) + ".root"
		self._dbCursor.executemany("UPDATE mDSTs SET filteredFileId=(?) WHERE id=(?)", update)
		self._dbCursor.execute('''UPDATE filteredFiles
		                              SET histogramPath=(?),
		                                  outputFilePath=(?)
		                              WHERE id=(?)''', (histogramPath, outputFilePath, newJobId))
		return (newJobId, histogramPath, outputFilePath)

	def setFilteringJobStatus(self, jobId, status):
		self._dbCursor.execute("UPDATE filteredFiles SET status=(?) WHERE id=(?)", (status, jobId))

	def getFilteringJobsWithStatus(self, status):
		retval = []
		for row in self._dbCursor.execute("SELECT id FROM filteredFiles WHERE status=(?)", (status,)):
			retval.append(row[0])
		return retval

	def getFilteringJobInfo(self, jobId):
		data = self._dbCursor.execute("SELECT * FROM filteredFiles WHERE id=(?)", (jobId,)).fetchone()
		if data is None:
			self._logger.error("FilteringJobID '" + str(jobId) + "' not found in database.")
			return {}
		retval = {"id": data[0],
		          "histogramPath": data[1],
		          "outputFilePath": data[2],
		          "creationDate": data[3],
		          "status": data[4]
		         }
		return retval
