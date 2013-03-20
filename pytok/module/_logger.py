
import datetime
import sys

import pytok

class Logger(object):

	dateFormat = "%Y-%m-%d %H:%M:%S"
	_level = 0
	_printingBookkeeping = []

	def __init__(self, level):
		_level = level.lower()
		self._printingBookkeeping = [0, 0, 0, 0, 0]
		if _level == "debug":
			self._level = 0
		elif _level == "info":
			self._level = 1
		elif _level == "warning":
			self._level = 2
		elif _level == "error":
			self._level = 3
		elif _level == "fatal":
			self._level = 4
		else:
			self._level = 0
			self.error("Logging level '" + str(level) + "' not supported.\nLogging switched to level 'Debug'.")

	def _stringHandler(self, msg, level):
		msgs = msg.split('\n')
		string = "[" + str(datetime.datetime.now().strftime(self.dateFormat)) + "] " + level
		string = '%-30s %s' % (string, (': ' + msgs[0]))
		for line in msgs[1:]:
			string += '%-33s %s' % ('\n', line)
		string = string.rstrip('\n').rstrip(' ').rstrip('\n')
		return string

	def printSummary(self):
		string = ""
		string += "Debug messages   : " + str(self._printingBookkeeping[0])
		if self._level > 0:
			string += " (not printed)"
		string += "\n"
		string += "Info messages    : " + str(self._printingBookkeeping[1])
		if self._level > 1:
			string += " (not printed)"
		string += "\n"
		string += "Warning messages : " + str(self._printingBookkeeping[2])
		if self._level > 2:
			string += " (not printed)"
		string += "\n"
		string += "Error messages   : " + str(self._printingBookkeeping[3])
		if self._level > 3:
			string += " (not printed)"
		string += "\n"
		string += "Fatal messages   : " + str(self._printingBookkeeping[4])
		if self._level > 4:
			string += " (not printed)"
		string += "\n"
		self.info(string)
		nBads = (self._printingBookkeeping[2] +
		         self._printingBookkeeping[3] +
		         self._printingBookkeeping[4])
		return (nBads == 0)

	def debug(self, msg):
		self._printingBookkeeping[0] += 1
		if self._level < 1:
			print(self._stringHandler(msg, 'DEBUG'))

	def info(self, msg):
		self._printingBookkeeping[1] += 1
		if self._level < 2:
			print(self._stringHandler(msg, 'INFO'))

	def warning(self, msg):
		self._printingBookkeeping[2] += 1
		if self._level < 3:
			print(self._stringHandler(msg, 'WARNING'))

	def error(self, msg):
		self._printingBookkeeping[3] += 1
		if self._level < 4:
			print(self._stringHandler(msg, 'ERROR'))

	def fatal(self, msg):
		self._printingBookkeeping[4] += 1
		if self._level < 5:
			print(self._stringHandler(msg, 'FATAL'))

	def warningPrompt(self, msg):
		self.warning(msg)
		self._printingBookkeeping[2] -= 1
		if not pytok.promptYesOrNo("Continue?"):
			self.fatal("Aborting because of user decision...")
			self.printSummary()
			sys.exit(100)
