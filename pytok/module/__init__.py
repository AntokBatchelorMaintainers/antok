
import libPytok as cpp
import subprocess as _subprocess
import sys as _sys

from _logger import Logger


def exit(exitCode, logger):
	if exitCode == 0 and not logger.printSummary():
		exitCode = 255
	_sys.exit(exitCode)


def runCommand(name, commandString, logger):

	logger.info(name + " command String: \n" + commandString)
	commandString = "errHandler() { (( errcount++ )); }; trap errHandler ERR\n" + commandString.rstrip('\n') + "\nexit $errcount"
	process = _subprocess.Popen(commandString,
	                            shell=True,
	                            stdout=_subprocess.PIPE,
	                            stderr=_subprocess.STDOUT,
	                            executable="/bin/bash")
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
		_sys.stdout.write(question + prompt)
		choice = raw_input().lower()
		if default is not None and choice == '':
			return retvals[default]
		elif choice.lower() in retvals:
			return retvals[choice.lower()]
		else:
			print("Please respond with 'yes' or 'no' (or 'y' or 'n').")
