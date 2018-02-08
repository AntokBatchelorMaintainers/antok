#!/usr/bin/env python
# coding: utf-8
'''
Created on Sat 09 Jan 2016 12:05:35 AM CET


@author: Stefan Wallner
'''

program_description = '''
    execute the tree reader
'''

# std includes
import sys
from optparse import OptionParser
import os
import subprocess as sp
import tempfile
import atexit
import signal
import time

# add path to our own pytok
sys.path.insert(0,os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))),"pyLib"))
import pytok.yamltok as yaml


def main():
	optparser = OptionParser( usage="Usage:run <args> [<options>]", description = program_description );
	optparser.add_option('-c', '--configfile', dest='configfile', action='store', type='str', default="", help="Config file for antok.")
	optparser.add_option('-o', '--outfile', dest='outfile', action='store', type='str', default="antok.root",
	                     help="Merge files to the given output file. If given, all sub-files will be stored in a folder named according to the output file.")

	( options, args ) = optparser.parse_args();

	if not options.configfile:
		print "No configfile given"
		print optparser.usage
		exit( 100 );

	if not os.path.isfile( options.configfile ):
		print "Config file '{0}' not found!".format(options.configfile)
		print optparser.usage
		exit( 100 );

	if options.outfile and os.path.isfile( options.outfile ):
		print "Output file '{0}' exists found!".format(options.outfile)
		print optparser.usage
		exit( 100 );

	options.outfile = os.path.realpath(options.outfile)

	# make my own copy of the config file and use it to be protected from changes during the execution
	# load and store it to handle !include statements in the config file
	local_configfile = tempfile.mktemp(suffix=".yaml", prefix="antok_run")
	try:
		with open(options.configfile, "r") as fin:
			config = yaml.load(fin, yaml.LoaderTok)
			with open(local_configfile, "w") as fout:
				# merge lists
				yaml.merge_lists(config)
				yaml.dump(config, fout)
	except Exception as e:
		print e
		print "Cannot read config"
		sys.exit(1)

	options.configfile = local_configfile
	atexit.register(lambda: os.remove(options.configfile))

	print "==========================================================================="
	print "========================  START  =========================================="
	print "==========================================================================="

	bin = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'treereader')
	infiles = "' '".join(args)
	proc = sp.Popen("{bin} '{infiles}' '{outfile}' '{config}'".format(bin=bin,
	                                                                  infiles = infiles,
	                                                                  outfile = options.outfile,
	                                                                  config = options.configfile),
	               shell=True)


	def signalHandlerKill(sig, frame):
		proc.kill()
		sys.exit(1)
	def signalHandlerWait(sig, frame):
		signal.signal(signal.SIGINT, signalHandlerKill)
		proc.communicate() # wait to finish
	signal.signal(signal.SIGINT, signalHandlerWait)

	proc.communicate()

	if proc.poll() > 0:
		print "==========================================================================="
		print "treereader stopped with exit status: ", proc.poll()
		sys.exit(proc.poll())

	print "==========================================================================="
	print "========================  FINISHED  ======================================="
	print "==========================================================================="

	sys.exit(0)

if __name__ == '__main__':
	main()




