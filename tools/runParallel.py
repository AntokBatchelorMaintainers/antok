#!/usr/bin/env python
# coding: utf-8
'''
Created on Sat 09 Jan 2016 12:05:35 AM CET


@author: Stefan Wallner
'''

program_description = '''
    Execute the tree reader in parallel on the compute cluster.
'''

# std includes
import sys
from optparse import OptionParser
from collections import defaultdict
import os
import subprocess as sp
import shutil
import re
import multiprocessing
import tempfile
import math
import copy
import time
import batchelor

# add path to our own pytok
sys.path.insert(0,os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))),"pyLib"))
import pytok.yamltok as yaml


def getRunNumber(filename):
	'''
	@return: run number, extracted from mDST filename
	'''
	runnbr = None;
	parsed = re.findall("-([0-9]+)-[0-9]-[0-9]\.root", filename)
	if(parsed):
		runnbr = int( parsed[0])
	else:
		msg = "Can not get run number from file name '{0}'. Using run number 0.".format( filename )
		print msg
		runnbr = 0
	return runnbr

def getSlot(filename):
	'''
	@return: run number, extracted from mDST filename
	'''
	slot = None;
	parsed = re.findall("-[0-9]+-([0-9]-[0-9])\.root", filename)
	if(parsed):
		slot = parsed[0]
	else:
		msg = "Can not get slot from file name '{0}'. Using slot '0-0'".format( filename )
		print msg
		slot='0-0'
	return slot



def mergeRootFiles(output_file, input_files, merge_trees, parallel = True, n_jobs_max = None):

	n_files = len( input_files )

	n_jobs = max( 1, min( n_files/2, int(math.sqrt(n_files)), multiprocessing.cpu_count() ) ) if parallel else 1
	if n_jobs_max:
		n_jobs = min(n_jobs, n_jobs_max)

	q = n_files / n_jobs;
	r = n_files % n_jobs;

	output_files = []
	input_subfiles_list = []
	processes = []
	for i in xrange(n_jobs):
		input_subfiles = input_files[i*q + min(r,i): (i+1)*q + min(r,i+1) ]
		output_subfile = tempfile.mktemp() if n_jobs > 1 else output_file;
		output_files.append( output_subfile );
		input_subfiles_list.append( input_subfiles );
		processes.append( sp.Popen("hadd {opts} {outfile} {infiles} ".format( opts = "" if merge_trees else "-T",
		                                                                      outfile = output_subfile,
		                                                                      infiles = "'" + "' '".join(input_subfiles) + "'" )
		                           ,shell=True, stdout=sp.PIPE, stderr = sp.STDOUT  )
		                 )

	ok = True
	for i,p in enumerate(processes):
		pout, _ = p.communicate()
		status = p.returncode
		if status != 0:
			print "Process exited with exit status", status, "when merging files"
			for f in input_subfiles_list[i]:
				print '\t', f
			for l in pout.split('\n'):
				print '\t', l
			ok = False;

	if ok and n_jobs > 1:
		# merge subfiles
		ok = mergeRootFiles(output_file, output_files, merge_trees, parallel = False)

	if n_jobs > 1: # delete temporary files
		for f in output_files:
			os.remove(f)

	return ok

def pfunc(args):
	return mergeRootFiles(**args)

def mergeSubfolders(output_files, merge_trees):
	subfolders = defaultdict(list)
	for f in output_files:
		subfolders[ os.path.basename(os.path.dirname(f))].append( f )


	output_files = []


	n_cores = multiprocessing.cpu_count()
	n_subfolders = len( subfolders.keys() )
	n_parallel_subfolders = max(1, n_cores/6)
	pool = multiprocessing.Pool( n_parallel_subfolders )

	args = []
	for subname, subfiles in subfolders.iteritems():
		outfile = os.path.join( os.path.dirname( subfiles[0] ), "{0}.root".format(subname) )
		output_files.append( outfile )
		args.append( {'output_file': outfile, "input_files": subfiles, "merge_trees": merge_trees,
		              "n_jobs_max": max( 1, multiprocessing.cpu_count() / n_parallel_subfolders  )} )

	oks = pool.map(pfunc, iterable = args)
	ok = not (False in oks );


	return output_files if ok else [];


def filterExcludes( fileslist, excludes_file):
	'''
	Removes files from args which are listed in the excludes_file
	'''
	out_fileslist = []

	exclude_filelist = []
	if excludes_file:
		if os.path.isfile( excludes_file ):
			with open( excludes_file ) as fin:
				for line in fin:
					line = line.rstrip('\n').strip().strip("'").strip('"')
					if line:
						exclude_filelist.append( line );

			print exclude_filelist

			for filename in fileslist:
				filename_real = os.path.realpath( os.path.abspath(filename) )
				print filename_real
				print [ exclude in filename_real for exclude in exclude_filelist ]
				if not True in [ exclude in filename_real for exclude in exclude_filelist ]:
					out_fileslist.append( filename )
		else:
			print "Can not open excludes file '{0}'.".fromat(excludes_file)
			exit(1)
	else:
		out_fileslist = copy.copy( fileslist )
	return out_fileslist;


def splitFilesToJobs(all_files, n_files_per_job, doNotMixRuns):
	# group input files to jobs
	# do no nix files from different folders in one job
	grouped_inputfiles = defaultdict(list)
	for arg in all_files:
		dir = os.path.dirname( os.path.abspath(arg))
		grouped_inputfiles[dir].append( arg )

	if doNotMixRuns:
		grouped_inputfiles_byrun = defaultdict(list)
		for g, files in grouped_inputfiles.iteritems():
			for f in files:
				run = getRunNumber(f)
				grouped_inputfiles_byrun[(g,run)].append(f)
		grouped_inputfiles = grouped_inputfiles_byrun;

	files = []
	for group in sorted(grouped_inputfiles.keys()):
		group_files = grouped_inputfiles[group]
		n_files = len(group_files);
		n_files_job = min(n_files_per_job, n_files);
		n_jobs = n_files / n_files_job;
		q = n_files / n_jobs
		m = n_files % n_jobs
		files += [ group_files[i*q + min(m,i): (i+1)*q + min(m, i+1)] for i in xrange(n_jobs)]
	return files

def main():
	optparser = OptionParser( usage="Usage:%prog <args> [<options>]", description = program_description );
	optparser.add_option('-b', '--batchelorconfigfile', dest='batchelorconfigfile', action='store', type='str', default="", help="Config file for batchelor.")
	optparser.add_option('-c', '--configfile', dest='configfile', action='store', type='str', default="", help="Config file for antok.")
	optparser.add_option('-o', '--outfile', dest='outfile', action='store', type='str', default="hist.root",
	                     help="Merge files to the given output file. If given, all sub-files will be stored in a folder named according to the output file.")
	optparser.add_option('-t', '--merge-trees', dest='merge_trees', action='store_true', help="Merge not only histograms but also all trees in one file.")
	optparser.add_option('-l', '--local', dest='local', action='store_true', help="Run local.")
	optparser.add_option('-n', '--n-files-per-job', dest='n_files_job', action='store', default=1, help="Number of input files per parallel job [default: %default].")
	optparser.add_option('-j', '--n-jobs', dest='n_jobs', action='store', type='int', default=0, help="Number of parallel jobs [default: %default].")
	optparser.add_option('-s', '--subfolders', dest='subfolders' , action='store_true', help="Distribute output files to different folders, according to the input folders.")
	optparser.add_option('',   '--not-mix-runs', dest='not_mix_runs', action='store_true', help="Do not mix phast files from different runs in one job (base on filename).")
	optparser.add_option('-e', '--excludes-file', dest='excludes_file' , action='store', default="", help='''File which contains a list of files which should be excluded from the processing (One line for each excluded file / Real paths have to be used ). Also full folders can be excluded by giving the folder path.''')
	optparser.add_option('',   '--memory', dest='memory', default=None, help="Set memory requirements of jobs")

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

	outfilePath = os.path.realpath(options.outfile)
	if options.outfile == 'hist.root':
		outfolder = os.path.dirname(outfilePath)
	else:
		outfolder = os.path.join(os.path.dirname(outfilePath), os.path.splitext(os.path.basename(outfilePath))[0])
	options.outfile = outfilePath

	if options.batchelorconfigfile and not os.path.isfile( options.batchelorconfigfile ):
		print "Batchelor configfile '{0}' not found!".format(options.batchelorconfigfile)
		print optparser.usage
		exit( 100 );

	if options.n_jobs != 0 and options.n_files_job != 1:
		print "Cannot set number of files per job and number of jobs at the same time"
		exit( 100 )

	args = filterExcludes( args, options.excludes_file)

	options.n_files_job = min( int(options.n_files_job), len(args))
	if options.n_jobs != 0:
		options.n_files_job = max( int( float(len(args)) / options.n_jobs ), 1)

	files = splitFilesToJobs(args, options.n_files_job, doNotMixRuns = options.not_mix_runs)

	defaultBatchelorConfig = os.path.expanduser( "~/.batchelorrc");
	if not os.path.isfile( defaultBatchelorConfig ):
		print "Could not find default config file '{0}' for batchelor!".format(defaultBatchelorConfig)
		print optparser.usage
		exit( 100 );

	antok = os.path.join(os.path.dirname(os.path.realpath(__file__)), "treereader" )
	handler = batchelor.BatchelorHandler(configfile=defaultBatchelorConfig if not options.batchelorconfigfile else options.batchelorconfigfile,
	                                     systemOverride="local" if options.local else "",
	                                     memory=options.memory)

	log_files = [];
	out_files = [];

	log_dir = os.path.join( outfolder, 'log' )
	if not os.path.isdir(log_dir):
		os.makedirs( log_dir )
		if not options.local:
			# wait 3 seconds such that the log folder is appearing
			# I know it is stupid!
			time.sleep(3)


	# make my own copy of the config file and use it to be protected from changes during the execution
	# load and store it to handle !include statements in the config file
	local_configfile = os.path.join( log_dir, os.path.basename(options.configfile) )
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


	for i_job, in_files in enumerate(files):
		try:
			run_first_file = getRunNumber( in_files[0] )
			run_last_file = getRunNumber( in_files[-1] )
			slot = getSlot( in_files[0] )
		except Exception:
			handler.shutdown();
			exit(1)

		out_folder_root = outfolder
		if options.subfolders:
			out_folder_root = os.path.join( out_folder_root, os.path.basename( os.path.dirname( in_files[0] ) ) );
		if not os.path.isdir(out_folder_root):
			os.makedirs( out_folder_root );
		out_file = os.path.join( out_folder_root, "mDST-{0}-{1}-{2}.root.{3:03d}".format(run_first_file, run_last_file, slot, i_job) )
		out_files.append( out_file );

		cmd = "echo \"Start: $(date)\""
		in_filenames = ""
		for i_infile, in_file in enumerate(in_files):
			if not os.path.isfile( in_file ):
				print "File '{0}' not found!".format(in_file);
				handler.shutdown();
				exit(1)
			in_filenames += " '{0}'".format( in_file );
		cmd += " && {antok} {in_files} {out_file} {config_file}".format( antok = antok,
		                                                                out_file = out_file,
    	                                                                in_files = in_filenames,
    	                                                                config_file = options.configfile )



		cmd += " && echo \"STATUS: OK\" || echo \"STATUS: ERROR\""

		log_file = os.path.join( log_dir , os.path.basename(out_file));
		log_files.append( log_file )
		if os.path.isfile(log_file):
			with open(log_file + '.old', 'a') as fout:
				with open(log_file) as fin:
					fout.write(fin.read())
			os.remove(log_file)

		if os.path.isfile( out_file ) :
			print "File '{0}' already exists!".format(out_file);
			handler.shutdown();
			exit(1)


		print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
		print "Process files: "
		for f in in_files:
			print "\t'{0}'".format(f)
		print cmd
		print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
		handler.submitJob( cmd, log_file, wd = os.getcwd(), jobName='Antok' );




	try:
		handler.wait(2)
	except batchelor.CancelException:
		handler.shutdown();
		exit(1)
	handler.shutdown();

	print "==========================================================================="
	print "========================  FINISHED  ======================================="
	print "==========================================================================="

	errors = [];
	for l in log_files:
		for i in xrange(10):
			if os.path.exists(l):
				break
			time.sleep(1)
		if os.path.exists(l):
			with open(l) as fin:
				log_content = fin.read()
				if not "STATUS: OK\n" in log_content:
					errors.append(l)
					print "***************************************************************************"
					print "ERROR in this process: "
					print "***************************************************************************"
					print log_content
					print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
		else:
			errors.append(l)
			print "***************************************************************************"
			print "ERROR in this process: "
			print "***************************************************************************"
			print "Can not open logfile"
			print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"


	if options.outfile:
		if errors:
			options.outfile = os.path.join(os.path.dirname(options.outfile), 'error_{0}'.format(os.path.basename(options.outfile)))
		print "==========================================================================="
		print "========================  MERGING  ========================================"
		print "==========================================================================="

		if options.subfolders:
			files_to_merge = mergeSubfolders(out_files, merge_trees = options.merge_trees)
			mergeRootFiles(options.outfile, files_to_merge, merge_trees = options.merge_trees, parallel = False)
		else:
			mergeRootFiles(options.outfile, out_files, merge_trees = options.merge_trees)




#	 if options.clear and not errors:
#		 sp.check_call( "rm -f {0}".format(" ".join(out_files)), shell=True )

	for error in errors:
		print "***************************************************************************"
		print "ERROR in logfile '{0}' ".format(error)
		print "***************************************************************************"



	exit(0)


if __name__ == '__main__':
	main()




