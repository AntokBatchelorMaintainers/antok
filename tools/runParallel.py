#!/usr/bin/env python
# coding: utf-8
'''
Created on Sat 09 Jan 2016 12:05:35 AM CET


@author: Stefan Wallner
'''

program_description = ''' 
    <description>
'''

# std includes
from optparse import OptionParser
from collections import defaultdict
import os
import subprocess as sp
import shutil
import re

# root includes

# own includes
import batchelor


def getRunNumber(filename):
    '''
    @return: run number, extracted from mDST filename
    '''
    runnbr = None;
    parsed = re.findall("-([0-9]+)-[0-9]-[0-9]\.root", filename) 
    if(parsed):
        runnbr = int( parsed[0]) 
    else:
        msg = "Can not get run number from file name '{0}'".format( filename )
        print msg
        raise Exception(msg)
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
        msg = "Can not get slot from file name '{0}'".format( filename )
        print msg
        raise Exception(msg)
    return slot



def mergeSubfolders(output_files, merge_trees):
    subfolders = defaultdict(list)
    for f in output_files:
        subfolders[ os.path.basename(os.path.dirname(f))].append( f )
        
        
    processes = []
    output_files = []
    
    for subname, subfiles in subfolders.iteritems():
        outfile = os.path.join( os.path.dirname( subfiles[0] ), "{0}.root".format(subname) )
        output_files.append( outfile )
        processes.append( sp.Popen("hadd {opts} {outfile} {infiles} ".format( opts = "" if merge_trees else "-T",
                                                                              outfile = outfile, 
                                                                              infiles = "'" + "' '".join(subfiles) + "'" )
                                   ,shell=True, stdout=sp.PIPE, stderr = sp.STDOUT  )
                         )
        
    ok = True
    for p in processes:
        p.communicate()
        status = p.returncode
        if status != 0:
            print "Process exited with exit status", status
            ok = False;
    
    return output_files if ok else [];

def main():
    optparser = OptionParser( usage="Usage:%prog <args> [<options>]", description = program_description );
    optparser.add_option('-c', '--configfile', dest='configfile', action='store', type='str', default="", help="Config file for antok.")
    optparser.add_option('-o', '--outfile', dest='outfile', action='store', type='str', default="hist.root", help="Merge files to the given output file.")
    optparser.add_option('-t', '--merge-trees', dest='merge_trees', action='store_true', help="Merge not only histograms but also all trees in one file.")
    optparser.add_option('-l', '--local', dest='local', action='store_true', help="Run local.")
    optparser.add_option('-n', '--n-files-per-job', dest='n_files_job', action='store', default=1, help="Number of input files per parallel job [default: %default]")
    optparser.add_option('-s', '--subfolders', dest='subfolders' , action='store_true', help="Distribute output files to different folders, according to the input folders")
    optparser.add_option('', '--no-clear', dest='clear', action='store_false', default=True, help="Do not clear temporarily created histogram files after merging them to one" );

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
        
    options.n_files_job = min( int(options.n_files_job), len(args))
    
    # group input files to jobs
    # do no nix files from different folders in one job
    grouped_inputfiles = defaultdict(list)
    for arg in args:
        dir = os.path.dirname( os.path.realpath(arg))
        grouped_inputfiles[dir].append( arg )

    files = []
    for group in sorted(grouped_inputfiles.keys()):
        group_files = grouped_inputfiles[group]
        print group, group_files
        n_files = len(group_files);
        n_files_job = min(options.n_files_job, n_files);
        n_jobs = n_files / n_files_job;
        q = n_files / n_jobs
        m = n_files % n_jobs
        files += [ group_files[i*q + min(m,i): (i+1)*q + min(m, i+1)] for i in xrange(n_jobs)]
        

    antok = os.path.join( os.path.dirname( __file__), "treereader" )
    handler = batchelor.BatchelorHandler(configfile="~/.batchelorrc", 
                                         systemOverride="local" if options.local else "", 
                                         memory='2.5G',
                                         n_threads=3)


    log_files = [];
    out_files = [];
    
    log_dir = os.path.join( os.path.dirname(options.outfile), 'log' )
    if not os.path.isdir(log_dir):
        os.makedirs( log_dir )

    for i_job, in_files in enumerate(files):
        try:
            run_first_file = getRunNumber( in_files[0] )
            run_last_file = getRunNumber( in_files[-1] )
            slot = getSlot( in_files[0] )
        except Exception:
            handler.shutdown();
            exit(1)

        out_folder = os.path.dirname( options.outfile );
        if options.subfolders:
            out_folder = os.path.join( out_folder, os.path.basename( os.path.dirname( in_files[0] ) ) );
        if not os.path.isdir(out_folder):
            os.makedirs( out_folder );
        out_file = os.path.join( out_folder, "mDST-{0}-{1}-{2}.root.{3:03d}".format(run_first_file, run_last_file, slot, i_job) ) 
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

        if os.path.isfile( out_file ):
            print "File '{0}' already exists!".format(out_file);
            handler.shutdown();
            exit(1)


        print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
        print "Process files: "
        for f in in_files:
            print "\t'{0}'".format(f)
        print cmd
        print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
        handler.submitJob( cmd, log_file, wd = os.getcwd() );




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
        with open(l) as fin:
            log_content = fin.read()
            if not "STATUS: OK\n" in log_content:
                errors.append(l)
                print "***************************************************************************"
                print "ERROR in this process: "
                print "***************************************************************************"
                print log_content
                print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"

    if options.outfile:
        print "==========================================================================="
        print "========================  MERGING  ========================================"
        print "==========================================================================="

        files_to_merge = out_files
        if options.subfolders:
            files_to_merge = mergeSubfolders(out_files, merge_trees = options.merge_trees)

        sp.check_call( "hadd {merge} {outfile} {out_files}".format( merge = "" if options.merge_trees else "-T",  
                                                                    outfile = options.outfile, 
                                                                    out_files = " ".join(files_to_merge) ), 
                      shell=True, stdout=sp.PIPE, stderr = sp.STDOUT  );
            
        
#     if options.clear and not errors:
#         sp.check_call( "rm -f {0}".format(" ".join(out_files)), shell=True )
        
    for error in errors:
        print "***************************************************************************"
        print "ERROR in logfile '{0}' ".format(error)
        print "***************************************************************************"
        


    exit(0)


if __name__ == '__main__':
    main()




