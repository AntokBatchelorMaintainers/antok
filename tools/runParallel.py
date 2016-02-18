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

# root includes

# own includes
import batchelor


def main():
    optparser = OptionParser( usage="Usage:%prog <args> [<options>]", description = program_description );
    optparser.add_option('-c', '--configfile', dest='configfile', action='store', type='str', default="", help="Config file for antok.")
    optparser.add_option('-o', '--outfile', dest='outfile', action='store', type='str', default="hist.root", help="Merge files to the given output file.")
    optparser.add_option('-l', '--local', dest='local', action='store_true', help="Run local.")
    optparser.add_option('-n', '--n-files-per-job', dest='n_files_job', action='store', default=1, help="Number of input files per parallel job [default: %default]")
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
    n_jobs = len(args) / int(options.n_files_job)

    antok = os.path.join( os.path.dirname( __file__), "treereader" )
    
    q = len(args) / n_jobs
    m = len(args) % n_jobs
    files = [ args[i*q + min(m,i): (i+1)*q + min(m, i+1)] for i in xrange(n_jobs)]


    handler = batchelor.BatchelorHandler(configfile="~/.batchelorrc", 
                                         systemOverride="local" if options.local else "", 
                                         n_threads=3)


    log_files = [];
    out_files = [];

    for i_job, in_files in enumerate(files):
        out_file = os.path.join( os.path.dirname( options.outfile ), "{0}.root.{1:03d}".format(os.path.basename(in_files[0]).split('.root')[0], i_job) )
        out_files.append( out_file );
        cmd = "echo \"Start: $(date)\""
        local_out_files= []
        for i_infile, in_file in enumerate(in_files):
            if not os.path.isfile( in_file ):
                print "File '{0}' not found!".format(in_file);
                handler.shutdown();
                exit(1)
            local_out_file = out_file + '.{0:03d}'.format(i_infile)
            local_out_files.append( local_out_file )
            cmd += " && {antok} {in_file} {out_file} {config_file}".format( antok = antok, 
                                                                           out_file = local_out_file, 
                                                                           in_file = in_file, 
                                                                           config_file = options.configfile ) 
        
        
        if len(local_out_files) > 1:
            cmd += " && hadd {0} {1} >> /dev/null".format(out_file, " ".join(local_out_files))
        else:
            cmd += " && mv {0} {1}".format( local_out_files[0], out_file )
        
        if options.clear:
            cmd += " && rm -f " + " ".join(local_out_files)

        cmd += " && echo \"STATUS: OK\" || echo \"STATUS: ERROR\""

        log_file = os.path.join( os.path.dirname(out_file), "log-{0}".format(os.path.basename(out_file)))
        log_files.append( log_file )

        if os.path.isfile( out_file ):
            print "File '{0}' already exists!".format(out_file);
            handler.shutdown();
            exit(1)


        print "Process file", in_file
        print cmd
        handler.submitJob( cmd, log_file, wd = os.getcwd() );




    handler.wait(2)
    handler.shutdown();

    print "==========================================================================="
    print "========================  FINISHED  ======================================="
    print "==========================================================================="

    errors = [];
    for l in log_files:
        with open(l) as fin:
            log_content = fin.read()
            if not "STATUS: OK" in log_content:
                errors.append(log_file)
                print "***************************************************************************"
                print "ERROR in this process: "
                print "***************************************************************************"
                print log_content
                print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"

    if options.outfile:
        print "==========================================================================="
        print "========================  MERGING  ========================================"
        print "==========================================================================="

        sp.check_call( "hadd -T {outfile} {out_files}".format( outfile = options.outfile, out_files = " ".join(out_files) ), shell=True, stdout=sp.PIPE, stderr = sp.STDOUT  );
        
#     if options.clear and not errors:
#         sp.check_call( "rm -f {0}".format(" ".join(out_files)), shell=True )
        
    for error in errors:
        print "***************************************************************************"
        print "ERROR in logfile '{0}' ".format(error)
        print "***************************************************************************"
        


    exit(0)


if __name__ == '__main__':
    main()




