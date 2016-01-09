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
    optparser.add_option('-o', '--outfile', dest='outfile', action='store', type='str', default="", help="Merge files to the given output file.")

    ( options, args ) = optparser.parse_args();

    if not options.configfile:
        print "No configfile given"
        print options.usage
        exit( 100 );

    if not os.path.isfile( options.configfile ): 
        print "Config file '{0}' not found!".format(options.configfile)
        print options.usage
        exit( 100 );


    if options.outfile and os.path.isfile( options.outfile ): 
        print "Output file '{0}' exists found!".format(options.outfile)
        print options.usage
        exit( 100 );

    antok = os.path.join( os.path.dirname( __file__), "treereader" )


    handler = batchelor.BatchelorHandler(configfile="~/.batchelorrc", systemOverride="local")


    log_files = [];
    out_files = [];

    for in_file in args:
        if not os.path.isfile( in_file ):
            print "File '{0}' not found!".format(in_file);
            break;

        out_file = os.path.join( os.getcwd(), os.path.basename( in_file ) );
        out_files.append( out_file );

        log_file = os.path.join( os.getcwd(), os.path.splitext(os.path.basename( in_file ))[0] ) + ".log";
        log_files.append( log_file )

        if os.path.isfile( out_file ):
            print "File '{0}' already exists!".format(out_file);
            break;

        cmd = "{antok} {in_file} {out_file} {config_file}".format( antok = antok, out_file = out_file, in_file = in_file, config_file = options.configfile ) 

        print "Process file", in_file
        print cmd
        handler.submitJob( cmd, log_file, wd = os.getcwd() );




    handler.wait(2)
    handler.shutdown();

    print "==========================================================================="
    print "========================  FINISHED  ======================================="
    print "==========================================================================="

    for l in log_files:
        with open(l) as fin:
            print fin.read()
            print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"

    if options.outfile:
        print "==========================================================================="
        print "========================  MERGING  ========================================"
        print "==========================================================================="

        sp.check_call( "hadd {output} {out_files}".format( output = options.outfile, out_files = " ".join(out_files) ), shell=True  );



    exit(0)


if __name__ == '__main__':
    main()




