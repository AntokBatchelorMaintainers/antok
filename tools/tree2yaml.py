#!/usr/bin/env python
# coding: utf-8
'''
Created on Fri 02 Oct 2020 14:13:12 AM CET


@author: Christian Dreisbach
'''

program_description = '''
    <description>
'''

from optparse import OptionParser
import ROOT
import os
import yaml
import time

def createDataTypeNodes(infile, treename):
	nodes = dict()
	rootFile = ROOT.TFile(infile, "READ")
	tree = rootFile.Get(treename)
	branches = tree.GetListOfBranches()
	print "Found {0} branches".format(len(branches))
	for branch in branches:
		branchName = branch.GetName()
		branchTitle = branch.GetTitle().strip()
		print " Processing branch '{0}'".format(branchName)
		className = branch.GetClassName()
		if className is not "":
			dataType = "std::"+className
		else:
			type = branchTitle[-1]
			if type is "D":
				dataType = "double"
			elif type is "I":
				dataType = "int"
			else:
				dataType = "UNKNOWN"
		if dataType not in nodes:
			nodes[dataType] = dict()
			print " -> Created new data type node '{0}'".format(dataType)

		if dataType.find("std::vector") > -1:
			nodes[dataType]["- &"+branchName] = branchName
		else:
			nodes[dataType]["- &" + branchName] = branchName
	return nodes

def main():
	optparser = OptionParser( usage="Usage:%prog <args> [<options>]", description = program_description )
	optparser.add_option('-i', '--infile' , dest='infile' , action='store', type='str'                       , help="Input ROOT file")
	optparser.add_option('-o', '--outfile', dest='outfile', action='store', type='str', default="output.yaml", help="YAML outputfile")
	optparser.add_option('-t', '--tree'   , dest='tree'   , action='store', type='str',                        help="Tree to read from input file.")

	( options, args ) = optparser.parse_args()

	if not options.tree:
		print "No tree name given"
		print optparser.usage
		exit( 100 )

	if not options.infile:
		print "No input file name given"
		print optparser.usage
		exit( 100 )

	if not options.outfile:
		print "Using default output file '{0}'".format( options.outfile)

	if not os.path.isfile( options.infile ):
		print "Input file '{0}' not found!".format(options.infile)
		print optparser.usage
		exit( 100 )

	options.outfile = os.path.realpath(options.outfile)
	options.infile  = os.path.realpath(options.infile)
	print "Input file set to:  '{0}'".format(options.infile)
	print "Tree set to:        '{0}'".format(options.tree)
	print "Output file set to: '{0}'".format(options.outfile)

	# - Create structure -
	nodes = []
	treeNameNode     = dict( TreeName = options.tree)
	dataTypesNode    = createDataTypeNodes(options.infile, options.tree)
	onePerEventNode  = dict( onePerEvent  = dataTypesNode)
	treeBranchesNode = dict( TreeBranches = onePerEventNode )

	nodes.append(treeNameNode)
	nodes.append(treeBranchesNode)

	# - Write to tmp yaml file -
	tempFileName = "tmp_{0}.yaml".format(int(time.time()))
	with open(tempFileName, "w") as tempFile:
		for node in nodes:
			yaml.dump(node, tempFile, default_flow_style=False)

	# - Refactor tmp yaml file into ANTOK "format" -
	with open(tempFileName, "r") as tempFile:
		with open(options.outfile, "w") as outfile:
			lines = tempFile.readlines()
			longestLine = 0
			for line in lines:
				if line.find("- &") > -1:
					line = line.split(":")
					if len(line[0]) > longestLine:
						longestLine = len(line[0])
			for line in lines:
				line = line.replace("'", "")
				line = line.replace('"', "")
				if line.find("- &") > -1:
					data = line.split(":")
					formatString = "{0[0]:<10} {0[1]:>"+str(longestLine)+"}"
					line =  formatString.format(data)
				if line.find(":") > -1:
					outfile.write("\n")
				outfile.write(line)
	os.remove(tempFileName)

	print "Finished! Data saved in: '{0}'".format(options.outfile)

	exit(0)

if __name__ == '__main__':
	main()
