#!/usr/bin/env python3
# coding: utf-8

import argparse
import coloredlogs
import numpy as np
import os
import sys
import uproot
import verboselogs

from collections import namedtuple, Counter
from collections.abc import Iterable


Pair = namedtuple('Pair', 'orig new')


def compareHists(TKey):
	# get histograms
	keyName = TKey.orig.object_path
	className = TKey.orig.classname()
	fileName = Pair(TKey.orig.file.file_path,
	                TKey.new.file.file_path)
	hist = Pair(TKey.orig.get(),
	            TKey.new.get())
	# compare number of entries
	entries = Pair(hist.orig.member('fEntries'),
	               hist.new.member('fEntries'))
	if entries.orig != entries.new:
		logger.error(f"histograms for key '{keyName}' have different number of entries: "
		             f"'{entries.orig}' in file '{fileName.orig}' vs. "
		             f"'{entries.new}' in file '{fileName.new}'")
		return False
	# compare binning
	axes = []
	if className.startswith('TH1'):
		axes = ['x']
	elif className.startswith('TH2'):
		axes = ['x', 'y']
	elif className.startswith('TH3'):
		axes = ['x', 'y', 'z']
	for axis in axes:
		edges = Pair(hist.orig.axis(axis).edges(),
		             hist.new.axis(axis).edges())
		if not np.array_equal(edges.orig, edges.new):
			logger.error(f"histograms for key '{keyName}' have different binning of {axis} axis")
			logger.verbose(f"     '{edges.orig!r}' in file '{fileName.orig}'")
			logger.verbose(f" vs. '{edges.new!r}' in file '{fileName.new}'")
	# compare histogram contents
	histContents = Pair([hist.orig.values(), hist.orig.errors()],
	                    [hist.new.values(),  hist.new.errors()])
	if not np.array_equal(histContents.orig, histContents.new):
		logger.error(f"histograms for key '{keyName}' have different contents")
		logger.verbose(f" orig - new values = '{histContents.orig[0] - histContents.new[0]!r}'")
		logger.verbose(f" orig - new errors = '{histContents.orig[1] - histContents.new[1]!r}'")
		return False
	return True


def compareLeafs(branch):
	branchName = branch.orig.name
	keyName = branch.orig.tree.object_path
	# compare leaf type
	if branch.orig.typename != branch.new.typename:
		logger.error(f"branches '{branchName}' in trees '{keyName}' have different types: "
		             f"'{branch.orig.typename}' vs. '{branch.new.typename}'")
		return False
	# compare leaf values
	for itOrig, itNew in zip(branch.orig.iterate(), branch.new.iterate()):
		assert len(itOrig) == len(itNew), (f"branches '{branchName}' in trees '{keyName}' have different number of entries: "
		                                   f"'{len(itOrig)}' vs. '{len(itNew)}'")
		# np.array_equal(np.array(itOrig), np.array(itNew)) does not work for leafs with variable sized arrays
		# hence we have to loop over the entries in Python making things unfortunately much slower
		for entry, valOrig in enumerate(itOrig):
			if not np.array_equal(np.array(valOrig[branchName]), np.array(itNew[entry][branchName]), equal_nan=True):
				logger.error(f"values in branches '{branchName}' of trees '{keyName}' differ: "
				             f"{itOrig} vs. {itNew}")
				return False
	return True


def compareTrees(TKey):
	result = True
	# get trees
	keyName = TKey.orig.object_path
	fileName = Pair(TKey.orig.file.file_path,
	                TKey.new.file.file_path)
	tree = Pair(TKey.orig.get(),
	            TKey.new.get())
	# compare number of entries
	if tree.orig.num_entries == tree.new.num_entries:
		logger.success(f"trees '{keyName}' have the same number of entries")
	else:
		logger.error(f"trees '{keyName}' have different number of entries: "
		             f"{tree.orig.num_entries} vs. {tree.new.num_entries}")
		return False
	# compare branches
	branchKeys = Pair(sorted(tree.orig.keys()),
	                  sorted(tree.new.keys()))
	branchKeysSet = Pair(set(branchKeys.orig),
	                     set(branchKeys.new))
	# check that branches in each tree are unique
	if len(branchKeysSet.orig) == len(branchKeys.orig):
		logger.success(f"branches in tree '{keyName}' in file '{fileName.orig}' are unique")
	else:
		logger.warning(f"branches '{branchKeys.orig}' in tree '{keyName}' in file '{fileName.orig}' are not unique")
	if len(branchKeysSet.new) == len(branchKeys.new):
		logger.success(f"branches in tree '{keyName}' in file '{fileName.new}' are unique")
	else:
		logger.warning(f"branches '{branchKeys.new}' in tree '{keyName}' in file '{fileName.new}' are not unique")
	# check tast the two tress have the same branches
	if branchKeysSet.orig == branchKeysSet.new:
		logger.success(f"trees for key '{keyName}' have identical branches")
	else:
		result = False
		logger.error(f"trees for key '{keyName}' have different branches:")
		onlyIn = Pair(branchKeysSet.orig - branchKeysSet.new,
		              branchKeysSet.new  - branchKeysSet.orig)
		if onlyIn.orig:
			logger.error(f"    the following branches are found only in tree '{keyName}' in file '{fileName.orig}': {onlyIn.orig}")
		if onlyIn.new:
			logger.error(f"    the following branches are found only in tree '{keyName}' in file '{fileName.new}': {onlyIn.new}")
	# check leafs
	count = Counter()
	for branchKey in branchKeys.orig:
		count['BranchesTotal'] += 1
		logger.info(f"checking branch '{branchKey}' in tree '{keyName}' [{count['BranchesTotal']} of {len(branchKeys.orig)}]")
		try:
			branchNew = tree.new[branchKey]
		except uproot.KeyInFileError:
			result = False
			logger.error(f"cannot find branch '{branchKey}' in tree '{keyName}' file '{fileName.new}'")
			continue
		count['BranchesBoth'] += 1
		branch = Pair(tree.orig[branchKey], branchNew)
		if compareLeafs(branch):
			count['BranchesSuccess'] += 1
			logger.success(f"branch '{branchKey}' in trees for key '{keyName}' is identical in both files")
		else:
			result = False
		# break

	# print summary for tree
	logger.info(f"Comparison summary after checking {len(branchKeys.orig)} branches in trees for key '{keyName}':")
	if count['BranchesSuccess'] == count['BranchesBoth']:
		logger.success(f"all of the {count['BranchesBoth']} branches that are found in both trees for key '{keyName}' are identical")
	else:
		logger.warning(f"{count['BranchesBoth'] - count['BranchesSuccess']} branches that are found in both trees for key '{keyName}' differ")

	return result


def compareFiles(fileName=Pair('orig.root', 'new.root')):
	result = True
	logger.info(f"comparing files '{fileName.orig}' and '{fileName.new}'")
	file = Pair(uproot.open(fileName.orig),
	            uproot.open(fileName.new))
	keys = Pair(sorted(file.orig.keys(recursive=True)),
	            sorted(file.new.keys(recursive=True)))

	# check that keys in each file are unique
	keySet = Pair(set(keys.orig),
	              set(keys.new))
	assert len(keySet.orig) == len(keys.orig), f"keys in file '{fileName.orig}' are not unique"
	assert len(keySet.new) == len(keys.new), f"keys in file '{fileName.new}' are not unique"

	# check each key
	count = Counter()
	for keyName in keys.orig:
		count['KeysTotal'] += 1
		logger.info(f"checking key '{keyName}' [{count['KeysTotal']} of {len(keys.orig)}]")
		# check that object with same key exists in new file
		try:
			TKeyNew = file.new.key(keyName)
		except uproot.KeyInFileError:
			logger.error(f"cannot find key '{keyName}' in file '{fileName.new}'")
			continue
		TKey = Pair(file.orig.key(keyName), TKeyNew)
		# check that objects are of the same class
		className = Pair(TKey.orig.classname(), TKey.new.classname())
		if className.orig != className.new:
			logger.error(f"classes for key '{keyName}' differ: "
			             f"'{className.orig}' in file '{fileName.orig}' vs. "
			             f"'{className.new}' in file '{fileName.new}'")
			continue
		# check histograms
		elif className.orig.startswith(('TH1', 'TH2', 'TH3')):
			count['HistsTotal'] += 1
			if compareHists(TKey):
				count['HistsSuccess'] += 1
				logger.success(f"histograms for key '{keyName}' are identical in both files")
		# check trees
		elif className.orig == 'TTree':
			count['TreesTotal'] += 1
			if compareTrees(TKey):
				count['TreesSuccess'] += 1
				logger.success(f"trees for key '{keyName}' are identical in both files")
		else:
			logger.notice(f"cannot compare object of class '{className.orig}'")

	# print summary
	logger.info(f"Comparison summary after checking {len(keys.orig)} keys:")
	if keySet.orig == keySet.new:
		logger.success(f"files '{fileName.orig}' and '{fileName.new}' have identical keys")
	else:
		result = False
		logger.error('files have different keys:')
		onlyIn = Pair(keySet.orig - keySet.new,
		              keySet.new  - keySet.orig)
		if onlyIn.orig:
			logger.error(f"    the following keys are found only in file '{fileName.orig}': {onlyIn.orig}")
		if onlyIn.new:
			logger.error(f"    the following keys are found only in file '{fileName.new}': {onlyIn.new}")

	if count['HistsSuccess'] == count['HistsTotal']:
		logger.success(f"all {count['HistsTotal']} histograms that are found in both files are identical")
	else:
		result = False
		logger.error(f"{count['HistsTotal'] - count['HistsSuccess']} histograms that are found in both files differ")
	if count['TreesSuccess'] == count['TreesTotal']:
		logger.success(f"all {count['TreesTotal']} trees that are found in both files are identical")
	else:
		result = False
		logger.error(f"{count['TreesTotal'] - count['TreesSuccess']} trees that are found in both files differ")

	file.orig.close()
	file.new.close()
	return result


if __name__ == '__main__':
	# setup logger
	# see https://pypi.org/project/coloredlogs/
	logger = verboselogs.VerboseLogger(__name__)
	coloredlogs.DEFAULT_FIELD_STYLES['levelname'] = {'bold': True}
	coloredlogs.DEFAULT_LEVEL_STYLES['success'] = {'color': 'green'}
	coloredlogs.install(level='INFO', logger=logger, fmt='%(levelname)s %(message)s')
	# logger.spam("this is a spam message")
	# logger.debug("this is a debugging message")
	# logger.verbose("this is a verbose message")
	# logger.info("this is an informational message")
	# logger.notice("this is a notice message")
	# logger.warning("this is a warning message")
	# logger.success("this is a success message")
	# logger.error("this is an error message")
	# logger.critical("this is a critical message")

	# parse command line arguments
	parser = argparse.ArgumentParser(
		description="Compares the contents of histograms and trees of two ROOT files."
	)
	parser.add_argument("filePathA", type=str, help="first input file")
	parser.add_argument("filePathB", type=str, help="second input file")
	args = parser.parse_args()

	if compareFiles(Pair(os.path.abspath(args.filePathA),
	                     os.path.abspath(args.filePathB))):
		logger.success(f"files '{args.filePathA}' and '{args.filePathB}' are identical")
		sys.exit(os.EX_OK)
	else:
		logger.error(f"files '{args.filePathA}' and '{args.filePathB}' differ")
		sys.exit(os.EX_DATAERR)
