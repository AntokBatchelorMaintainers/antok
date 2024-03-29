#!/usr/bin/env python

import compileall
import os
import sys


if __name__ == "__main__":

	indir = sys.argv[1]
	destdir = sys.argv[2]

	# legacy is required to maintain directory structure for python3, not compatible with python2
	try:
		compileall.compile_dir(dir=indir, ddir=destdir, legacy=True, quiet=True)
	except:
		compileall.compile_dir(dir=indir, ddir=destdir, quiet=True)

	os.system('rm -rf ' + destdir + '/*')
	os.system('cd ' + indir + '; for dir in `find * -type d`; do [[ "$dir" != *"__pycache__" ]] && mkdir ' + destdir + '/$dir; done; cd - > /dev/null')
	os.system('cd ' + indir + '; for pycfile in `find * -name "*.pyc"`; do mv $pycfile ' +  destdir + '/$pycfile; done; cd - > /dev/null')
