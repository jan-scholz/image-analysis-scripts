#!/usr/bin/python
#
# weighted_vote.py -r bdnf_reg_all_nlin-3.mnc -i bdnf-mincants atlantes/*average.mnc -v
# /projects/mice/share/arch/linux-x86_64-eglibc2_11_1/lib/python2.6/site-packages/pyminc/volumes/*py

import os
import sys
import glob
from optparse import OptionParser

from pyminc.volumes.factory import *
from numpy import *
import tempfile
import functools


def get_basename(filepath):
	"""returns file basename with extension removed, e.g. /foo/bar.txt -> bar"""
	return os.path.basename(os.path.splitext(filepath)[0])

def get_MAGeTfiles(reference, averages, magetdir, suffix):
	"""returns files created by MAGeT.py by aligning averages to a target reference
	reference --  input file to MAGeT
	averages  --  files in atlas directory fed into MAGeT
	magetdir  --  the directory that MAGeT output its files too
	suffix    --  either 'labels' or 'resampled'"""
	outfiles = []
	for f in averages:
		f = glob.glob(os.path.join(magetdir,get_basename(reference),suffix,get_basename(reference) + '*_to_' + get_basename(f) + '*' + suffix + '.mnc'))
		if len(f)>1: sys.exit('multiple matches for input file ' + f)
		if not(os.path.isfile(f[0])): sys.exit('could not open MAGeT file ' + f[0])
		outfiles.append(f[0])
	return outfiles


def get_weights(reference, inputimages, verbose=False):

	# radius of neighbourhood cube, excluding centre voxel
	r = 2

	refvol = volumeFromFile(reference)

	for f in inputimages:
		if verbose: print 'reading file', f
		vol = volumeFromFile(f)

		# order of dimensions as written on disk, see mincinfo
		for a in range(r,vol.data.shape[0]-r):
			for b in range(r,vol.data.shape[1]-r):
				for c in range(r,vol.data.shape[2]-r):
					

	return 1
	#return weightfiles




if __name__ == "__main__":
	usage = "usage: %prog [options] -i MAGETDIR -r REFERNECE AVERAGES.."
	parser = OptionParser(usage)
	parser.add_option("-i", "--magetdir", help='directory containing output of MAGeT.py', type="string")
	parser.add_option("-r", "--reference", help='the reference that the averages got aligned to', type="string")
	parser.add_option("-m", "--mask", help='restrict vote to masked region', type="string")
	parser.add_option("-v", "--verbose", help='verbose output', action='store_true', default=False)
	(options, args) = parser.parse_args()

	if not(options.magetdir) or not(options.reference):
		parser.print_usage()
		sys.exit(1)

	averages = args

	labelfiles = get_MAGeTfiles(options.reference, averages, options.magetdir, 'labels')
	resampledfiles = get_MAGeTfiles(options.reference, averages, options.magetdir, 'resampled')

	#for l,f in zip(labelfiles,averages): print l, '\n', f, '\n\n'

	weightfiles = get_weights(options.reference,resampledfiles,options.verbose)


	
	# estimate SSD at each location in averages
	# create weight maps

	# use weight maps to do label fusion
	# output fused labels



# bdnf_reg_all_nlin-3_to_NRXN1a_version_1_average_base_0-resampled.mnc

#	print labelsuffix

	#source = args[0]
	#target = args[1]

	#sourcevol = volumeFromFile(source)
	#targetvol = volumeFromFile(target)

	#result = compute_xcorr(sourcevol, targetvol, options.mask)

	#print result



