#!/usr/bin/python
#

import os
import sys
import glob
from optparse import OptionParser

from pyminc.volumes.factory import *
from numpy import *
#import tempfile
#import functools


def get_neighbour_volume(inputfile, outputfile, order=1, lowerthr=0.5, upperthr=inf, verbose=False):

	if not(os.path.exists(inputfile)):
		raise IOError('File does not exist: %s' % inputfile)

	if (lowerthr > upperthr):
		raise ValueError('upper threshold must be greater than lower threshold')

	invol = volumeFromFile(inputfile)
	if verbose: print 'calculating neighbours for:', inputfile

	outvol = volumeFromInstance(invol, outputfile, volumeType='int')
	invol.data = ((invol.data>=lowerthr) & (invol.data<=upperthr)).astype(int)

	# order of dimensions as written on disk, see mincinfo
	for i in range(1,order+1):
		if verbose: print 'order', i
		if i==1:
			outvol.data[1:-1,1:-1,1:-1] = (
				invol.data[0:-2,1:-1,1:-1] +
				invol.data[1:-1,0:-2,1:-1] +
				invol.data[1:-1,1:-1,0:-2] +
				invol.data[2:  ,1:-1,1:-1] +
				invol.data[1:-1,2:  ,1:-1] +
				invol.data[1:-1,1:-1,2:  ] +
				invol.data[0:-2,0:-2,1:-1] +  #
				invol.data[0:-2,1:-1,0:-2] +
				invol.data[1:-1,0:-2,0:-2] +
				invol.data[0:-2,0:-2,0:-2] +
				invol.data[2 : ,2:  ,1:-1] +
				invol.data[2:  ,1:-1,2:  ] +
				invol.data[1:-1,2:  ,2:  ] +
				invol.data[2:  ,2:  ,2:  ] +
				invol.data[0:-2,0:-2,2:  ] +  #
				invol.data[0:-2,2:  ,0:-2] +
				invol.data[2:  ,0:-2,0:-2] +
				invol.data[0:-2,0:-2,0:-2] +
				invol.data[2 : ,0:-2,1:-1] +
				invol.data[0:-2,1:-1,0:-2] +
				invol.data[1:-1,0:-2,0:-2] +
				invol.data[0:-2,0:-2,0:-2] +
				invol.data[0:-2,2:  ,1:-1] +  #
				invol.data[0:-2,1:-1,2:  ] +
				invol.data[2:  ,0:-2,1:-1] +
				invol.data[1:-1,0:-2,2:  ] +
				invol.data[2:  ,1:-1,0:-2] +
				invol.data[1:-1,2:  ,0:-2]
			)
			outvol.data = outvol.data * invol.data
		elif i>1:
			outvol.data[1:-1,1:-1,1:-1] = (
				outvol.data[0:-2,1:-1,1:-1] +
				outvol.data[1:-1,0:-2,1:-1] +
				outvol.data[1:-1,1:-1,0:-2] +
				outvol.data[2:  ,1:-1,1:-1] +
				outvol.data[1:-1,2:  ,1:-1] +
				outvol.data[1:-1,1:-1,2:  ] +
				outvol.data[0:-2,0:-2,1:-1] +  #
				outvol.data[0:-2,1:-1,0:-2] +
				outvol.data[1:-1,0:-2,0:-2] +
				outvol.data[0:-2,0:-2,0:-2] +
				outvol.data[2 : ,2:  ,1:-1] +
				outvol.data[2:  ,1:-1,2:  ] +
				outvol.data[1:-1,2:  ,2:  ] +
				outvol.data[2:  ,2:  ,2:  ] +
				outvol.data[0:-2,0:-2,2:  ] +  #
				outvol.data[0:-2,2:  ,0:-2] +
				outvol.data[2:  ,0:-2,0:-2] +
				outvol.data[0:-2,0:-2,0:-2] +
				outvol.data[2 : ,0:-2,1:-1] +
				outvol.data[0:-2,1:-1,0:-2] +
				outvol.data[1:-1,0:-2,0:-2] +
				outvol.data[0:-2,0:-2,0:-2] +
				outvol.data[0:-2,2:  ,1:-1] +  #
				outvol.data[0:-2,1:-1,2:  ] +
				outvol.data[2:  ,0:-2,1:-1] +
				outvol.data[1:-1,0:-2,2:  ] +
				outvol.data[2:  ,1:-1,0:-2] +
				outvol.data[1:-1,2:  ,0:-2]
			)
			outvol.data = outvol.data * invol.data
	outvol.writeFile()
	outvol.closeVolume()
	invol.closeVolume()
	if verbose: print 'saved output to:', outputfile


if __name__ == "__main__":
	usage = "usage: %prog [options] -i MAGETDIR -r REFERNECE AVERAGES.."
	parser = OptionParser(usage)
	parser.add_option("-n", "--order", help='neighbourhood order (>=1)', type="int", default=0)
	parser.add_option("-o", "--output", help='output filename', type="string", default='')
	parser.add_option("--lowerthreshold", help='threshold input with this lower threshold', type="float", default=0.5)
	parser.add_option("--upperthreshold", help='threshold input with this upper threshold', type="float", default=inf)
	parser.add_option("-v", "--verbose", help='verbose output', action='store_true', default=False)
	(options, args) = parser.parse_args()

	if len(args) < 1:
		parser.print_usage()
		sys.exit(1)

	input = args[0]

	if options.output:
		output = options.output
	else:
		output = '%s_connectivity_order%i%s' % (os.path.splitext(input)[0],options.order,os.path.splitext(input)[1])
		
	get_neighbour_volume(input,output,order=options.order,lowerthr=options.lowerthreshold,upperthr=options.upperthreshold,verbose=options.verbose)
	
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



