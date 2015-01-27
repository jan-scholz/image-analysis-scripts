#!/usr/bin/env python
#
# author: jan.scholz@mouseimaging.ca @ MICe, SickKids, Toronto, Canada
#
# ssh -Y topolina
# . /home/matthijs/./scripts/topolina_environment
# cd /projects/mice/jscholz/tmp/def/simulation
# simulate_sublabels.py -i foo002.mnc -o foo -v 500

from optparse import OptionParser,OptionGroup
from pyminc.volumes.factory import *
from scipy import ndimage,stats
import numpy as np
from os import path
import sys
import csv
import math


def simulateData(labelfilename,outcsv,replications,seed,verbose=False):
# creates sample data based on labelfile and model
#
# outcsv         csv containing filenames and group association
# replications   dict of groupname and nr of subjects: {'controls':10,'mutants':10}
# modelfiles     simulate hierachy based on text file, named according to group
#                  label,depth_for_simulation,mu,mu_sd,sd_gamma
#                  from top to bottom
#                  mu's that have a parent get added to the the value from their parent

# set seed
# read in labels file
# sample top+1 level
#   and work the way down the hierachy
# assing values to array
# write to disk
# create csv with filename & group columns
	pass





def checkFiles(table, filenamecolumn='filename', verbose=False):
    # check that table contains column: filename
    filenames = sorted(list(set([e[filenamecolumn] for e in table])))
    for f in filenames:
        if not path.exists(f):
            sys.exit('Could not open file: %s' % f)
        else:
            pass #if verbose: print 'file exists: %s' % f
    return filenames


def getTable(tablefilename, verbose=False):
    f = open(tablefilename, "rb")
    reader = csv.DictReader(f)
    table = [r for r in reader]
    f.close()
    return table





		#outvol = volumeFromInstance(invol, '%s%03d.mnc' % (outbase,depth), dtype='float32', volumeType='float')
#outvol.data[:l,:l,:l] = 1


	# read input volume
	# grow/shrink until desired factor
	# write volume



###############################################################################
# MAIN
###############################################################################
if __name__ == "__main__":
	usage = """usage: %prog [-h/--help] [options] -i INPUT -o OUTPUT"""
	description = """median filter with different kernel sizes. """
	parser = OptionParser(usage=usage, description=description)

	parser.add_option("-i", "--in", dest="input", help="input", type='string', default="")
    parser.add_option("-t", "--table", dest="table", help="table file name, csv file with filename associations", type='string', default="")
	parser.add_option("-o", "--out", dest="outbase", help="outbase", type='string', default="")
	parser.add_option("-f", "--factor", dest="factor", help="growth/shrinkage factor (1 is no change)", type='float', default=1.0)
	parser.add_option("-v", "--volume", dest="volume", help="total cube/brain volume [mm^3]", type='float', default=1)
	parser.add_option("-n", "--subdivisions", dest="subdivisions", help="subdivide cube/brain into N pieces", type='int', default=1)
	parser.add_option("--verbose", dest="verbose", help="more verbose output", action="store_true", default=False)

	(options, args) = parser.parse_args()

	if not path.exists(options.input):
		raise IOError('Could not open input: %s' % options.input)

	if not options.outbase:
		parser.print_usage()
		sys.exit(1)
	
	simulateData(labelfilename,outcsv,replications,seed,verbose=options.verbose)

