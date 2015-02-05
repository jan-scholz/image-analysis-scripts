#!/usr/bin/env python
#
# author: jan.scholz@mouseimaging.ca @ MICe, SickKids, Toronto, Canada
#
# find resonable parameters for deformation simulation
# reads in real data and outputs parameters: means, variance, ...
#   - for each level separately
#   - for multilevel model
#
# ssh -Y topolina
# . /home/matthijs/./scripts/topolina_environment
# cd /projects/mice/jscholz/tmp/def/simulations


from optparse import OptionParser,OptionGroup
from pyminc.volumes.factory import *
#from scipy import ndimage,stats
import numpy as np
from os import path
import sys
import csv
import math


def getTable(tablefilename, verbose=False):
    f = open(tablefilename, "rb")
    reader = csv.DictReader(f)
    table = [r for r in reader]
    f.close()
    return table



def writeSimfilesCSV(table, filename):
	with open(filename, 'wb') as f:
		writer = csv.DictWriter(f,fieldnames=['structure','right label','parent'])
		#writer.writeheader()  # only works from python 2.7 on
		writer.writerows(table)

d
