#!/usr/bin/env python
#
# author: jan.scholz@mouseimaging.ca @ MICe, SickKids, Toronto, Canada
#
# TODO:
#   - extend dictionaries to 2d (filename, kernelrad)
#   - speed up correlation function
#
# sys.path.append(os.path.abspath('/micehome/jscholz/bin/'))
# import corrclust
# k = corrclust.cubicKernel([2])
# filenames = ["rot29_hr_lsq6_flip-log-determinant-fwhm0.1-200um.mnc"]
# f = filenames[0]
# invol = volumeFromFile(f)
# for kr,k in enumerate(kernels):
# ndimage.filters.median_filter(invol.data, footprint=k, mode='constant', cval=0.0)

from optparse import OptionParser,OptionGroup
from pyminc.volumes.factory import *
from scipy import ndimage,stats
import numpy as np
from os import path
import sys
import csv

#import resource
# print memory consumption
#def Using(point):
#    usage=resource.getrusage(resource.RUSAGE_SELF)
#    usage=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
#    return '''%s: usertime=%s systime=%s mem=%s mb'''%(point,usage[0],usage[1],(usage[2]*resource.getpagesize())/1000000.0)


def appendSuffix(filename, suffix=''):
    (root, ext) = path.splitext(filename)
    return ''.join((root,suffix,ext))


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
#with open("input", "r") as inp, open("output", "w") as out: out.write(inp.read())


def cubicKernel(radia):
    d = dict.fromkeys(radia)
    for r in radia:
        d[r] = np.ones((r,r,r), dtype=np.int)
    return d


def readImages(table, filenamecolumn='filename', verbose=False):
    # reads image volumes and puts them into dictionary with filename key
    filenames = checkFiles(table, filenamecolumn, verbose)
    vols =  dict.fromkeys(filenames)
    for f in filenames:
        if verbose: print "processing %s" % f
        vols[f] = volumeFromFile(f)
    return vols


def kernelFilter(vols, kernel, kernelradius, suffix='_kernel', verbose=False):
    # creates filtered versions, writes them to disk (because volFromInst cannot be prevented from doing that)
    outvols = dict.fromkeys(vols.keys())
    for f in vols.keys():
        kr,k = kernelradius,kernel
        outvols[f] = volumeFromInstance(vols[f], appendSuffix(f,suffix='%s%02i' % (suffix,kr)), dtype='float32', volumeType='float')
        ndimage.filters.median_filter(vols[f].data, footprint=k, output=outvols[f].data, mode='constant', cval=0.0, origin=0)
        outvols[f].writeFile()
        outvols[f].closeVolume()
    # maybe filter should return 4d volume already? but that's not split left/right
    return outvols


# works only with one kernel (more kernels need to be address in 2d)
def optimalKernel(vols, table, kernels, suffix='_kernel', verbose=False):

    corrfiletable = [] # can be removed. only to ggplot r vs. kernel radius

    for kr,k in kernels.iteritems():
        if verbose: print 'filtering with radius %0.2f' % kr
        filtervols = kernelFilter(vols, k, kr, suffix=suffix, verbose=verbose)

###############################################################################
# move to read fucntion
        # combine to 4d volume
        volshape = filtervols.values()[0].data.shape
        n = len(filtervols)
        volsall   = np.zeros(volshape + (n,))   #, dtype=np.int)

        # intreleave left and right in 4d volume
        ids = sorted(list(set([e['id'] for e in table])))
        for i,mouseid in enumerate(ids):
            rightfile = [e['filename'] for e in table if e['id']==mouseid and e['side']=='right'][0]
            leftfile  = [e['filename'] for e in table if e['id']==mouseid and e['side']=='left'][0]
            volsall[:,:,:,i*2  ] = filtervols[rightfile].data 
            volsall[:,:,:,i*2+1] = filtervols[leftfile].data 
            if verbose: print 'read mouseid %s, kernel radius %0.2f' % (mouseid,kr) 

###############################################################################
# move into correlation function
        def corfunc(a):
            np.seterr(invalid='ignore')  # turn of NAN produces error
            #return np.corrcoef(a[::2],a[1::2])[1,1] # doesn't work?
            return stats.pearsonr(a[::2],a[1::2])[0] # this works and returns p-value, too!

        rr = np.apply_along_axis(corfunc, 3, volsall)
        
        outputname = 'corr%s%02i.mnc' % (suffix,kr)
        corrfiletable.append({'filename' : outputname, 'radius' : kr})
        v = volumeFromInstance(vols.values()[0], outputname, dtype='float32', volumeType='float')
        v.data = np.squeeze(rr)
        # DATA SHOULD ALSO BE WRITTEN INTO 4D ARRAY TO DETERMIN MAx
        # corrvols[:,:,:,counter] = rr
        #print 'v.data.dtype:', v.data.dtype, 'v.data at 20,20,20:', v.data[20,20,20]
        v.writeFile()
        v.closeVolume()

    with open('table_corrfiles.csv', 'wb') as f:
        writer = csv.DictWriter(f,corrfiletable[0].keys())
        #writer.writeheader()  # only works from python 2.7 on
        writer.writerows(corrfiletable)

    # maximum kernel number for each voxel from corrvol4d
    return 1



###############################################################################
# MAIN
###############################################################################
if __name__ == "__main__":
    usage = """usage: %prog [-h/--help] [options] MINCFILE.."""
    description = """median filter with different kernel sizes. """
    parser = OptionParser(usage=usage, description=description)
    
    parser.add_option("-m", "--mask", dest="mask",
                      help="mask name, find maximum within masked area",
                      type='string', default="")
    parser.add_option("-t", "--table", dest="table",
                      help="table file name, csv file with filename associations",
                      type='string', default="")
    parser.add_option("--suffix", dest="suffix",
                      help="suffix appended to input files for output file names",
                      type='string', default='_tfce')
    parser.add_option("-v", "--verbose", dest="verbose",
                      help="more verbose output",
                      action="store_true", default=False)

    (options, args) = parser.parse_args()
    
    if options.mask and not path.exists(options.mask):
        raise IOError('Could not open mask: %s' % options.mask)
    
    if not options.table:
        parser.error("Table file not supplied")

    table = getTable(options.table, verbose=options.verbose)

    kernels = cubicKernel([3,5])
    #kernels = cubicKernel([1,2,3,4,5,6])
    vols = readImages(table, filenamecolumn='filename', verbose=options.verbose)
    r = optimalKernel(vols, table, kernels, verbose=options.verbose)


#    for f in args:
#        if not path.exists(f):
#            raise IOError('Could not open file: %s' % f)

#    evs = getTable('foo',args, verbose=options.verbose)
#    kernels = cubicKernel([5])
#    medianvols = readImages(args, kernels, verbose=options.verbose)
#    optimalKernel(medianvols, kernels, evs)




###############################################################################

# for f in corr_kernel??.mnc; do mincmath -quiet -clobber -constant 1 -mult $f `basename $f .mnc`s.mnc; done && mincconcat -clobber -concat_dimension time corr_kernel??s.mnc out.mnc && mincstats out.mnc

#    if verbose: print 'reading table file: %s' % tablefilename
#    d = dict.fromkeys(filenames)
#    for f in filenames[::2]:
#        d[f] = {'side' : 'left'}
#    for f in filenames[1::2]:
#        d[f] = {'side' : 'right'}
#    return d


#    if len(args) < 1:
#        parser.error("Incorrect number of arguments")

# get all ids, unique, sorted
# sorted(list(set([e['id'] for e in out])))
# get filenames for certain id
# [e['filename'] for e in out if e['id']==i]
# restrict to certain side/hemisphere
# [e['filename'] for e in out if e['id']==i and e['side']=='left']

#    # arrange according to ids, to correlate l/r of the same individual
#    for i,k in enumerate(vols.keys()):
#        side = evs[k]['side']
#        if side == 'left':
#            vols4d[:,:,:,i] = vols[k].data
#        elif side == 'right':
#            vols4d[:,:,:,i] = vols[k].data
#        else:
#            raise ValueError('Illeagel side value %s' % side)



#    for kr,k in kernels.iteritems():
#        outvol[f] = volumeFromInstance(invol[f], appendSuffix(f,suffix='%s%02i' % (suffix,kr)), dtype='ushort')
#        ndimage.filters.median_filter(invol[f].data, footprint=k, output=outvol[f].data, mode='constant', cval=0.0, origin=0)
#        outvol[f].writeFile()
#        outvol[f].closeVolume()

# library(RMINC)
# t <- read.csv('table_all_localqqk.csv')
# p <- c(20,20,20)
# x <- mincGetVoxel(t$filename,p)
# cor.test(x[seq(1,8,by=2)],x[seq(2,8,by=2)])
# 
# o<- array(NA,c(6,6)); for (i in 1:6) { v <- mincGetVolume(t$filename[i]); for (j in 1:6) { o[i,j] <- which.min(abs(v-coords[j])); }}

#        volsleft  = np.zeros(volshape + (n/2,)) #, dtype=np.int)
#        volsright = np.zeros(volshape + (n/2,)) #, dtype=np.int)
#            volsright[:,:,:,i] = filtervols[rightfile].data 
#            volsleft[:,:,:,i]  = filtervols[leftfile].data 

#        x = volsleft
#        y = volsright

        #xm = np.expand_dims(np.mean( volsleft,  axis=3), axis=3) # across subjects/time
        #ym = np.expand_dims(np.mean( volsright, axis=3), axis=3)
        #xs = np.expand_dims(np.std(  volsleft,  axis=3), axis=3)
        #ys = np.expand_dims(np.std(  volsright, axis=3), axis=3)

        #volshape4d = volshape + (1,)

        #nom = ((x-xm)*(y-ym)).sum(axis=3)
        #nom.shape = volshape4d
        #denom = (xs*ys)

        #np.seterr(invalid='ignore')  # turn of NAN produces error
        #r = nom / (denom + 0.00001)
        #r = np.nan_to_num(r)
        #r = np.squeeze(r)
        #sum( (x-xm)*(y-ym) ) / (xs*ys)


    #     concatenate((a1, a2, ...)[, axis])    # tuple(d.values())
    # calc correlation                       # left volumes          right volumes
    #   - this is the correct one: np.corrcoef(left_ijk,right_ijk)  for each voxel
    #   - scipy.ndimage.filters.correlate1d(vol4d[:2:],        weights=vol4d[:2:])  # middle value idx or idx+1
    #     could use 4dvols_left and weights 4dvols_right and then middle slice contains all the values
#yunbiased = y-np.mean(y)
#ynorm = np.sum(yunbiased**2)
#acor = np.correlate(yunbiased, yunbiased, "same")/ynorm
## use only second half
#acor = acor[len(acor)/2+1:len(acor)/2+2]   # autocorrelatio with lag one

#numpy.mean(a, axis=)
#numpy.std(a, axis=)

#sum( x-mean(x)*y-mean(y) ) / std(x)*std(y)
