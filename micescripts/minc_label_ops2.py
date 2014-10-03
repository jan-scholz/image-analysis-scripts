#!/usr/bin/python

from pyminc.volumes.factory import *
from numpy import *
from optparse import OptionParser



if __name__ == "__main__":

    usage = "usage: %prog [options] input_labels.mnc output_labels.mnc"
    description = """write description here"""

    parser = OptionParser(usage=usage, description=description)

    parser.add_option("--merge", dest="merge",
                      help="Second label file to whose labels are to be merged with the input_labels",
                      type="string")
    parser.add_option("--remap", dest="remap",
                      help="Label numbers to replace with new label numbers. For example, to turn label 3 into label 11, it would be --remap 3:11. More than one can be specified at once, separated by commas, i.e. --remap 3:11,21:55.",
                      type="string")
    parser.add_option("--select", dest="select",
                      help="Label numbers to select (comma separated, i.e. 3,7,211) - outfile will be just these labels",
                      type="string")
    parser.add_option("--remove", dest="remove",
                      help="Label numbers to remove (comma separated, i.e. 3,7,211) - outfile will exclude these labels",
                      type="string")
    #parser.add_option("--clean", help="convert to round integers", action="store_true", default=False)
    parser.add_option("--binarize", help="threshold >0.5 and binarize output", action="store_true", default=False)
    parser.add_option("-v", "--verbose", help='verbose output', action='store_true', default=False)

    (options, args) = parser.parse_args()

    if len(args) != 2:
        parser.error("Incorrect number of arguments")


    inlabels = volumeFromFile(args[0], dtype='ubyte')
    outlabels = volumeLikeFile(args[0], args[1])

    if options.verbose:
        print("Labels found in input: %s" % unique(inlabels.data[::]))

    if options.remap:
        print 'v', options.clean
        outlabels.data[::] = inlabels.data[::]
        labelpairs = options.remap.split(',')
        for pair in labelpairs:
            labels = pair.split(':')
            if options.verbose:
                print("Mapping: %d to %d" % (int(labels[0]), int(labels[1])))
            outlabels.data[inlabels.data == int(labels[0])] = int(labels[1])

    if options.select:
        labels = options.select.split(',')
        for l in labels:
            outlabels.data[inlabels.data == int(l)] = int(l)

    if options.remove:
        outlabels.data[::] = inlabels.data[::]
        labels = options.remove.split(',')
        for l in labels:
            outlabels.data[inlabels.data == int(l)] = 0

    if options.merge:
        addlabels = volumeFromFile(options.merge, dtype='ubyte')
        # should add check to make sure dimensions are the same

        # blank one set of labels - i.e. they replace rather than add
        inlabels.data[addlabels.data > 0] = 0
        # now add them together
        outlabels.data[::] = inlabels.data[::] + addlabels.data[::]
        addlabels.closeVolume()

    if options.binarize:
        outlabels.data[outlabels.data > 0.5] = 1

#    if options.clean:
#        outlabels.data[::] = inlabels.data[::]
#        #print('foo')
#        #for l in range(min(inlabels.data[::]),max(inlabels.data[::])+1):
#        #   if options.verbose:
#        #       print "Mapping: %d" % l
#        #   outlabels.data[outlabels.data > 0.5] = 1

    # write to file
    outlabels.writeFile()
    outlabels.closeVolume()
    inlabels.closeVolume()

