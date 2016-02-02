#!/usr/bin/python
# formally known as distortion_correction_november2007.pl
# modified by jan.scholz@phenogenomics.ca

# Takes a brain MR scan and will apply the november 2007 distortion
# correction transformation to it:
#
# -MR to CT transform
#
# Will determine the coil the mouse was scanned in
# using the vnmr:coil entry in the MINC header

# TODO:
#   - return queue number
#   - job array

version = "0.1";

import sys
import os
import optparse
import subprocess


basedir = "/projects/mice/matthijs/distortion-correction-november-2007/"

coil1_xfm = basedir + "/coil1-phantom-17oct2007/coil1-mri-to-CAD-lsq6-nlin-extra-large-super-inverted-lsq6.xfm"
coil2_xfm = basedir + "/coil2-phantom-18oct2007/coil2-mri-to-CAD-lsq6-nlin-extra-large-super-inverted-lsq6.xfm"
coil3_xfm = basedir + "/coil3-phantom-16oct2007/coil3-mri-to-CAD-lsq6-nlin-extra-large-super-inverted-lsq6.xfm"

coil1_likefile = basedir + "/coil1-phantom-17oct2007/img_17oct07.0-reshaped-coil1-likefile.mnc"
coil2_likefile = basedir + "/coil2-phantom-18oct2007/img_18oct07.1-reshaped-coil2-likefile.mnc"
coil3_likefile = basedir + "/coil3-phantom-16oct2007/img_16oct07.2-reshaped-coil3-likefile.mnc"


def validcoil(x):
	try:
		coil = int(x)
	except ValueError:
		raise

	if not coil in range(1,17):
		raise ValueError

	return coil


def splitarg(a):
	filename = a.split(":")[0]
	coil = float("nan")

	try:
		coil = validcoil(a.split(":")[1])
	except IndexError:
		pass
	except ValueError:
		print >> sys.stderr, 'ERROR: not a valid coil number:', a.split(":")[1]
		sys.exit(1)

	return (filename, coil)


def getcoil(filename):
	try:
		p = subprocess.Popen(['mincinfo','-attvalue','vnmr:coil',filename], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	except OSError:
		pass

	if p.wait():
		if p.communicate()[1].startswith('Error reading file.'):
			coil = float("nan")
		else:
			print >> sys.stderr, 'ERROR: could not read file', filename
			sys.exit(1)
	else:
		try:
			coil = int(p.communicate()[0].strip())
		except ValueError:
			print >> sys.stderr, 'ERROR: read invalid coil number from file:', filename
			sys.exit(1)

	return (filename, coil)



###############################################################################
# MAIN
###############################################################################

def main():
	usage = "%prog [options] --output-dir OUTPUTDIR MINCFILE [MINCFILE..]\n" + 9*" " + "-h for help"
	p = optparse.OptionParser(usage=usage, version="%prog " + version)
	p.add_option('--output-dir', '-o',                   dest="outdir",        default="discorr", help="output directory (default: discorr)")
	p.add_option('-n', '--dry-run', action="store_true", dest="dryrun",        default=False,     help="Print commands without executing them")
	p.add_option('--interp', '-i',                       dest="interp",        default="sinc",    help="resample interpolation (sinc,tricubic,trilinear)")
	p.add_option('--likefile', '-l',                     dest="likefile",      default={1:coil1_likefile,2:coil2_likefile,3:coil3_likefile}, help="likefiles")
	p.add_option('--queues', '-q',                       dest="queues",        default='',        help="queues (default: any)")
	p.add_option('--verbose', '-v', action="store_true", dest="verbose",       default=False,     help="Prints more information.")

	options, arguments = p.parse_args()

	if len(arguments) < 1:
		p.error("no minc file to process")

	options.interp = '-' + options.interp

	for arg in arguments:
		(filename,coil) = splitarg(arg)
		if coil == coil and getcoil(filename)[1] != coil:
			print >> sys.stderr, 'WARNING: coil in file header "%s" does not match requested coil "%s"' % (getcoil(filename)[1],coil)

		print 'sge_batch', options.queues, '-l vf=4G -J dc-4G-rot36_b0 mincresample -2', options.interp, 'transform /projects/mice/matthijs/distortion-correction-november-2007/coil3-phantom-16oct2007/coil3-mri-to-CAD-lsq6-nlin-extra-large-super-inverted-lsq6.xfm -like ../raw/rot03_b0.mnc', filename, 'discorr/rot36_b0.november_2007_distortion_corrected.mnc'

if __name__ == '__main__':
	main()
	sys.exit(0)

