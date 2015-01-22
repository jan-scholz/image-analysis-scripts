#!/usr/bin/python
#
# /projects/mice/jscholz/rot/reg_hr_new2/hr_all_lsq6s+flip_resampled_atlas
#


import sys

from pyminc.volumes.factory import *
from numpy import *
from optparse import OptionParser

if __name__ == "__main__":
	usage = "usage: %prog [options] "
	parser = OptionParser(usage)
	parser.add_option("-i", "--magetdir", help='directory containing output of MAGeT.py', type="string")
	parser.add_option("-v", "--verbose", help='verbose output', action='store_true', default=False)
	(options, args) = parser.parse_args()

	if len(args):
		parser.print_usage()
		sys.exit(1)

	print 'foo'

	sys.exit(0)

# build up library
#  - input: lsq12 images, labels
#  - normalize images
#  - normalize patches?

# subject selection
#  sum of the squared difference (SSD) across the initialization mask (union of all subjects for specific labels),
#  N closest subjects are retained during the entire segmentation process

# seach volume definition
#  limited search volume Vi, defined as a cube centered on the voxel xi under study. Thus, within each of the N selected subjects, we search for similar patches in a cubic region around the location under study. This search volume can be viewed as the intersubject variability of the structure of interest in stereotaxic space.

# patch preselection
#  structural similarity measure (SSIM) 
#  ss = (2*m_i * m_s,j) / (m_i^2 + m_s,j^2)  *   (2*s_i * s_s,j) / (s_i^2 + s_s,j^2)
#  @ voxel i; subject s, location j 
#  value of ss is greater than a given threshold th, the intensity distance between patches i and j is computed. The threshold th was set to 0.95 for all the experiments.
#  precompute Patch mean and variance?

# 

prep="""
cd /projects/mice/jscholz/rot/reg_hr_new2/hr_all_lsq6s+flip_processed

LIKE=../hr_all_lsq6s+flip_nlin/nlin-3.mnc
ATLAS=../hr_all_lsq6s+flip_resampled_atlas/resampled_atlas_new.mnc
DIRS="rot0?_hr_lsq6"

# align to native space (here lsq6?)
for d in $DIRS; do
IN=$ATLAS
OUT=lr_12/${d}_`basename $IN .mnc`-resampled-final-nlin
A="mincresample -2 -clobber -nearest_neighbour -like $LIKE -invert_transformation -transform $d/transforms/${d}-final-nlin.xfm $IN ${OUT}.mnc";
IN=${d}/resampled/${d}-resampled-lsq12.mnc
OUT=lr_12/${d}_`basename $IN .mnc`-resampled-final-nlin
C="mincresample -2 -clobber -tricubic -like $LIKE -invert_transformation -transform $d/transforms/${d}-final-nlin.xfm $IN ${OUT}.mnc";
echo "$A; $C"
done | parallel -v

# downsample to speed up testing
for d in $DIRS; do
IN=$ATLAS
OUT=lr_12/${d}_`basename $IN .mnc`-resampled-final-nlin
B="`autocrop -noexecute -clobber -isostep 0.1 ${OUT}.mnc ${OUT}_step0-100.mnc| cut -d' ' -f5-` -nearest_neighbour"
IN=${d}/resampled/${d}-resampled-lsq12.mnc
OUT=lr_12/${d}_`basename $IN .mnc`-resampled-final-nlin
D=`autocrop -noexecute -clobber -isostep 0.1 ${OUT}.mnc ${OUT}_step0-100.mnc| cut -d' ' -f5-`
echo "$B; $D"
done | parallel -v


#D="mincresample -clobber -like lr_12/${d}_`basename $IN .mnc`-resampled-final-nlin_step0-100.mnc ${OUT}.mnc ${OUT}_step0-100.mnc"

"""
