#! /bin/bash
# converts a series of tiff slices to a minc volume
#
# 2012-05-28 jan.scholz@phenogenomics.ca at MICe
# based on Leila's code

# ADDARGS="-swap_bytes -2"
ADDARGS="-mri -short -signed"
# ORIENTATION=-transverse
ORIENTATION=-xyz
#ORIENTATION="-dimorder x,y,z"

usage ()
{
	echo "Usage: $(basename $0) OUTFILE XSTEP YSTEP ZSTEP TIFFFILE.."
	echo "  OUTFILE   output file"
	echo "  XSTEP     the step size in x-direction [mm]"
	echo "  YSTEP     the step size in y-direction [mm]"
	echo "  ZSTEP     the step size in z-direction [mm]"
	echo "  TIFFFILE  one or more TIFF files"
	echo ""
	echo "./convertTIFF2MINC.sh volume.mnc 0.1 0.1 0.1 *.tif"
}

[ $# -lt 4 ] && { usage; exit 1; }
[ -z "$TMPDIR" ] && { echo "ERROR: \$TMPDIR not set, run: export TMPDIR=/tmp"; exit 1; }


OUTFILE=$1
XSTEP=" $2"
YSTEP="-$3"
ZSTEP="-$4"
shift 4


DIMS=`identify -format '%h %w' $1`
TMP=($@)
N=${#TMP[@]}

TMPFILE=$TMPDIR/tmp_tiff2minc_$$.mnc

for i in $@; do
	convert $i GRAY:-;
done | rawtominc $ADDARGS $ORIENTATION -zstep $ZSTEP -ystep $YSTEP -xstep $XSTEP -xstart 0 -ystart 0 -zstart 0 -clobber $TMPFILE $N $DIMS

if [ -f "$TMPFILE" ]; then
	cp $TMPFILE $OUTFILE
else
	echo "ERROR: conversion failed"
	exit 1
fi

rm -f $TMPFILE
exit 0


