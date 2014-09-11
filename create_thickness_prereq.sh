#! /bin/bash
# creates necessary files for cortical thickness estimation
#
# 2013-03-25 jan.scholz@phenogenomics.ca @ MICe, Toronto, Canada
#


usage ()
{
	echo "Usage: $(basename $0) OUTBASE ATLAS"
	echo "OUTBASE  e.g. foo -> foo_cortex.mnc, foo_inside.mnc ..."
	echo "ATLAS    labels"
}

select_and_map ()
{
	IN=$1; OUT=$2; CONST=$3; shift 3
	T=/tmp/$$.mnc
	minc_label_ops2.py $IN $T --select=$@ --binarize
	mincmath -quiet -clobber -constant $CONST -mult $T $OUT
}


###############################################################################
# MAIN
###############################################################################

[ $# -lt 2 ] && { usage; exit 1; }

OUTBASE=$1
ATLAS=$2

TDIR=/tmp
TBASE=$TDIR/`basename $OUTBASE`


RCORTEX="64,181,130,209,133"
LCORTEX="190,180,164,230,131"
CORTEX="$RCORTEX $LCORTEX"
INSIDE="8,16,211,122,6,16,13,112,207,57,59,7,51,55,52,115,23,146,194,66,11,22,106,66,63,12,107,77,159,17,151,155,152,215,103,68"
CC="8,68"
#OLFBULB="5 105"

###############################################################################
# create masks: inside, cortex, outside, resistence
select_and_map $ATLAS ${TBASE}_inside.mnc          1 $INSIDE  
select_and_map $ATLAS ${OUTBASE}_cortex_left.mnc   4 $LCORTEX 
select_and_map $ATLAS ${OUTBASE}_cortex_right.mnc  6 $RCORTEX 
mincmath -quiet -clobber -add ${OUTBASE}_cortex_left.mnc ${OUTBASE}_cortex_right.mnc ${TBASE}_cortex.mnc
select_and_map $ATLAS ${TBASE}_resist.mnc         19 $CC      
miceroi 0 -1 0 -1 $((114-9)) $((8*2)) ${TBASE}_resist.mnc ${TBASE}_resist.mnc

# careful, slice scaling in byte volume results in non-integer label values
mincmath -quiet -clobber -add ${TBASE}_cortex.mnc ${TBASE}_inside.mnc ${TBASE}_resist.mnc ${TBASE}_thickness_labels_tmp.mnc
minc_label_ops2.py --remap=0:10,1:0 ${TBASE}_thickness_labels_tmp.mnc ${OUTBASE}_thickness_labels.mnc

rm ${TBASE}_{cortex,inside,resist}.mnc ${TBASE}_thickness_labels_tmp.mnc

###############################################################################
# lapalacian

minc_label_ops2.py -v --remap=4:5,6:10 ${OUTBASE}_thickness_labels.mnc ${OUTBASE}_thickness_labels_left.mnc
laplacian_thickness -potential_only -use_third_boundary  -from_grid ${OUTBASE}_thickness_labels_left.mnc ${TBASE}_potential_left.mnc
mincreshape -quiet -clobber -dimorder zspace,yspace,xspace ${TBASE}_potential_left.mnc ${OUTBASE}_potential_left.mnc

minc_label_ops2.py -v --remap=6:5,4:10 ${OUTBASE}_thickness_labels.mnc ${OUTBASE}_thickness_labels_right.mnc
laplacian_thickness -potential_only -use_third_boundary  -from_grid ${OUTBASE}_thickness_labels_right.mnc ${TBASE}_potential_right.mnc
mincreshape -quiet -clobber -dimorder zspace,yspace,xspace ${TBASE}_potential_right.mnc ${OUTBASE}_potential_right.mnc

rm ${TBASE}_thickness_labels_{left,right}.mnc

###############################################################################
# marching cubes

#minccalc -quiet -clobber -expression "A[1]>0 ? A[0] : 0"  ${OUTBASE}_potential.mnc  ${OUTBASE}_cortex_left.mnc  ${OUTBASE}_potential_left.mnc
#minccalc -quiet -clobber -expression "A[1]>0 ? A[0] : 0"  ${OUTBASE}_potential.mnc  ${OUTBASE}_cortex_right.mnc  ${OUTBASE}_potential_right.mnc

marching_cubes  ${OUTBASE}_potential_left.mnc  ${OUTBASE}_cortex_left.obj  3.5 7
marching_cubes  ${OUTBASE}_potential_right.mnc ${OUTBASE}_cortex_right.obj 3.5 7

bicobj2iv ${OUTBASE}_cortex_left.obj  > ${OUTBASE}_cortex_left.iv
bicobj2iv ${OUTBASE}_cortex_right.obj > ${OUTBASE}_cortex_right.iv

echo "now simplify/smooth in amira"


