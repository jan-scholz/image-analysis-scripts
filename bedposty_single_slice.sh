#!/bin/sh

subjdir=$1
slice=$2
shift
shift
opts=$*

slicezp=`${FSLDIR}/bin/zeropad $slice 4`

COMMAND="${FSLDIR}/bin/xfibres --data=$subjdir/data_slice_$slicezp --mask=$subjdir/nodif_brain_mask_slice_$slicezp -b $subjdir/bvals -r $subjdir/bvecs --forcedir --logdir=$subjdir.bedpostX/diff_slices/data_slice_$slicezp $opts  > $subjdir.bedpostX/logs/log$slicezp" 

echo "command = $COMMAND"
$COMMAND


#COMMAND="echo \"\`date +\"%y/%m/%d %H:%M:%S\"\`\" \$HOSTNAME"; ${FSLDIR}/bin/xfibres --data=$subjdir/data_slice_$slicezp --mask=$subjdir/nodif_brain_mask_slice_$slicezp -b $subjdir/bvals -r $subjdir/bvecs --forcedir --logdir=$subjdir.bedpostX/diff_slices/data_slice_$slicezp $opts  > $subjdir.bedpostX/logs/log$slicezp"
#$COMMAND && echo Done 
