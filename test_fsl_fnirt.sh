#! /bin/bash
# testing for FNRIT (build 418) bug
# 2011-08-26

FNIRTOPS="--logout=/dev/null --applyrefmask=1 --applyinmask=1 --subsamp=2 --miter=2 --infwhm=0.3 --reffwhm=0.3 --lambda=100 --estint=1 --warpres=10,10,10 --regmod=bending_energy --intmod=global_non_linear_with_bias --intorder=5 --biasres=100,100,100 --biaslambda=500"


mse ()
{
	# calculate mean square error
	${FSLDIR}/bin/imrm tmp
	${FSLDIR}/bin/fslmaths $2 -sub $1 -sqr -mas ref_mask tmp
	${FSLDIR}/bin/fslstats tmp -m
	${FSLDIR}/bin/imrm tmp
}


echo "testing for FNIRT (build 418) bug"
echo
echo -n "date "; date
echo "FSLDIR $FSLDIR"
echo -n "FNIRT build "; fnirt 2>&1 | grep build
flirt -version
echo


mkdir -p fnirtbug
cd fnirtbug/
${FSLDIR}/bin/imcp $FSLDIR/data/standard/MNI152_T1_1mm_brain ref
${FSLDIR}/bin/imcp $FSLDIR/data/standard/MNI152_T1_1mm_brain_mask ref_mask
printf "0.984808  0  -0.173648  0.828463\n0  1  0  -1\n0.173648  0  0.984808  -0.51938\n0  0  0  1\n" > aff_orig.mat
${FSLDIR}/bin/flirt -in ref -ref ref -applyxfm -init aff_orig.mat -out in || exit 1
echo "finished setting up, running tests"

${FSLDIR}/bin/flirt -in in -ref ref -omat aff.mat -out aff -nosearch
MSEAFF=`mse ref aff`
printf "affine    to ref MSE %10.0f\n" $MSEAFF

${FSLDIR}/bin/fnirt --in=in --ref=ref --aff=aff.mat --iout=nlin1 --cout=nlin1_warp $FNIRTOPS
MSENLIN1=`mse ref nlin1`
printf "nlin1     to ref MSE %10.0f\n" $MSENLIN1

${FSLDIR}/bin/fnirt --in=in --ref=ref --inwarp=nlin1_warp --iout=nlin2 --cout=nlin2_warp $FNIRTOPS
MSENLIN2=`mse ref nlin2`
printf "nlin2     to ref MSE %10.0f\n" $MSENLIN2

${FSLDIR}/bin/applywarp -i in -r ref --premat=aff.mat -w nlin2_warp -o nlin2b
MSENLIN2b=`mse ref nlin2b`
printf "aff+nlin2 to ref MSE %10.0f\n" $MSENLIN2b

[ 1 -eq `echo "$MSENLIN2 > $MSENLIN2b * 100" | bc` ] && { echo; echo "*** FOUND FNIRT BUG ***"; exit 17; }

exit 0

