#! /bin/bash
# create slice images for rotarod publication
# 2013-10-30,  jan.scholz@mouseimaging.ca

LABELSCC=raytrace_colorcoding.txt 
POSCC=(hotred_new  cc_green cc_yellow cc_brown)
NEGCC=(hotblue_new cc_brown cc_brown  cc_blue)
LOWERTHRESH=0
UPPERTHRESH=12

usage ()
{
	echo "Usage: cat coords.txt | $(basename $0) -b BGFILE -l LABELS -s SURFACE -o OUTBASE  STAT1 THRESH1.."
	echo "  -b BGFILE           background "
	echo "  -l LABELS           labels  "
	echo "  -s SURFACE          surfaces, i.e. outline "
	echo "  -o OUTBASE          basename of output, coordinates get appended "
	echo "  STAT1 THRESH1       pairs of statistical overlay maps with associated (lower) threshold"
}
ycoord_to_png ()
{
	ADD="-size 800 800 -bg black -crop -front -sup 3"
	COORDS=`python -c "print ' '.join(['%.3f' % (${1}+i*0.056) for i in range(-1,2)])"`

	for COORD in $COORDS; do
		[ -z "$VERBOSE" ] || echo "processing $COORD"
		SUFF="_y$COORD"
		OBASE=$TDIR/slice$SUFF
		make_slice $BGFILE ${OBASE}.obj y w $COORD

		plane_polygon_intersect $LABELSURFS $TDIR/lines.obj y $COORD
		set_object_colour $TDIR/lines.obj  $TDIR/lines.obj lightgray
		#set_object_opacity $TDIR/lines.obj  $TDIR/lines.obj 0.5

		COMMAND="-output ${OBASE}.rgb -nolight
					-line_width 0.03 $TDIR/lines.obj
					-gray 0 1600 $BGFILE 2 1
					-under transparent -usercc $LABELSCC 1 301 $LABELS -1 0.5"


#echo "ARGS ${ARGS[@]}"

		for ((i=0;i<${#ARGS[@]}/2;i++)); do
			STATS=${ARGS[i*2]}
			LOWERTHRESH=${ARGS[i*2+1]}
#echo "$i"
#echo "STATS $STATS"
#echo "LOWERTHRESH $LOWERTHRESH"
			COMMAND="$COMMAND
				-under transparent -usercc ${POSCC[i]}     $LOWERTHRESH   $UPPERTHRESH $STATS  -1 0.8
				-under transparent -usercc ${NEGCC[i]}    -$LOWERTHRESH  -$UPPERTHRESH $STATS  -1 0.8"
		done

		ray_trace $COMMAND ${OBASE}.obj $ADD
		convert  ${OBASE}.rgb -flip ${OUTBASE}${SUFF}.png
	done
}


###############################################################################
# MAIN
###############################################################################

while getopts b:l:s:T:o:v opt
do
	case "$opt" in
		b)  BGFILE="$OPTARG";;
		l)  LABELS="$OPTARG";;
		s)  LABELSURFS="$OPTARG";;
		#i)  STATS="$OPTARG";;
		#t)  LOWERTHRESH="$OPTARG";;
		T)  UPPERTHRESH="$OPTARG";;
		o)  OUTBASE="$OPTARG";;
		v)  VERBOSE=1;;
		\?)  usage; exit 1;;
	esac
done
shift $(expr $OPTIND - 1)

[ -z "$TMPDIR" ] && { echo "ERROR: \$TMPDIR not set"; exit 1; }
[ "$( tty )" = 'not a tty' ] || { usage; exit 1; }

[ -f "$BGFILE" ]     || { echo "ERROR: could not find BACKGROUND_FILE: $BGFILE"; exit 1; }
[ -f "$LABELS" ]     || { echo "ERROR: could not find LABELS: $LABELS"; exit 1; }
[ -f "$LABELSURFS" ] || { echo "ERROR: could not find LABELSURFS: $LABELSURFS"; exit 1; }
#[ -f "$STATS" ]      || { echo "ERROR: could not find STATS: $STATS"; exit 1; }

ARGS=($@)

TDIR=$TMPDIR/$$;
mkdir -p $TDIR
trap "{ rm -fr $TDIR; echo 'cleaned up temp. directory'; exit 255; }" SIGINT SIGTERM


while read YCOORD; do
	ycoord_to_png $YCOORD
done

rm -rf $TDIR
exit 0




# cd /projects/mice/jscholz/rot/stats_hr_new2+flip/dti_to_basket

# Website with all the information:
# 
# http://www.bic.mni.mcgill.ca/~david/Ray_trace/ray_trace_tutorial.html

# Content of /tmp/colour_coding
# 
# 1.000 1.0000 0.0000 0.0000 
# 2 0 1 0
# 180  0 0 0
# 190 0 0 1
# 255 1 1 1



# create cluster-extend thresholded signed images
# minccalc -expression 'A[0]>56?A[1]:0' clusters_dti_meanCmeanFA_lmerp_mean-p-sign_to_basket_thr0.005.mnc dti_meanCmeanFA_lmerp_mean-p-sign_to_basket.mnc clusters_dti_meanCmeanFA_lmerp_mean-p-sign_to_basket_thr0.005_gt56.mnc
# minccalc -expression 'A[0]>56?A[1]:0' clusters_dti_faCmeanFA_lmerp_grouprotarod-p-sign_to_basket_thr0.005.mnc  dti_faCmeanFA_lmerp_grouprotarod-p-sign_to_basket.mnc  clusters_dti_faCmeanFA_lmerp_grouprotarod-p-sign_to_basket_thr0.005_gt56.mnc
# minccalc -expression 'A[0]>56?A[1]:0' clusters_final_groupPsexRsideRid_lmerp_grouprotarod-p-corr-sign_thr0.1.mnc  final_groupPsexRsideRid_lmerp_grouprotarod-p-corr-sign.mnc clusters_final_groupPsexRsideRid_lmerp_grouprotarod-p-corr-sign_thr0.1_gt56.mnc
# minccalc -expression 'A[0]>56?A[1]:0' clusters_final_meanPsexRsideRid_lmerp_mean-p-corr-sign_thr0.1.mnc  final_meanPsexRsideRid_lmerp_mean-p-corr-sign.mnc clusters_final_meanPsexRsideRid_lmerp_mean-p-corr-sign_thr0.1_gt56.mnc 


#awk 'NR>1&&NR<5 {print $4}' clusters_dti_faCmeanFA_lmerp_grouprotarod-p-sign_to_basket_thr0.005.mm.txt | ./raytrace.sh -v -b template_basket_masked.mnc -l resampled_atlas_new_stead+ull_right.mnc -s blender_constrained_smooth.obj -i clusters_dti_faCmeanFA_lmerp_grouprotarod-p-sign_to_basket_thr0.005_gt56.mnc -t 0.995 -o slices/dti_group

#awk 'NR>1&&NR<44 {print $4}' clusters_dti_meanCmeanFA_lmerp_mean-p-sign_to_basket_thr0.005.mm.txt | ./raytrace.sh -v -b template_basket_masked.mnc -l resampled_atlas_new_stead+ull_right.mnc -s blender_constrained_smooth.obj -i clusters_dti_meanCmeanFA_lmerp_mean-p-sign_to_basket_thr0.005_gt56.mnc -t 0.995 -o slices/dti_behav_mean

#awk 'NR>1&&NR<11 {print $4}' clusters_final_meanPsexRsideRid_lmerp_mean-p-corr-sign_thr0.1.mm.txt | ./raytrace.sh -v -b template_basket_masked.mnc -l resampled_atlas_new_stead+ull_right.mnc -s blender_constrained_smooth.obj -o slices/vol_beahv_mean clusters_final_meanPsexRsideRid_lmerp_mean-p-corr-sign_thr0.1_gt56.mnc 0.9
######awk 'NR>1&&NR<11 {print $4}' clusters_final_meanPsexRsideRid_lmerp_mean-p-corr-sign_thr0.1.mm.txt | ./raytrace.sh -v -b template_basket_masked.mnc -l resampled_atlas_new_stead+ull_right.mnc -s blender_constrained_smooth.obj -i clusters_final_meanPsexRsideRid_lmerp_mean-p-corr-sign_thr0.1_gt56.mnc -t 0.9 -o slices/vol_beahv_mean

#awk 'NR>1&&NR<21 {print $4}' clusters_final_groupPsexRsideRid_lmerp_grouprotarod-p-corr-sign_thr0.1.mm.txt | ./raytrace.sh -v -b template_basket_masked.mnc -l resampled_atlas_new_stead+ull_right.mnc -s blender_constrained_smooth.obj -o slices/vol_group clusters_final_groupPsexRsideRid_lmerp_grouprotarod-p-corr-sign_thr0.1_gt56.mnc 0.9
######awk 'NR>1&&NR<21 {print $4}' clusters_final_groupPsexRsideRid_lmerp_grouprotarod-p-corr-sign_thr0.1.mm.txt | ./raytrace.sh -v -b template_basket_masked.mnc -l resampled_atlas_new_stead+ull_right.mnc -s blender_constrained_smooth.obj -i clusters_final_groupPsexRsideRid_lmerp_grouprotarod-p-corr-sign_thr0.1_gt56.mnc -t 0.9 -o slices/vol_group

# paste this into svg file:
# i=0; for f in `ls slices-src/dti_behav_mean_y*png | sort -r`; do x=$(($i%3*210-500)); y=$(($i/3*150)); echo "<image xlink:href=\"$f\" x=\"$x\" y=\"$y\" height=\"170.8\" width=\"223.3\" />"; let i+=1;  done`

# 3.40
# z = 3.14


# hippocampus hotspot of changes
# printf " -0.63\n-0.69\n" | ./raytrace.sh -v -b template_basket_masked.mnc -l resampled_atlas_new_stead+ull_right.mnc -s blender_constrained_smooth.obj -o slices/hotspot clusters_final_groupPsexRsideRid_lmerp_grouprotarod-p-corr-sign_thr0.1_gt56.mnc 0.9 clusters_dti_faCmeanFA_lmerp_grouprotarod-p-sign_to_basket_thr0.005_gt56.mnc 0.995 clusters_dti_meanCmeanFA_lmerp_mean-p-sign_to_basket_thr0.005_gt56.mnc 0.995


