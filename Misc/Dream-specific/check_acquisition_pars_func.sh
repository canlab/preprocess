#!/bin/bash

REPORTNAME=acqpars.txt
RANDNUM=$RANDOM
DCMDIR=/data/archive/human/dicom/triotim/twager
[ $(hostname) == "dream" ] && DREAM=set

PARAMETERS="\
AcquisitionOrder \
SliceThickness \
RepetitionTime \
EchoTime \
NumberOfAverages \
EchoNumbers \
AcquisitionMatrix \
FlipAngle \
ScanningSequence \
SequenceVariant \
ScanOptions \
MRAcquisitionType \
SequenceName \
AngioFlag \
PixelSpacing"

usage() {
    echo "Usage: ./check_acquisition_pars_func.sh [options]"
    echo "Description:"
    echo "   This script will check acquisition parameters for functional scans in a study."
    echo "   It will check parameters ${PARAMETERS// /, }, ImagingFrequency, SpacingBetweenSlices, and Acquisition Order"
    echo "      in the first dicom in each functional series in the study."
    echo "   A report (<study>_acqpars.txt) is printed containing the number of instances of each unique value for each parameter."
    echo "   If any parameter has more than one unique value in the study, a separate report (<study>_acqpars_<parameter>.txt) is generated"
    echo "      with a line for each series and the value for the 'offending' parameter."
    echo "Input options:"
    [ $DREAM ] && echo "   SINGLE ARGUMENT: <study>"
    [ $DREAM ] && echo "      suited for looking at all functional images in a twager study on dream"
    [ $DREAM ] && echo "      will report on all functional series dicom directories (*_v0*) in $DCMDIR/<study>*"
    echo "   MULTIPLE ARGUMENTS: <studyname> <functional series dcm dir> [...]"
    echo "      first string is study name (for use in output file naming)"
    echo "      subsequent strings are dicom directories for functional series to report on"
    echo "      allows user to specify inputs instead of looking in $DCMDIR"
    echo
    echo "Notes:"
    echo "   Acquisition order is pulled from the readable strings in the dicom and not the header."
    echo
    echo "Examples:"
    [ $DREAM ] && echo "   this will report on all functional series in $DCMDIR/ilcp_200008:"
    [ $DREAM ] && echo "   $ ./check_acquisition_pars_func.sh ilcp"
    [ $DREAM ] && echo
    echo "   this will report on all functional series for task runs in ilcp subjects 30-39:"
    echo "   $ ./check_acquisition_pars_func.sh ilcp ilcp/Imaging/ilcp3*/raw/*TASK*"

    exit 1
}

mktempfail() {
    echo "ERROR: failed to create temporary file with mktemp; exiting" >&2
    exit 1
}


### PARSE INPUTS
[[ "$*" =~ help ]] && usage

if [ $1 ]&&[ $1 == "parameters" ]; then
    printf "%s\n" $PARAMETERS
    exit 1
fi

if [ $1 ]&&[ $2 ]; then
    STUDY=$1
    shift
    while [ $1 ]; do
	if [ -d $1 ] && [ $(ls $1/*dcm 2>/dev/null | wc -l) -ge 1 ]; then
	    SERIES+=($1)
	else
	    echo "ERROR: $1 is not a valid dicom directory."
	fi
	shift
    done
elif [ $1 ]; then
    STUDY=($(readlink -f $DCMDIR/$1*))
fi
    
if [ ! $STUDY ]||[ ! -d $STUDY ]||[ ${#STUDY[*]} -gt 1 ]; then
    if [ $DREAM ]; then
	echo "POSSIBLE STUDIES:"
	ls -d1 $DCMDIR/* | sed -e 's|.*/||' -e 's|_[^_]*$||' -e 's|^|   |'
	printf "\n%s\n" "TO SEE USAGE STATEMENT: ./check_acquisition_pars_func.sh --help"
	exit 0
    else
	usage
    fi
else
    SERIES=($(ls -d $STUDY/*/*/*_v0*_[0-9][0-9][0-9][0-9]))
fi

# prep
if [ $DREAM ]; then
    td=$(mktemp -d -p /data/tmp) || mktempfail
else
    td=acqparswork$$
    mkdir $td
fi
#echo "Temporary Directory: $td"

# get dicom headers
echo "... STUDY: $STUDY"
for series in ${SERIES[*]}; do
    firstdcm=$(ls $series/*dcm 2>/dev/null | head -1)
    if [ "$firstdcm" ]; then
	subj=$(dirname $series | sed -e "s|${STUDY}/||" -e 's|/|_|g')
	dcmdump $firstdcm > $td/${subj}_$(basename $series)
	strings $firstdcm | grep sSliceArray.ucMode | sed -e 's|sSliceArray\.ucMode| AcquisitionOrder|' -e 's|0x1|Ascending_Sequential|' -e 's|0x2|Descending_Sequential|' -e 's|0x4|Ascending_Interleaved|' >> $td/${subj}_$(basename $series)
    else
	echo "ERROR: no dicoms in $series"
    fi
done

# check for consistency in parameters
report=$(basename ${STUDY%_*})_${REPORTNAME}
rm $report 2> /dev/null
echo "... PREPARING REPORT: $report"

# add acquisition order (needs extra parsing)
#echo "AcquisitionOrder" >> $report
#if [ $(awk '/sSliceArray\.ucMode/ {if ($3=="0x1") p="Ascending_Sequential"; if ($3=="0x2") p="Descending_Sequential"; if ($3=="0x4") p="Ascending_Interleaved"; print p}' $td/* | | sort -n | uniq -c | tee -a $report | wc -l) -gt 1 ]; then
#    (cd $td; awk '/sSliceArray\.ucMode/ {if ($3=="0x1") p="Ascending_Sequential"; if ($3=="0x2") p="Descending_Sequential"; if ($3=="0x4") p="Ascending_Interleaved"; print p,FILENAME}' *) | sort -n -k1 > ${report%.txt}_AcquisitionOrder.txt
#    echo "WARNING! detected more than one parameter for AcquisitionOrder (check ${report%.txt}_AcquisitionOrder.txt)"
#fi
    
# loop through standardly-formatted parameters
for par in $PARAMETERS; do
    echo $par >> $report
    if [ $(awk '/ '${par}'/ {print $3}' $td/* | sed -e 's|[][]||g' -e 's|\\| |g' | sort -n | uniq -c | tee -a $report | wc -l) -gt 1 ]; then
	(cd $td; awk '/ '$par'/ {print $3,FILENAME}' *) | sort -n -k1 > ${report%.txt}_$par.txt
	echo "WARNING! detected more than one parameter for $par (check ${report%.txt}_$par.txt for details)"
    fi
    echo >> $report
done

# add fields that require formatted print
echo ImagingFrequency >> $report
if [ $(sed 's|[][]||g' $td/* | awk '/ImagingFrequency/ {printf("%.2f\n",$3)}' | sort -n | uniq -c | tee -a $report | wc -l) -gt 1 ]; then
    (cd $td; awk '/ImagingFrequency/ {printf ("%.2f %s\n",$3,FILENAME)}' *) | sort -n -k1 > ${report%.txt}_ImagingFrequency.txt
    echo "WARNING! detected more than one parameter for ImagingFrequency (check ${report%.txt}_ImagingFrequency.txt for details)"
fi
echo >> $report

echo SpacingBetweenSlices >> $report
if [ $(sed 's|[][]||g' $td/* | awk '/SpacingBetweenSlices/ {printf("%.2f\n",$3)}' | sort -n | uniq -c | tee -a $report | wc -l) -gt 1 ]; then
    (cd $td; awk '/SpacingBetweenSlices/ {printf ("%.2f %s\n",$3,FILENAME)}' *) | sort -n -k1 > ${report%.txt}_SpacingBetweenSlices.txt
    echo "WARNING! detected more than one parameter for SpacingBetweenSlices (check ${report%.txt}_SpacingBetweenSlices.txt for details)"
fi


# clean up
echo "... CLEANING UP"
rm -r $td

exit 0
