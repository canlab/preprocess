#!/bin/bash

# STUDY-SPECIFIC VARIABLES:
TRIM=6
FUNC="BASELINE PARTNER RELIEF EMPATHY RESTING"
STRUCT="mprage_5e_RMS"

# STUDY-SPECIFIC FUNCTIONS:
get_expected_run_length() {
    echo $1 | sed 's|.*_\([0-9]*\)_v.*|\1|'
}

get_series_name() {
    echo ${1##*_}
}

get_run_name() {
    echo ${1%%_[0-9][0-9]*}
}


usage() {
    echo "USAGE: ./convert_and_distribute.sh <subdir> [ -s ]"
    echo "Options:"
    echo "   -s    skip dicom conversion"
    echo "<subdir> must contain raw/"
    echo "convert with dcm2nii"
    echo "populate canlab-style directory structure with raw 4D .nii image files (re-named)"
    
    exit 1
}

[[ "$*" =~ "help" ]] && usage

CONVERT=set
SUBDIR=$(readlink -f $1); shift

[ $1 ]&&[ $1 == "-s" ] && unset CONVERT

if [ ! $SUBDIR ] || [ ! -d $SUBDIR/raw ]; then
    usage
fi

cd $SUBDIR


### CONVERT
if [ $CONVERT ]; then
    echo "... DELETING ANY PREEXISTING NIFTIs"
    for f in $FUNC $STRUCT; do
	rm $SUBDIR/raw/$f*/*.nii 2>/dev/null
    done
    
    echo "... CONVERTING DICOMS (with dcm2nii)"
    for f in $STRUCT $FUNC; do
	for i in $(ls -d $SUBDIR/raw/$f*); do
	    cd $i
	    j+=($(qup dcm2nii -b /data/projects/wagerlab/.dcm2nii/dcm2nii.ini `ls *.dcm | head -1`))
	done
    done
    cd $SUBDIR
    snooze -s -j ${j[*]}
    unset j
fi

### DISTRIBUTE
echo "... DISTRIBUTING IMAGES"
# distribute structural
mkdir -p Structural/SPGR
mv -i raw/$STRUCT*/s[0-9]*.nii Structural/SPGR/mprage_rms.nii

# distribute functionals
mkdir -p Functional/{Raw,Preprocessed}
for f in $FUNC; do
    for i in raw/$f*/s[0-9]*.nii; do
	d=$(basename $(dirname $i))

	# check run length
	actualvol=$(fslval $i dim4)
	expectedvol=$(get_expected_run_length $d)
	if [ $actualvol -ne $expectedvol ]; then
	    echo "... WARNING: Ignoring short run $d"
	    continue
	fi

	eval $SETSERIES
	series=$(get_series_name $d)
	runname=$(get_run_name $d)
	
	mkdir -p Functional/Raw/s$series
	mv -i $i Functional/Raw/s$series/$runname.nii
	
	cd ..
	echo "... DISDAQS: trimming $TRIM volumes from $SUBDIR/Functional/Raw/s$series/$runname.nii"
	jorbs+=($(qup ./disdaqs.sh $SUBDIR/Functional/Raw/s$series/$runname.nii $TRIM))
	cd $SUBDIR
    done
done

echo "... WAITING: for disdaqs to complete"
snooze -j ${jorbs}

exit 0
