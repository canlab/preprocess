#!/bin/bash

grp=wagerlab
archivedir=/dreamio/archive/human/dicom/triotim/twager
studiesdir=/data/projects/wagerlab/labdata/current
corrections_fname=accession_id_corrections.txt
ionode=$(readlink -f . | cut -d/ -f2)

[ "$ionode" == "dream" ] && echo "ERROR: ionode isn't getting set right (it's ${ionode}). Ask ruzic@colorado.edu for help!"
# ideally this spawns rsync jobs that run on the io node appropriate for your data
# but if it's finding the wrong node (getting the head node instead) then 
#    manually set ionode (for wagerlab this should be dreamio3)
# or remove the ssh statement from the rsync command below to run on the current working node

usage() {
    echo "Usage: ./dicomsdicomsdicoms.sh <STUDY> [options]"
    echo "Tabulates a list of all subjects/scanning sessions stored in /data/archive under <STUDY>"
    echo "  prints with acquisition date, number of series, and accession number."
    echo "If the file \"${corrections_fname}\" exists in <STUDY>/Imaging, it will be"
    echo "  used to substitute correct accession IDs for bad or missing ones"
    echo "  (whitespace-delimited text file containing 2 columns: URSI correct_ID)"
    echo "ASSUMPTION: dicom directories are $archivedir/*/M*/Study*/*"
    echo "Options"
    echo "  -i            interactively offer to copy directories into your study directory (uses rsync on io node)"
    echo "  -u            update: try to copy all scanning sessions/subjects that don't appear to already exist"
    echo "  -a            automatic: proceed without asking for permission"
    echo "  -b <basename> prepend <basename> to IDs when making directory names in interactive mode"
    echo "  --noheader    don't print a header for the table"
    echo
    echo "POSSIBLE STUDIES:"
    ls -d1 $archivedir/* 2> /dev/null | sed -e 's|.*/||' -e 's|_[^_]*$||'
    echo

    exit 1
}


### SET UP
[[ "$1" =~ help ]] && usage

[ $1 ] && STUDY=$1 || usage
shift
HEADER=set
while [ $1 ]; do
    if [ "$1" == "-i" ]; then
	INTERACTIVE=set
	unset UPDATE
    elif [ "$1" == "-u" ]; then
	UPDATE=set
	unset INTERACTIVE
    elif [ "$1" == "-a" ]; then
	AUTOMATIC=set
    elif [ "$1" == "-b" ]; then
	shift
	BASE=$1
    elif [ "$1" == "--noheader" ]; then
	unset HEADER
    else
	echo "UNRECOGNIZED OPTION: $1"
	exit 1
    fi
    shift
done

dicomdirs=($(ls -d $archivedir/${STUDY}* 2> /dev/null))
studydir=($(ls -d $studiesdir/* | grep -i $STUDY))

### ERROR CHECKING
if [ ${#studydir[*]} -gt 1 ]; then
    echo "$STUDY does not pick out a single study directory:"
    echo "$studydir"
    exit 1
fi

if [ ${#dicomdirs[*]} -gt 1 ]; then
    echo "ERROR: $STUDY does not pick out a unique study:"
    echo "${dicomdirs[*]}"
    exit 1
fi
if [ ! $dicomdirs ]; then
    echo "POSSIBLE STUDIES:"
    ls -d1 $archivedir/* 2> /dev/null | sed -e 's|.*/||' -e 's|_[^_]*$||'
    exit 1
fi

corrections_file=$studydir/Imaging/$corrections_fname
if [ -f $corrections_file ]; then
    echo "... USING corrections file: $corrections_file"
fi

subs=($(ls -d1 ${dicomdirs}/M*/Study* 2> /dev/null))

if [ ${#subs[*]} -eq 0 ]; then
    echo "no subjects found in $dicomdirs"
    exit 1
fi


### GET INFO
n=0
for s in ${subs[*]}; do
    l[n]=${#s}

    onedicom=$(ls $s/*RMS*/*.dcm 2>/dev/null | head -1)
    if [ ! $onedicom ]; then
	for i in $(ls -d $s/*); do onedicom=$(ls $i/*dcm 2> /dev/null | head -1); [ $onedicom ] && break; done
    fi
    if [ ! $onedicom ]; then
	echo "ERROR: can't find dicom in $s"
    fi

    
    accession[n]="$(dcmdump $onedicom | grep -i accession | grep -v "no value available" | sed 's|.*\[\(.*\)\].*|\1|')"
    if [ -f $corrections_file ]; then
	u=$(basename $(dirname $s))
	while read ursi corrected; do
	    [ $ursi == $u ] && accession[n]=$(echo $corrected | sed 's|.*\[\(.*\)\].*|\1|')
	done < $corrections_file
    fi
    [ ! "${accession[n]}" ] && accession[n]="NOT_FOUND"

    ndir[n]=$(ls -d $s/*_[0-9][0-9][0-9][0-9] | wc -l)

    date[n]="$(dcmdump $onedicom | grep -i seriesdate | sed 's|.*\[\([0-9]*\)\].*|\1|')"
    
    let n++
done

# determine longest filename
l=$(echo ${l[*]} | tr ' ' '\n' | sort -n | tail -1)



if [ ! $UPDATE ]; then
    ### PRINT
    if [ $HEADER ]; then    
	printf "%${l}s   %10s  %3s   %s\n" directory date n ID
	echo "-----------------------------------------------------------------------------------------------------------"
    fi
    
    for ((n=0;n<${#subs[*]};n++)); do
	printf "%${l}s   %10s  %3s   %s\n" ${subs[n]} "${date[n]}" ${ndir[n]} "${accession[n]}"
    done | sort -n -k2
fi


### COPY if desired
if [ $INTERACTIVE ]; then
    d=0
    echo
    td=$(mktemp -d -p /data/tmp)
    while :; do
	echo
	read -p "ID (leave blank to $([ $d -eq 0 ] && echo quit || echo begin copying)): " s
	
	[ "$s" ] || break

	if [ $s == "NOT_FOUND" ]; then
	    echo "NOT_FOUND unallowable"
	    continue
	fi

	targdir=$studydir/Imaging/${BASE}${s}/raw
	unset found
	for ((i=0;i<=${#subs[*]};i++)); do
	    if [ "${accession[i]}" == "$s" ]; then
		found=set
		if [ ! -d $targdir ]; then
		    echo "SOURCE: ${subs[i]}" 
		    echo "TARGET: ${targdir}"
		    [ $AUTOMATIC ] && ans=Y || read -p "PROCEED? (y/[n])  " ans
		    if [[ "$ans" =~ ^[Yy]$ ]]; then
			id[d]=$s
			src[d]=$(readlink -f ${subs[i]})
			trg[d]=$targdir
			let d++
		    fi
		else
		    echo "ERROR: $targdir ALREADY EXISTS"
		fi
	    fi
	done

	[ ! $found ] && echo "ID not found: $s"
    done
elif [ $UPDATE ]; then
    d=0
    echo
    td=$(mktemp -d -p /data/tmp)
    for ((i=0; i<${#subs[*]}; i++)); do
	if [ "${accession[i]}" == "NOT_FOUND" ]; then
	    echo
	    echo "SKIPPING \"NOT_FOUND\""
	    continue
	fi

	targdir=$studydir/Imaging/${BASE}"${accession[i]}"/raw

	if [ -d ${targdir/Imaging/Imaging\/Unused} ]; then
	    echo "${targdir/Imaging/Imaging\/Unused} ALREADY EXISTS"
	    continue
	fi

	if [ ! -d $targdir ]; then
	    echo
	    echo "SOURCE: ${subs[i]}"
	    echo "TARGET: ${targdir}"
	    [ $AUTOMATIC ] && ans=Y || read -p "PROCEED? (y/[n])  " ans
	    if [[ "$ans" =~ ^[Yy]$ ]]; then
		id[d]="${accession[i]}"
		src[d]=$(readlink -f ${subs[i]})
		trg[d]=$targdir
		let d++
	    fi
	else
	    echo "$targdir ALREADY EXISTS"
	fi
    done
fi

# run rsync commands
if [ "$d" ]&&[ $d -gt 0 ]; then
    echo
    echo -n "COPYING DIRECTORY:"
    for ((i=0;i<d;i++)); do
	echo -n " "${id[i]}
        mkdir -p ${trg[i]}
	ssh $ionode "rsync -a ${src[i]}/ ${trg[i]}; chgrp -R $grp ${trg[i]}"
    done
    echo
fi

exit 0
