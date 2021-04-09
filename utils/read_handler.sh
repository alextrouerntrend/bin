#!/bin/bash
# =======================================
#
# Usage: bash read_handler.sh -e [experiment number] -d [destination format] -c [True or False] -g [./dir/to/harvest/files] -o [dropbox folder for upload]
#
#        -e | --experiment (e.g. Exp444)
#        -d | --destinationformat (fastq, fastq.gz, dsrc) NOTE: Can only handle 1 per run.
#        -c | --concatenate (if true, will concatenate multiple lanes of same sample (keeps paired ends separate for PE data).
#        -h | --harvest (if none given, new files will not be copied)
#        -o | --dropbox 
#  
# =======================================
# NOTES:
#     Need to build master list of samples
#     Try out SPRING fastq compression, a state-of-the-art lossless compressor
#     add single-end support
#     MUST ENTER ABSOLUTE PATH! NO TILDA SHORTCUT!
#     add multithreading support
#     add gather function
#     add export function
# =======================================

printf "\n\n\n"
DIVI='============================='

# Argument Parsing

POSITIONAL=()
while [[ $# -gt 0 ]]
do
key="$1"

case $key in
    -e|--experiment)
    EXP="$2"
    shift
    shift
    ;;
    -d|--destinationformat)
    DEST="$2"
    shift
    shift
    ;;
    -c|--concatenate)
    CONCAT="$2"
    shift
    shift
    ;;
    -h|--harvest)
    HAR="$2"
    shift
    shift
    ;;
    -o|--dropbox)
    DBOX="$2"
    shift
    shift
    ;;
    --default)
    DEFAULT=YES
    shift
    ;;
    *)
    POSITIONAL+=("$1")
    shift
esac
done
set -- "${POSITIONAL[@]}"


# Confirm harvest is dir

if [[ -d $HAR ]]; then
    cd $HAR
    HAR=$(pwd)
    cd -
        :
elif [[ -z $HAR ]]
    then
        :
else
    echo Make sure -h is set to a directory.
    exit 1
fi

echo $DIVI
echo "EXPERIMENT         = ${EXP}"
echo "DESTINATION FORMAT = ${DEST}"
echo "HARVEST LOCATION   = ${HAR}"
echo "CONCATENATE LANES  = ${CONCAT}"
echo $DIVI

# Confirm experiment path exists 

if test -d /gpfs/group/flits/ProjectData/$EXP
    then
        printf "Found Experiment "$EXP"\n"
        DIR=/gpfs/group/flits/ProjectData/$EXP
elif test -d $EXP
     then
        printf "Found Experiment at "$EXP"\n"
        DIR=$EXP
else
    read -p 'Experiment ID not found in /gpfs/group/flits/ProjectData/. If it exists elsewhere, please indicate the absolute path: ' DIR
    if test -d $DIR
        then  
            echo "Found Experiment "$EXP
    fi
fi


# test argument 2
ARG2=(fastq fastq.gz dsrc)
if [[ "${ARG2[@]}" =~ "$DEST" ]]; then
    :
else
    echo "Destination format not understood! Recognized formats are fastq, fastq.gz and dsrc"
    exit 1
fi

cd $DIR
echo $DIVI

# Recursively search directory structure for each filetype, produce array of filenames and measure length
FQGZ=()
while IFS=  read -r -d $'\0'; do
    FQGZ+=("$REPLY")
done < <(find . -name "*.fastq.gz" -print0)
FQGZ_C=${#FQGZ[@]}

FQ=()
while IFS=  read -r -d $'\0'; do
    FQ+=("$REPLY")
done < <(find . -name "*.fastq" -print0)
FQ_C=${#FQ[@]}

DS=()
while IFS=  read -r -d $'\0'; do
    DS+=("$REPLY")
done < <(find . -name "*.dsrc" -print0)
DS_C=${#DS[@]}

TOT="$(( $DS_C + $FQ_C + $FQGZ_C ))"

# Report numbers
printf "Found $TOT files.\n\n"
echo of these:
echo $FQGZ_C are gzipped fastq
echo $FQ_C are uncompressed fastq
printf "$DS_C are dsrc \n\n"
echo $DIVI

# Identify attributes
attribute_test () {
    printf "Identifying library attributes...\n\n"
    if find . -name "*_R2_*" &> /dev/null; then
        PE=true
        printf "Found paired-end reads\n"; else
        PE=false
        printf "Found single-end reads\n"
    fi

    if [[ -n $(find -name "*_L001_*") ]] && [[ -n $(find -name "*_L002_*") ]]; then
        UNCAT=true
        printf "Found unconcatenated fastq files.\n\n"; else
        UNCAT=false
        printf "Did not find unconcatenated fastq files.\n\n"
        echo $DIVI
    fi
}

attribute_test





# Concatenate (NEEDS TO BE FINISHED)
concatenate () {
    if $UNCAT && $PE; then
        ENDS=(R1 R2) # Each end concatenated separately
        for end in "${ENDS[@]}"
        do
            for name in $(find . -name "*L001_R1_*" | awk 'BEGIN { FS = "_L001_" } {print $1}')
            do
                if [[ $DEST == "fastq.gz" ]]; then
                    echo Concatenating .fastq.gz files by sample $name and read $end.
                    cat $name*$end* > ${name}_${end}.fastq.gz
                fi
                if [[ $DEST == "fastq" ]]; then
                    echo Concatenating .fastq files by sample
                    cat $name*$end* > ${name}_${end}.fastq
                fi
            done
        done
    fi
    if $UNCAT && ! $PE; then
        echo "Single-end support not yet available"; exit
    fi
}

concatenate

# Four base functions for formatting

_gunzip () {
    echo "Unzipping zipped files: "
    echo ${FQGZ[*]}              # list files being decompressed
    gunzip --verbose ${FQGZ[*]}
    FQGZ=("${FQGZ[@]/fastq.gz/fastq}")    # rename elements of array to reflect unzipping   
    for i in ${FQGZ[@]};
    do
        if [[ $FQ[*] =~ $i ]]; then
            :
        else
            FQ+=("$i")
        fi
    done 
    FQGZ=()                      # gzipped array is now empty
    echo $DIVI
}

_gzip () {
    echo "Gzipping fastq files: "
    echo ${FQ[*]}
    gzip --verbose ${FQ[*]}
    FQ=("${FQ[@]/fastq/fastq.gz}")    
    for i in $FQ[@];
    do
        if [[ $FQGZ[*] =~ $i ]]; then
            :
        else
            FQGZ+=("$i")
        fi
    done

    FQGZ=("${FQGQ[@]}" "${FQ[@]}")
    FQ=()
    echo $DIVI
}

_dsrc_c () { 
    echo "DSRC compressing fastq files: "
    echo ${FQ[*]}
    for i in ${FQ[@]}
    do
        NP="${i##*/}"     # Remove path from element
        NE="${NP%.fastq}" # Remove extension from element
        LOC="${i%/*}"
       /gpfs/group/flits/atrouern/tools/dsrc c $i ${LOC}/${NE}.dsrc
    done
    FQ=("${FQ[@]/fastq/dsrc}")       
    for i in ${FQ[@]};
    do
        if [[ $DS[*] =~ $i ]]; then
            :
        else
            DS+=("$i")
        fi
    done
    FQ=()
    echo $DIVI
}

_dsrc_d () {
    echo "DSRC decompressing .drsc files: "
    echo ${DS[*]}
    for i in ${DS[@]}
    do
        NP="${i##*/}"
        NE="${NP%.dsrc}"
        LOC="${i%/*}"
        /gpfs/group/flits/atrouern/tools/dsrc d $i ${LOC}/${NE}.fastq
    done
    DS=("${DS[@]/dsrc/fastq}")
    for i in ${DS[@]};
    do
        if [[ $FQ[*] =~ $i ]]; then
            :
        else
            FQ+=("$i") 
        fi
    done
    DS=()
    echo $DIVI
}

_dest_fq () {
    if [[ ${#DS[@]} -gt 0 ]]; then
        _dsrc_d
    fi
    if [[ ${#FQGZ[@]} -gt 0 ]]; then
        _gunzip
    fi
    if [[ -d $HAR ]]; then
        for i in "${FQ[@]}"
        do
            cp -v "${i}" "${HAR}"
        done
    fi
}

_dest_fqgz () {
    if [[ ${#DS[@]} -gt 0 ]]; then
        _dsrc_d
    fi
    if [[ ${#FQ[@]} -gt 0 ]]; then
        _gzip
    fi
    if [[ -d $HAR ]]; then
        for i in ${FQGZ[@]}
        do
            cp -v "${i}" "${HAR}"
        done
    fi
}

_dest_dsrc () {
    if [[ ${#FQGZ[@]} -gt 0 ]]; then
        _gunzip
    fi
    if [[ ${#FQ[@]} -gt 0 ]]; then
        _dsrc_c
    fi
    if [[ -d $HAR ]]; then
        for i in ${DS[@]}
        do
            cp -v ${i} ${HAR}
        done
    fi
}

_main () {
    if [[ $DEST == "fastq" ]]; then
        _dest_fq
    elif [[ $DEST == "fastq.gz" ]]; then
        _dest_fqgz
    elif [[ $DEST == "dsrc" ]]; then
        _dest_dsrc
    fi
}

_main

# Dropbox Upload
if [[ -v DBOX ]]; then
    echo $DIVI
    echo Uploading to Dropbox at $DBOX
    sh /gpfs/group/flits/atrouern/tools/Dropbox-Uploader/dropbox_uploader.sh upload $HAR $DBOX
    echo Dropbox upload complete!
fi

echo Done!
