#!/bin/bash

FILES=( $(ls */*.gz | tr '\n' ' ') )
TRS=( $(ls */TotalReads/*.dsrc | tr '\n' ' ') ) 
DIRF=( $(ls * | tr '\n' ' ') )

#for file in "${FILES[@]}"; do
#    echo "$file"
#    NEW=$( echo $file | sed 's/reseq_//g' )
#    echo $NEW
#    mv "$file" "$NEW"
#done

#declare -a SUBDIR
#
#for file in "${FILES[@]}"; do
#    SUBDIR+=( $(echo $file | cut -d_ -f1) )
#    #echo "${SUBDIR[*]}"
#done
#
#sorted_unique=($(echo "${SUBDIR[@]}" | tr ' ' '\n' | sort -u | grep '441' | tr '\n' ' '))
#echo "${sorted_unique[*]}"
#
#for sub in "${sorted_unique[@]}"; do
#    mkdir Sample_$sub
#    mv ${sub}* Sample_$sub/
#done

#for file in "${FILES[@]}"; do
#    echo "$file"
#    NEW=$( echo $file | sed 's/R1/R1_001/g' | sed 's/R2/R2_001/g' )
#    echo $NEW
#    mv "$file" "$NEW"
#done

# 441-12_S19_R2_001_Total.dsrc to 441-12_2_Total.dsrc

#for file in "${TRS[@]}"; do
#    NEW=$( echo $file | awk -F'[_/]' '{print $4"_"$6"_Total.dsrc"}' | sed 's/R1/1/g' | sed 's/R2/2/g')
#    echo "${NEW}"
#    mv "$file" "$NEW"
#done

for file in "${DIRF[@]}"; do
    NEW=$( echo $file | sed 's/R1/L001_R1/g' | sed 's/R2/L001_R2/g' )
    mv "$file" "$NEW"
done
