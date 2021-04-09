#! /usr/bin/bash

# Alexander J. Trouern-Trend (ATrouern-Trend@scripps.edu)
# April 2021

# This script takes two arguments, both fasta files. Each sequence
# of the first compared to each sequence of the second. 

# If sequence from file 1 exists as a subsequnce of any sequence
# in file 2, it is classified as a match and included in the
# output match.fasta. All sequences without matches are returned
#  in nomatch.fasta.

# get seqs only
grep -v '>' $1 | sed 's/\r//g' > shortseqs.txt
grep -v '>' $2 | sed 's/\r//g' > longseqs.txt

# make arrays to categorize short seqs as matching or otherwise
declare -a NO_MATCH
declare -a MATCH

# shortseqs into variable
LINES=$(cat shortseqs.txt)

# Classify each shortseq as either matching to longseq or not.
for line in $LINES; do
    if [ $(grep -c $line "longseqs.txt") -ne 0 ]; then
        MATCH+=($line)
    else
        NO_MATCH+=($line)
    fi
done

# Report number of seqs in each category
echo ${#NO_MATCH[@]} sequences without a match.
echo ${#MATCH[@]} sequences with a match.

# Remove previous attempts and reinitialize outfiles
rm nomatch.fasta
rm match.fasta
touch nomatch.fasta
touch match.fasta

# Create new fastas for each sequence classification.
for SEQ in ${NO_MATCH[@]}; do
    grep -w -B 1 "$SEQ" $1 >> nomatch.fasta
done

for SEQ in ${MATCH[@]}; do
    grep -w -B 1 "$SEQ" $1 >> match.fasta
done
