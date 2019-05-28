#!/bin/bash
# This is a simple script to remove sequence from a fasta file based on matching with blast.

# Basic Usage Info
if [ "$1" == "-h" ]; then
  echo "Usage: `basename $0` <target-sequence.fa> <candidate-sequence.fa>"
  exit 0
fi

source blast-2.6.0
blastn -query $1 -subject $2 -perc_identity 95 -outfmt 6 -max_target_seqs 1 > results-$1-$2.txt

# define some variables

var1=$(head -1 results-$1-$2.txt | awk '{ print $4 }')

declare -i match=75000

length=$(head -1 results-$1-$2.txt | awk '{ print $4 }')

start=$(head -1 results-$1-$2.txt | awk '{ print $9 }')

end=$(head -1 results-$1-$2.txt | awk '{ print $10 }')

#remain=$(tail -n +2 $2 | tr -d '[:space:]' | wc -c)

if [ "$var1" -gt "$match" ]; 
then
    echo "SGI4 Present";
	echo "start $start"
	echo "end $end"
	echo "length $length"
	echo "remain $remain"
	# make file 1 line
	if [ "$end" -gt "$start" ]; 
	then
	tail -n +2 $2 | tr -d '[:space:]' | cut --complement -c$start-$end | fold -w 80 > $2-remove-$1.txt
	echo ">"$2"-remove-"$1"" | cat - $2-remove-$1.txt > temp && mv temp $2-remove-$1.txt
    else
	tail -n +2 $2 | tr -d '[:space:]' | cut --complement -c$end-$start | fold -w 80 > $2-remove-$1.txt
	echo ">"$2"-remove-"$1"" | cat - $2-remove-$1.txt > temp && mv temp $2-remove-$1.txt
    fi;
else
    echo "SGI4 Absent";
fi;