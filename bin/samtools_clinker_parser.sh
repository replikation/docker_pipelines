#!/bin/bash        

FILE=$1
SEARCH=$2
RANGE=$3


grep "$SEARCH" $FILE | cut -f2-4 | awk -v myvar="$RANGE" 'BEGIN { FS=OFS="\t" } {print $1,$2-myvar,$3+myvar}' | sed 's/'$'\t''/:/' | sed 's/'$'\t''/-/' > seq_subsections.txt