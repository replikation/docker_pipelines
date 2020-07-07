#!/usr/bin/env bash

INPUT=$1
SCORE=$2    # 150
LENGTH=$3   # 250
OUTPUT=$4   

< $INPUT awk -v a="$LENGTH" '{if(NR < 2 || $4 >= a) {print}}' | awk -v b="$SCORE" '{if(NR < 2 || $6 >= b) {print}}' > $OUTPUT