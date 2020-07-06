#!/usr/bin/env bash

INPUT=$1
SCORE=$2    # 150
LENGTH=$3   # 250
OUTPUT=$4   

< $INPUT awk "{if(NR < 2 || $4 >= $LENGTH) {print}}" | awk "{if(NR < 2 || $6 >= $SCORE) {print}}" > $OUTPUT