#!/usr/bin/env bash
< centrifuge_results.out awk '{if(NR < 2 || $4 >= 250) {print}}' | awk '{if(NR < 2 || $6 >= 150) {print}}' \
> centrifuge_filtered.out