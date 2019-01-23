#!/usr/bin/env bash

## Docker ##
  type docker >/dev/null 2>&1 || { echo -e >&2 "${RED}Docker not found. Please run the installation skript, Aborting.${NC}"; exit 1; }

  # CPU cores
    CPU=$(lscpu -p | egrep -v '^#' | wc -l) # can be changed to e.g. CPU="16"
  # colours, needed for echos
    RED='\033[0;31m'
    YEL='\033[0;33m'
    NC='\033[0m'
    ## OPTIONS ##
      WORKDIRPATH=$(pwd) # for docker mountpoint (-v)
      WORKDIRNAME=${PWD##*/} # for docker mountpoint (-v)
      ${FASTQ}=FASTQ
      SCRIPTNAME=$(basename -- "$0")
###############
##  Modules  ##
###############
usage()
  {
    echo "Usage:    $SCRIPTNAME [-c fast5folder/ ] [-i fast5folder/]"
    echo -e "          [-c]    ${YEL}cpu based guppy${NC}"
    echo -e "          [-g]    ${GRE}gpu based guppy${NC}"
    #1>&2; exit 1;
    true>&2; exit 1;
  }

basecaller_ask()
  {
    echo -e "${RED}Input for Basecaller ${NC}"
    echo -e "Enter flowcell type (e.g. ${GRE}FLO-MIN106${NC} or ${GRE}FLO-MIN107${NC}) and hit [Enter]"
    read flowcell
    echo -e "Enter_library kit (e.g. ${GRE}SQK-LSK108${NC}, ${GRE}SQK-RBK004${NC} or ${GRE}SQK-RAD004${NC}) and hit [Enter]"
    read kittype
  }

guppy_cpu()
  {
    mkdir ${FASTQ}
    CPU_half=$(echo $(($CPU / 2)))
    docker run --rm -it -v ${WORKDIRPATH}:/${WORKDIRNAME} replikation/guppy \
    guppy_basecaller -r -t ${CPU_half} --runners 2 -i /${WORKDIRNAME}/${input} -s /${WORKDIRNAME}/${FASTQ}/ \
    --flowcell ${flowcell} --kit ${kittype} --enable_trimming on --trim_strategy dna -q 0
  }

guppy_gpu()
  {
    type guppy_basecaller >/dev/null 2>&1 || { echo -e >&2 "${RED}guppy not found, Aborting.${NC}"; exit 1; }
    guppy_basecaller -i $input -r --device auto -c dna_r9.4.1_450bps_flipflop.cfg -s ${FASTQ} \
    --enable_trimming on --trim_strategy dna
  }

############################
###   Start of script    ###
############################
echo "                                               _____________________________"
echo "______________________________________________/ Created by Christian Brandt \___"
echo " "

while getopts ":c:g:" arg; do
      case "${arg}" in
          c)
              input=${OPTARG}
              basecaller_ask; guppy_cpu
              ;;
          g)
              input=${OPTARG}
              guppy_gpu
              ;;
          *)
              usage
              ;;
      esac
done
shift $((OPTIND-1))
if [ -z "${input}" ]; then usage; fi
