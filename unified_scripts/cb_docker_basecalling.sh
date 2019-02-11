#!/usr/bin/env bash

  # colours, needed for echos
    RED='\033[0;31m'
    YEL='\033[0;33m'
    NC='\033[0m'
    GRE='\033[0;32m'
    DIM='\e[2m'
    ## OPTIONS ##
    WORKDIRPATH=$(pwd) # for docker mountpoint (-v)
    SCRIPTNAME=$(basename -- "$0")
      # variables
      nano_reads=''
      CPU=$(lscpu -p | egrep -v '^#' | wc -l)
      lable='results'
      config='dna_r9.4.1_450bps_flipflop.cfg'
      kittype=''
      flowcell=''
      cpu_mode=''
      batch_mode=''
      sum_mode=''

###############
##  Modules  ##
###############
usage()
  {
    echo "Usage:    $SCRIPTNAME -n fast5folder/ [OPTIONS] [HARDWARE] [APPROACH]"
    echo -e "          ${DIM}e.g. gpu flipflop:${NC}${GRE} $SCRIPTNAME -n fast5folder/ ${NC}"
    echo -e "          ${DIM}e.g. cpu no flipflop:${NC}${GRE} $SCRIPTNAME -n fast5folder/ -c -K SQK-LSK109 -F FLO-MIN106${NC}"
    echo -e "${YEL}Input:${NC}"
    echo -e "          [-n]    Nanopore fast5 folder; ${GRE}-n fast5folder/${NC}"
    echo -e "${YEL}Options:${NC}"
    echo -e "          [-t]    Default: ${GRE}-t ${CPU}${NC} - amount of cores"
    echo -e "          [-l]    Label added to output dir e.g. ${GRE}-l Seq_run_2018${NC}"
    echo ""
    echo -e "${YEL}Hardware:${NC} (default) GPU based local"
    echo -e "          [-c]    Change to CPU based via docker"
    echo ""
    echo -e "${YEL}Approach:${NC} (default) flipflop (dna_r9.4.1_450bps_flipflop.cfg)"
    echo -e "          For normal basecalling use -K and -F"
    echo -e "          [-K]    Kit used e.g. ${GRE}-K SQK-LSK109${NC}"
    echo -e "          [-F]    Flowcell used ${GRE}-F FLO-MIN106${NC}"
    echo -e "${DIM}Advanced"
    echo -e "          [-b]    gpu based batch run (overrides all options) Input: needs [-n]"
    echo -e "          [-s]    collect batch results (overrides all options)${NC}"
    true>&2; exit 1;
  }


guppy_cpu()
  {
    type docker >/dev/null 2>&1 || { echo -e >&2 "${RED}Docker not found. Aborting.${NC}"; exit 1; }
    output="fastq_$label"
    CPU_half=$(echo $(($CPU / 2)))
    flow_option=''
    kit_option=''
    config_option=''
    if [ ! -z "${flowcell}" ]; then flow_option="--flowcell ${flowcell}"; fi
    if [ ! -z "${kittype}" ]; then kit_option="--kit ${kittype}"; else config_option="-c $config"; fi
    mkdir -p ${output}
    docker run --rm -it \
      -v $nano_path:/input \
      -v $WORKDIRPATH/${output}:/output \
      replikation/guppy \
      guppy_basecaller -r -t ${CPU_half} --runners 2 -i /input/ -s /output \
      $flow_option $kit_option $config_option --enable_trimming on --trim_strategy dna -q 0
      exit 1
  }

guppy_gpu()
  {
    type guppy_basecaller >/dev/null 2>&1 || { echo -e >&2 "${RED}guppy locally not found, Aborting.${NC}"; exit 1; }
    output="fastq_$label"
    CPU_half=$(echo $(($CPU / 2)))
    flow_option=''
    kit_option=''
    config_option=''
    if [ ! -z "${flowcell}" ]; then flow_option="--flowcell ${flowcell}"; else config_option="-c $config"; fi
    if [ ! -z "${kittype}" ]; then kit_option="--kit ${kittype}"; else config_option="-c $config"; fi
    mkdir -p ${output}
    guppy_basecaller -r -i $nano_reads -s $output \
    $flow_option $kit_option $config_option --device auto --enable_trimming on --trim_strategy dna -q 0
    exit 1
  }

guppy_batch_gpu()
  {
   echo "Welcome to the batch basecaller for multi-fast5-files"
   echo " a.) Moves one fast5 file to fast5_batch/ and basecalls it."
   echo " b.) After each fast5 file you have 2 min to abort the script."
   echo "     As long as you abort during the 2 min you can always continue another time."
   echo " c.) Retype your command to continue (ignore the batch folders)"
   echo " d.) After all fast5 basecalled use [-s] option to pool everything"
   echo -e "${YEL}Found the following .fast5 files:$NC"
   ls $nano_path
   echo " "
   read -p "Proceed? [yes] or [no] " yn
   case $yn in
       [Yy]* ) echo "Starting Batch..."
                mkdir -p fastq_batch
                for fast5file in $nano_path/*.fast5; do
                  filename=$(basename $fast5file)
                  mkdir -p fast5_batch/${filename%.fast5}
                  mkdir -p fastq_batch/${filename%.fast5}
                  mv $fast5file fast5_batch/${filename%.fast5}
                  guppy_basecaller -i fast5_batch/${filename%.fast5} -s fastq_batch/${filename%.fast5} \
                    -c $config --device auto --enable_trimming on --trim_strategy dna -q 0
                  echo -e "${RED}Finished one fast5 file, waiting 2 min before starting the next one${NC}"
                  echo -e "${RED}Press ctrl + c to exit${NC}"
                  sleep 2m
                done
                ;;
       [Nn]* ) echo "Exiting..."; exit ;;
       * ) echo "  Please answer [y] or [n].";;
   esac
   exit 1
  }

collect()
{
  echo "Welcome to the COLLECTOR, use this after the batch option"
  echo "Your current Working has to contain fast5_batch/ and fastq_batch/"
  if [ -d fast5_batch/ ]; then echo -e "${GRE}fast5_batch found${NC}"; else echo -e "${RED}Can't find fast5_batch${NC}"; exit ; fi
  if [ -d fastq_batch/ ]; then echo -e "${GRE}fastq_batch found${NC}"; else echo -e "${RED}Can't find fastq_batch${NC}"; exit ; fi
  read -p "Proceed? [yes] or [no] " yn
  case $yn in
      [Yy]* ) echo "Starting Batch..."
              mkdir -p fast5/
              mv fast5_batch/*/*.fast5 fast5/
              mkdir -p fastq/
              cat fastq_batch/*/*.fastq > fastq/all_reads.fastq
              head -1 fastq_batch/*/sequencing_summary.txt -q | head -1 > fastq/sequencing_summary.txt
              tail -n+2 fastq_batch/*/sequencing_summary.txt >> fastq/sequencing_summary.txt
              mv fastq_batch/*/*.log fastq/
              # remove !!!empty!!! fast5_batch folder
              # remove fastq_batch
              ;;
      [Nn]* ) echo "Exiting..."; exit ;;
      * ) echo "  Please answer [y] or [n].";;
    esac
    exit 1
}

############################
###   Start of script    ###
############################
echo "                                               _____________________________"
echo "______________________________________________/ Created by Christian Brandt \___"
echo " "

while getopts "n:t:l:K:F:cbsh" arg; do
      case "${arg}" in

        n) nano_reads="${OPTARG}" ;;
        t) CPU="${OPTARG}" ;;
        l) label="${OPTARG}" ;;
        K) kittype="${OPTARG}" ;;
        F) flowcell="${OPTARG}";;
        c) cpu_mode='true';;
        b) batch_mode='true';;
        s) sum_mode='true';;
        h) usage;;
        *) usage
           exit 1 ;;
      esac
done

# getting absolute paths
  nano_path=$(cd "$nano_reads" 2>/dev/null && pwd)

if [ ! -z "${batch_mode}" ]; then
  if [ ! -z "${nano_reads}" ]; then guppy_batch_gpu ; else exit ; fi
fi

if [ ! -z "${sum_mode}" ]; then collect ; else exit ; fi

if [ ! -z "${cpu_mode}" ]; then
  if [ ! -z "${nano_reads}" ]; then guppy_cpu ; fi
fi

if [ -z "${cpu_mode}" ]; then
  if [ ! -z "${nano_reads}" ]; then guppy_gpu ; fi
fi

usage
