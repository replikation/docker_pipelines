#!/usr/bin/env bash

## Docker ##
  type docker >/dev/null 2>&1 || { echo -e >&2 "${RED}Docker not found. Please run the installation skript, Aborting.${NC}"; exit 1; }
  # use clonal pipeline by default
    meta=''
    fwd_reads=''
    rev_reads=''
    nano_reads=''
    L_wtdbg2='10000'
  # CPU cores
    CPU=$(lscpu -p | egrep -v '^#' | wc -l) # can be changed to e.g. CPU="16"
  # colours, needed for echos
    RED='\033[0;31m'
    YEL='\033[0;33m'
    NC='\033[0m'
    GRE='\033[0;32m'
    ## Parameters ##
      WORKDIRPATH=$(pwd) # for docker mountpoint (-v)
      WORKDIRNAME=${PWD##*/} # for docker mountpoint (-v)
      SCRIPTNAME=$(basename -- "$0")

###############
##  Modules  ##
###############

usage()
  {
    echo "Usage:    $SCRIPTNAME [-1 illumina_fwd.fastq ] [-2 illumina_rev.fastq] [-n nanopore.fastq] [-m] [-t <integer>]"
    echo "Inputs:"
    echo -e "          [-1]    ${YEL}Illumina fastq forward reads${NC}"
    echo -e "          [-2]    ${YEL}Illumina fastq reverse reads${NC}"
    echo -e "          [-n]    ${YEL}Nanopore fastq file${NC}"
    echo "Options:"
    echo -e "          [-m]    Default: no metagenome - add ${GRE}-m${NC} to change to metagenome assembly"
    echo -e "          [-t]    Default: ${GRE}-t ${CPU}${NC} - amount of cores"
    echo -e "          [-L]    Default: ${GRE}-L ${L_wtdbg2}${NC} - ONT only option, see wtdbg2"
    exit;
  }

wtdbg2_clonal()
{
echo "Starting wtdbg2 assembly"
# change to correct paths
# untested
  output="wtdgb2_assembly"
  mkdir -p $output
    docker run --rm -it -v ${WORKDIRPATH}:/${WORKDIRNAME} \
                        -v $nano_path:/input_nano \
                        -v ${WORKDIRPATH}/${output}:/output \
    replikation/wtdbg2 \
      wtdbg2 -t $CPU -i /input_nano/${nano_file} -o /output/${nano_file%.fastq} -L $L_wtdbg2
    #create contigs
    docker run --rm -it -v ${WORKDIRPATH}/${output}:/output \
    replikation/wtdbg2 \
      wtpoa-cns -t $CPU -i /output/${nano_file%.fastq}.ctg.lay.gz \
      -fo /output/${nano_file%.fastq}.fa
# remove ctg.lay.gz
# polish with medaka or unicycler_polish
exit 0
}

wtdbg2_meta()
{
  # change to correct paths
  # untested
  echo "nano: $nano_reads"
  echo "cpu: $CPU"
  echo "execute wtdbg2 assembler with special edges and stuff"
  echo "medaka polish would be usefull maybe"
  docker run --rm -it -v ${WORKDIRPATH}:/${WORKDIRNAME} replikation/wtdbg2 \
  wtdbg2 -t $CPU -x ont -e ${wtdbg2_edge_depth} -S ${wtdbg2_subsampling} --rescue-low-cov-edges -i /${WORKDIRNAME}/${fastqfile} -o /${WORKDIRNAME}/${FASTA_raw}/${filename%.fastq}_${wtdbg2_L1}  -L ${wtdbg2_L1}
  # L2
}

unicycler_illumina_only()
{
  echo "fwd: $fwd_reads"
  echo "rev: $rev_reads"
  echo "cpu: $CPU"
  echo "execute unicycler with just -1 and -2"
}

unicycler_hybrid()
{
  echo "fwd: $fwd_reads"
  echo "rev: $rev_reads"
  echo "nano: $nano_reads"
  echo "cpu: $CPU"
  echo "execute unicycler with all the flags for hybrid assembly -1 and -2"
}

meta_illumina_only()
{
  echo "fwd: $fwd_reads"
  echo "rev: $rev_reads"
  echo "check out whats the best illumina only assembly mybe meta_spades?"
  echo "cpu: $CPU"
}

meta_hybrid()
{
  # untested
  output="opera-ms_assembly"
  # unzip illumina if .gz - if fastq nothing happens
  gunzip $fwd_reads 2>/dev/null
  gunzip $rev_reads 2>/dev/null
  mkdir -p $output
  # create confige file
  mkdir -p config
  echo "OUTPUT_DIR /output" > config/config.file
  echo "ILLUMINA_READ_1 /input_fwd/${fwd_file%.gz}" >>  config/config.file
  echo "ILLUMINA_READ_2 /input_rev/${rev_file%.gz}" >>  config/config.file
  echo "LONG_READ /input_nano/$nano_file" >>  config/config.file
  echo "NUM_PROCESSOR $CPU" >>  config/config.file
  #OPTIONAL PARAMETERS
  echo "STRAIN_CLUSTERING YES" >>  config/config.file
  echo "CONTIG_LEN_THR 500" >>  config/config.file
  echo "CONTIG_EDGE_LEN 80" >>  config/config.file
  echo "CONTIG_WINDOW_LEN 340" >>  config/config.file
  echo "KMER_SIZE 60" >>  config/config.file
  echo "LONG_READ_MAPPER blasr" >>  config/config.file
  #echo "CONTIGS_FILE sample_files/sample_contigs.fasta"

  # run assembly
  docker run --rm -it \
    -v ${WORKDIRPATH}:/${WORKDIRNAME} \
    -v $fwd_path:/input_fwd \
    -v $rev_path:/input_rev \
    -v $nano_path:/input_nano \
    -v ${WORKDIRPATH}/${output}:/output \
    -v ${WORKDIRPATH}/config/:/config \
    replikation/opera_ms /config/config.file
  rm -fr config/config.file
}

#############################
###   Start of script    ####
#############################
echo "                                               _____________________________"
echo "______________________________________________/ Created by Christian Brandt \___"
echo " "

# you could add a output flag
while getopts '1:2:n:mt:L:' flag; do
    case "${flag}" in
      1) fwd_reads="${OPTARG}" ;;
      2) rev_reads="${OPTARG}" ;;
      n) nano_reads="${OPTARG}" ;;
      m) meta='true' ;;
      t) CPU="${OPTARG}" ;;
      L) L_wtdbg2="${OPTARG}" ;;
      *) usage
         exit 1 ;;
    esac
done

# getting dir names
  fwd_dir=$(dirname "$fwd_reads") 2>/dev/null
  rev_dir=$(dirname "$rev_reads") 2>/dev/null
  nano_dir=$(dirname "$nano_reads") 2>/dev/null
# getting absolute paths
  fwd_path=`cd "$fwd_dir"; pwd` 2>/dev/null
  rev_path=`cd "$rev_dir"; pwd` 2>/dev/null
  nano_path=`cd "$nano_dir"; pwd` 2>/dev/null
# getting filename
  fwd_file=${fwd_reads##*/} 2>/dev/null
  rev_file=${rev_reads##*/} 2>/dev/null
  nano_file=${nano_reads##*/} 2>/dev/null

# Deciding which assembly to use
# nanopore only clonal

if [ -z "${meta}" ]; then
  if [ -z "${fwd_reads}" ]; then
      if [ -z "${nano_reads}" ]; then usage; else wtdbg2_clonal; fi
  fi
fi

# Illumina only clonal
if [ -z "${meta}" ]; then
  if [ -z "${nano_reads}" ]; then
      if [ -z "${fwd_reads}" ]; then usage; else unicycler_illumina_only; fi
  fi
fi
# Hybrid assembly clonal
if [ -z "${meta}" ]; then
  if [ ! -z "${nano_reads}" ]; then
      if [ ! -z "${fwd_reads}" ]; then unicycler_hybrid; fi
  fi
fi
# nanopore only metagenome
if [ ! -z "${meta}" ]; then
  if [ -z "${fwd_reads}" ]; then
    if [ -z "${nano_reads}" ]; then usage; else wtdbg2_meta; fi
  fi
fi
# Illumina only metagenome
if [ ! -z "${meta}" ]; then
  if [ -z "${nano_reads}" ]; then
    if [ -z "${fwd_reads}" ]; then usage; else meta_illumina_only; fi
  fi
fi
# Hybridassembly metagenome
if [ ! -z "${meta}" ]; then
  if [ ! -z "${nano_reads}" ]; then
    if [ ! -z "${fwd_reads}" ]; then meta_hybrid; fi
  fi
fi
