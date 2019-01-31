#!/usr/bin/env bash

## Docker ##
  type docker >/dev/null 2>&1 || { echo -e >&2 "${RED}Docker not found.${NC} Please run: ${GRE}sudo apt install docker.io${NC}"; exit 1; }
  # use clonal pipeline by default
    meta=''
    opera=''
    fwd_reads=''
    rev_reads=''
    nano_reads=''
    L_wtdbg2='4.6m'
  # CPU cores
    CPU=$(lscpu -p | egrep -v '^#' | wc -l) # can be changed to e.g. CPU="16"
  # colours, needed for echos
    RED='\033[0;31m'
    YEL='\033[0;33m'
    NC='\033[0m'
    GRE='\033[0;32m'
    ## Parameters ##
      WORKDIRPATH=$(pwd) # for docker mountpoint (-v)
      SCRIPTNAME=$(basename -- "$0")

###############
##  Modules  ##
###############

usage()
  {
    echo "Usage:    $SCRIPTNAME [-1 illumina_fwd.fastq ] [-2 illumina_rev.fastq] [-n nanopore.fastq] [OPTIONS]"
    echo "Supports illumina paired only [-1] [-2], nanopore only [-n] and hybrid assembly [-1] [-2] [-n]"
    echo "Inputs:"
    echo -e "          [-1]    ${YEL}Illumina fastq forward reads${NC}; .fastq or .fastq.gz"
    echo -e "          [-2]    ${YEL}Illumina fastq reverse reads${NC}; .fastq or .fastq.gz"
    echo -e "          [-n]    ${YEL}Nanopore fastq file${NC}; .fastq"
    echo "Options:"
    echo -e "          [-t]    Default: ${GRE}-t ${CPU}${NC} - amount of cores"
    echo -e "          [-L]    Default: ${GRE}-L ${L_wtdbg2}${NC} - approx. genome length"
    echo "Metagenome options"
    echo -e "          [-m]    Default: Single organism - Add ${GRE}-m${NC} for metagenome sample"
    echo -e "          [-o]    Default: metaspades - add ${GRE}-o ${NC}to use opera-ms instead"
    exit;
  }

wtdbg2_clonal()
{
  # WTDBG2 - untested
  echo "Starting wtdbg2 assembly"
  output="wtdgb2_assembly"
  mkdir -p $output
    docker run --rm -it \
      -v $nano_path:/input_nano \
      -v ${WORKDIRPATH}/${output}:/output \
      replikation/wtdbg2 \
      wtdbg2 -x ont -t $CPU -g $L_wtdbg2 -i /input_nano/${nano_file} -o /output/${nano_file%.fastq}
    #create contigs
    docker run --rm -it \
      -v ${WORKDIRPATH}/${output}:/output \
      replikation/wtdbg2 \
      wtpoa-cns -t $CPU -i /output/${nano_file%.fastq}.ctg.lay.gz \
      -fo /output/${nano_file%.fastq}.fa
exit 0
}

wtdbg2_meta()
{
  # WTDBG2 - untested
  echo "canu is actually best for metagenome assembly"
}

unicycler_illumina_only()
{
  # UNICYCLER - untested
  echo "Starting unicycler assembly"
  output="unicycler_assembly"
  mkdir -p $output
  docker run --rm -it \
    -v $fwd_path:/input_fwd \
    -v $rev_path:/input_rev \
    -v ${WORKDIRPATH}/${output}:/output \
    replikation/unicycler \
    -1 /input_fwd/${fwd_file} -2 /input_rev/${rev_file} -o /output -t $CPU
    exit 0
}

unicycler_hybrid()
{
  # UNICYCLER - untested
  echo "Starting unicycler hybrid assembly"
  output="unicycler_hybrid_assembly"
  mkdir -p $output
  docker run --rm -it \
    -v ${fwd_path}:/input_fwd \
    -v ${rev_path}:/input_rev \
    -v $nano_path:/input_nano \
    -v ${WORKDIRPATH}/${output}:/output \
    replikation/unicycler \
    -1 /input_fwd/${fwd_file} -2 /input_rev/${rev_file} -l /input_nano/${nano_file} -o /output -t $CPU
    exit 0
}

meta_illumina_only()
{
  # METASPADES - untested
  echo "Starting metaspades assembly"
  output="metaspades_assembly"
  mkdir -p $output
  docker run --rm -it \
    -v ${fwd_path}:/input_fwd \
    -v ${rev_path}:/input_rev \
    -v ${WORKDIRPATH}/${output}:/output \
    replikation/spades metaspades.py \
    -1 /input_fwd/${fwd_file} -2 /input_rev/${rev_file} -o /output -t ${CPU}
    exit 0
}

meta_hybrid_assembly()
{
if [ -z "${opera}" ]; then
  # METASPADES - tested
  echo "Starting metaspades hybrid assembly"
  output="metaspades_assembly"
  mkdir -p $output
  docker run --rm -it \
    -v $fwd_path:/input_fwd \
    -v $rev_path:/input_rev \
    -v $nano_path:/input_nano \
    -v ${WORKDIRPATH}/${output}:/output \
    replikation/spades metaspades.py \
    -1 /input_fwd/${fwd_file} -2 /input_rev/${rev_file} --nanopore /input_nano/${nano_file} -o /output -t $CPU
    exit 0
else
  # OPERA-MS - tested
  echo "Starting opera-ms hybrid assembly"
  output="opera-ms_assembly"
  mkdir -p $output
  # unzip illumina if .gz - if .fastq nothing happens
  gunzip $fwd_reads 2>/dev/null
  gunzip $rev_reads 2>/dev/null
  # create config file
  mkdir -p config
  echo "OUTPUT_DIR /output" > config/config.file
  echo "ILLUMINA_READ_1 /input_fwd/${fwd_file%.gz}" >>  config/config.file
  echo "ILLUMINA_READ_2 /input_rev/${rev_file%.gz}" >>  config/config.file
  echo "LONG_READ /input_nano/$nano_file" >>  config/config.file
  echo "NUM_PROCESSOR $CPU" >>  config/config.file
  echo "STRAIN_CLUSTERING YES" >>  config/config.file
  echo "CONTIG_LEN_THR 500" >>  config/config.file
  echo "CONTIG_EDGE_LEN 80" >>  config/config.file
  echo "CONTIG_WINDOW_LEN 340" >>  config/config.file
  echo "KMER_SIZE 60" >>  config/config.file
  echo "LONG_READ_MAPPER blasr" >>  config/config.file
  #echo "CONTIGS_FILE sample_files/sample_contigs.fasta"
  # run assembly
  docker run --rm -it \
    -v $fwd_path:/input_fwd \
    -v $rev_path:/input_rev \
    -v $nano_path:/input_nano \
    -v ${WORKDIRPATH}/${output}:/output \
    -v ${WORKDIRPATH}/config/:/config \
    replikation/opera_ms /config/config.file
  rm -fr config/
  exit 0
fi
}

#############################
###   Start of script    ####
#############################
echo "                                               _____________________________"
echo "______________________________________________/ Created by Christian Brandt \___"
echo " "

# you could add a output flag
while getopts '1:2:n:mt:L:o' flag; do
    case "${flag}" in
      1) fwd_reads="${OPTARG}" ;;
      2) rev_reads="${OPTARG}" ;;
      n) nano_reads="${OPTARG}" ;;
      m) meta='true' ;;
      t) CPU="${OPTARG}" ;;
      L) L_wtdbg2="${OPTARG}" ;;
      o) opera='true';;
      *) usage
         exit 1 ;;
    esac
done

# getting dir names
  fwd_dir=$(dirname "$fwd_reads") 2>/dev/null
  rev_dir=$(dirname "$rev_reads") 2>/dev/null
  nano_dir=$(dirname "$nano_reads") 2>/dev/null
# getting absolute paths
  fwd_path=$(cd $fwd_dir && pwd) 2>/dev/null
  rev_path=$(cd "$rev_dir" && pwd) 2>/dev/null
  nano_path=$(cd "$nano_dir" && pwd) 2>/dev/null
# getting filename
  fwd_file=${fwd_reads##*/} 2>/dev/null
  rev_file=${rev_reads##*/} 2>/dev/null
  nano_file=${nano_reads##*/} 2>/dev/null
echo " "
# Deciding which assembly aproach to use
## nanopore only clonal
if [ -z "${meta}" ]; then
  if [ -z "${fwd_reads}" ]; then
      if [ -z "${nano_reads}" ]; then usage; else wtdbg2_clonal; fi
  fi
fi
## Illumina only clonal
if [ -z "${meta}" ]; then
  if [ -z "${nano_reads}" ]; then
      if [ -z "${fwd_reads}" ]; then usage; else unicycler_illumina_only; fi
  fi
fi
## Hybrid assembly clonal
if [ -z "${meta}" ]; then
  if [ ! -z "${nano_reads}" ]; then
      if [ ! -z "${fwd_reads}" ]; then unicycler_hybrid; fi
  fi
fi
## nanopore only metagenome
if [ ! -z "${meta}" ]; then
  if [ -z "${fwd_reads}" ]; then
    if [ -z "${nano_reads}" ]; then usage; else wtdbg2_meta; fi
  fi
fi
## Illumina only metagenome
if [ ! -z "${meta}" ]; then
  if [ -z "${nano_reads}" ]; then
    if [ -z "${fwd_reads}" ]; then usage; else meta_illumina_only; fi
  fi
fi
## Hybridassembly metagenome
if [ ! -z "${meta}" ]; then
  if [ ! -z "${nano_reads}" ]; then
    if [ ! -z "${fwd_reads}" ]; then meta_hybrid_assembly; fi
  fi
fi
