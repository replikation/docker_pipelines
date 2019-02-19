#!/usr/bin/env bash

## Docker ##
  type docker >/dev/null 2>&1 || { echo -e >&2 "${RED}Docker not found.${NC} Please run: ${GRE}sudo apt install docker.io${NC}"; exit 1; }
  # use clonal pipeline by default
    meta=''
    HybMetaAs='metaspades'
    fwd_reads=''
    rev_reads=''
    nano_reads=''
    gsize='4.6m'
    label='results'
  # CPU cores
    CPU=$(lscpu -p | egrep -v '^#' | wc -l) # can be changed to e.g. CPU=16
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
    echo "Wrapper script using docker to choose the best assembly approach for each input combination"
    echo " "
    echo "Usage:    $SCRIPTNAME [-1 illumina_fwd.fastq ] [-2 illumina_rev.fastq] [-n nanopore.fastq] [OPTIONS]"
    echo "          Supports illumina paired only [-1] [-2], nanopore only [-n] and hybrid assembly [-1] [-2] [-n]"
    echo "Inputs:"
    echo -e "          [-1]    ${YEL}Illumina fastq forward reads${NC}; .fastq or .fastq.gz"
    echo -e "          [-2]    ${YEL}Illumina fastq reverse reads${NC}; .fastq or .fastq.gz"
    echo -e "          [-n]    ${YEL}Nanopore fastq file${NC}; .fastq"
    echo " "
    echo "Options:"
    echo -e "          [-t]    Default: ${GRE}-t ${CPU}${NC} - amount of cores"
    echo -e "          [-L]    Default: ${GRE}-L ${gsize}${NC} - approx. genome length"
    echo -e "          [-l]    Label added to output dir e.g. ${GRE}-l Seq_run_2018${NC}"
    echo -e "          [-m]    Input reads are all metagenomic DNA"
    echo "Change Assembler manually (Metagenomic hybrid assembler only):"
    echo -e "          [-o]    Use opera-ms assembler instead"
    echo -e "          [-s]    Special: wtdbg2 assembly first, nanopore contigs polishing via illumina reads after"
    exit;
  }

wtdbg2_clonal()
{
  # WTDBG2 & polishing tested
  echo "Starting wtdbg2 assembly"
  output="wtdgb2_metagenome_${label}"
  mkdir -p $output
    docker run --rm -it --cpus="${CPU}"\
      -v ${nano_path}:/input_nano \
      -v ${WORKDIRPATH}/${output}:/output \
      replikation/wtdbg2_polish \
      wtdbg2 -g $gsize -x ont -L 7000 -t $CPU -i /input_nano/${nano_file} -o /output/draft
    # create contigs
    docker run --rm -it --cpus="${CPU}"\
      -v ${WORKDIRPATH}/${output}:/output \
      replikation/wtdbg2_polish \
      wtpoa-cns -t $CPU -i /output/draft.ctg.lay.gz -fo /output/draft.fa
  # Polishing
  echo "Starting wtdbg2 polishing"
  docker run --rm -it --cpus="${CPU}"\
    -v ${nano_path}:/input_nano \
    -v ${WORKDIRPATH}/${output}:/output \
    replikation/wtdbg2_polish \
    sh -c "minimap2 -t $CPU -x map-pb -a /output/draft.fa /input_nano/${nano_file} | samtools view -Sb - > draft.ctg.map.bam \
    && samtools sort draft.ctg.map.bam > /output/draft.ctg.map.srt.bam \
    && samtools view /output/draft.ctg.map.srt.bam | wtpoa-cns -t $CPU -d /output/draft.fa -i - -fo /output/polished.assembly.fa"
  exit 0
}

wtdbg2_meta()
{
  # WTDBG2 & polishing tested
  echo "Starting wtdbg2 metagenome assembly"
  output="wtdgb2_metagenome_${label}"
  mkdir -p $output
    docker run --rm -it --cpus="${CPU}"\
      -v ${nano_path}:/input_nano \
      -v ${WORKDIRPATH}/${output}:/output \
      replikation/wtdbg2_polish \
      wtdbg2 -p 23 -AS 2 -s 0.05 -e 3 -t $CPU -i /input_nano/${nano_file} -o /output/draft
    #create contigs
    docker run --rm -it --cpus="${CPU}"\
      -v ${WORKDIRPATH}/${output}:/output \
      replikation/wtdbg2_polish \
      wtpoa-cns -t $CPU -i /output/draft.ctg.lay.gz -fo /output/draft.fa
  # Polishing
  echo "Starting wtdbg2 polishing"
  docker run --rm -it --cpus="${CPU}"\
    -v ${nano_path}:/input_nano \
    -v ${WORKDIRPATH}/${output}:/output \
    replikation/wtdbg2_polish \
    sh -c "minimap2 -t $CPU -x map-pb -a /output/draft.fa /input_nano/${nano_file} | samtools view -Sb - > draft.ctg.map.bam \
    && samtools sort draft.ctg.map.bam > /output/draft.ctg.map.srt.bam \
    && samtools view /output/draft.ctg.map.srt.bam | wtpoa-cns -t $CPU -d /output/draft.fa -i - -fo /output/polished.assembly.fa"
  exit 0
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
  output="unicycler_hybrid_${label}"
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
case $HybMetaAs in
metaspades)
  echo "Starting metaspades hybrid assembly"
  output="metaspades_hybrid_${label}"
  mkdir -p $output
  docker run --rm -it \
    -v $fwd_path:/input_fwd \
    -v $rev_path:/input_rev \
    -v $nano_path:/input_nano \
    -v ${WORKDIRPATH}/${output}:/output \
    replikation/spades metaspades.py \
    -1 /input_fwd/${fwd_file} -2 /input_rev/${rev_file} --nanopore /input_nano/${nano_file} -o /output -t $CPU
    exit 0
;;
opera)
  echo "Starting opera-ms hybrid assembly"
  output="opera-ms_hybrid_${label}"
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
;;
wtdbg2_polish)
  echo "Starting wtdbg2 metagenome assembly"
  output="Hybrid-wtdgb2-illumina_metagenome_${label}"
  mkdir -p $output
  docker run --rm -it --cpus="${CPU}"\
    -v ${nano_path}:/input_nano \
    -v ${WORKDIRPATH}/${output}:/output \
    replikation/wtdbg2_polish \
    wtdbg2 -p 23 -AS 2 -s 0.05 -e 3 -t $CPU -i /input_nano/${nano_file} -o /output/draft
  #create contigs
  docker run --rm -it --cpus="${CPU}"\
    -v ${WORKDIRPATH}/${output}:/output \
    replikation/wtdbg2_polish \
    wtpoa-cns -t $CPU -i /output/draft.ctg.lay.gz -fo /output/draft.fa
  echo "Starting unicycler polishing"
  docker run --rm -it --cpus="${CPU}"\
    -v ${fwd_path}:/input_fwd \
    -v ${rev_path}:/input_rev \
    -v $nano_path:/input_nano \
    -v ${WORKDIRPATH}/${output}:/output \
    replikation/unicycler \
    -1 /input_fwd/${fwd_file} -2 /input_rev/${rev_file} -l /input_nano/${nano_file} \
    --existing_long_read_assembly /output/draft.fa  -o /output -t $CPU
;;
esac
}

#############################
###   Start of script    ####
#############################
echo "                                               _____________________________"
echo "______________________________________________/ Created by Christian Brandt \___"
echo " "

# you could add a output flag
while getopts '1:2:n:mt:L:l:osh' flag; do
    case "${flag}" in
      1) fwd_reads="${OPTARG}" ;;
      2) rev_reads="${OPTARG}" ;;
      n) nano_reads="${OPTARG}" ;;
      m) meta='true' ;;
      l) label="${OPTARG}" ;;
      t) CPU="${OPTARG}" ;;
      L) gsize="${OPTARG}" ;;
      o) HybMetaAs='opera';;
      s) HybMetaAs='wtdbg2_polish';;
      h) usage;;
      *) usage
         exit 1 ;;
    esac
done

# getting dir names
  fwd_dir=$(dirname "$fwd_reads" 2>/dev/null)
  rev_dir=$(dirname "$rev_reads" 2>/dev/null)
  nano_dir=$(dirname "$nano_reads" 2>/dev/null)
# getting absolute paths
  fwd_path=$(cd "$fwd_dir" 2>/dev/null && pwd)
  rev_path=$(cd "$rev_dir" 2>/dev/null && pwd)
  nano_path=$(cd "$nano_dir" 2>/dev/null && pwd)
# getting filename
  fwd_file=${fwd_reads##*/}
  rev_file=${rev_reads##*/}
  nano_file=${nano_reads##*/}

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
