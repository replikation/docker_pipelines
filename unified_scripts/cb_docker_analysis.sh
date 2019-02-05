#!/usr/bin/env bash

## Docker ##
  type docker >/dev/null 2>&1 || { echo -e >&2 "${RED}Docker not found.${NC} Please run: ${GRE}sudo apt install docker.io${NC}"; exit 1; }
  # use clonal pipeline by default
    error=''
    assembly_file=''
    bam_file=''
    nano_reads=''
    sour_DB=''
    CPU=''
    binning=''
    tax_read=''
    tax_assembly=''
    nanoQC=''
    centrif_DB=''
    label='results'
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
    echo "Usage:    $SCRIPTNAME [INPUT(S)]  [OPTIONS] [ANALYSIS OPTIONS]"
    echo "e.g.:     $SCRIPTNAME -a assembly.fa -b bamfile.bam -t 20 -l first_assembly -B"
    echo " "
    echo "Inputs:"
    echo -e "          [-a]    ${YEL}Assembly file${NC}; .fa or .fasta"
    echo -e "          [-b]    ${YEL}bam file${NC}; .bam"
    echo -e "          [-n]    ${YEL}Nanopore read file${NC}; .fastq"
    echo -e "          [-s]    ${YEL}sequencing_summary.txt location${NC};name must be: sequencing_summary.txt"
    echo -e "          [-D]    ${YEL}Sourmash database location${NC}; .sbt.json"
    echo "Options:"
    echo -e "          [-t]    Default: ${GRE}-t ${CPU}${NC} - amount of cores"
    echo -e "          [-l]    Label added to output dir e.g. ${GRE}-l Seq_run_2018${NC}"
    echo "Analysis options:"
    echo "  Add flag(s) to use one or more of these"
    echo -e "          [-B]    vamb, metagenome binning of contig file; Input: ${YEL}-a -b${NC}"
    echo -e "          [-L]    centrifuge, taxonomic classification of Long Reads (bacteria & archaea); Input: ${YEL}-n${NC}"
    echo -e "              [-k] to use bacteria, viruses, human and archaea database instead "
    echo -e "          [-C]    sourmash, taxonomic classification on contigs; Input: ${YEL}-a -D${NC}"
    echo -e "          [-Q]    nanopore QC, QC results for reads; Input: ${YEL}-s${NC}"
    exit;
  }

centrifuge_execute()
{
  # tested
  if [ -z "${centrif_DB}" ]; then
  echo "Starting centrifuge for bacteria and archaea"
  output="centrifuge_bac.arch_${label}"
  mkdir -p $output
  docker run --rm -i \
    -v $nano_path:/input \
    -v $WORKDIRPATH/${output}:/output \
    replikation/centrifuge \
    centrifuge -p $CPU -x /centrifuge/database/p_compressed -k 1 --min-hitlen 16 \
    -U /input/$nano_name -S /output/centrifuge_out.txt --report-file /output/centrifuge_out.log
    # create report for pavian
    docker run --rm -i -v $WORKDIRPATH/${output}:/output replikation/centrifuge \
    centrifuge-kreport -x /centrifuge/database/p_compressed --min-score 300 --min-length 500 /output/centrifuge_out.txt \
    > $WORKDIRPATH/${output}/${nano_name%.*}_pavian_report.csv
  echo "Results saved to $output"
else
  echo "Starting centrifuge for bacteria, viruses, human and archaea"
  output="centrifuge_bac.arch.vir.hum_${label}"
  mkdir -p $output
  docker run --rm -i \
    -v $nano_path:/input \
    -v $WORKDIRPATH/${output}:/output \
    replikation/centrifuge_all \
    centrifuge -p $CPU -x /centrifuge/database/p_compressed+h+v -k 1 --min-hitlen 16 \
    -U /input/$nano_name -S /output/centrifuge_out.txt --report-file /output/centrifuge_out.log
    # create report for pavian
    docker run --rm -i -v $WORKDIRPATH/${output}:/output replikation/centrifuge_all \
    centrifuge-kreport -x /centrifuge/database/p_compressed+h+v --min-score 300 --min-length 500 /output/centrifuge_out.txt \
    > $WORKDIRPATH/${output}/${nano_name%.*}_pavian_report.csv
  echo "Results saved to $output"
fi
}

sourmash_execute()
{
  # untested
  echo "Starting sourmash"
  output="sourmash_${label}"
  mkdir -p $output
  # create signature
    docker run --rm -it \
      -v $WORKDIRPATH/${output}:/output \
      -v $assembly_path:/input \
      replikation/sourmash \
    /bin/sh -c "cd /input/ && sourmash compute -n 5000 -k 31 $assembly_name -o /output/contig_signatures.sig"
    # compare
    docker run --rm -it \
      -v $WORKDIRPATH/${output}:/output \
      -v $sour_path:/DB_sour \
      replikation/sourmash \
      sourmash gather -k 31 /output/contig_signatures.sig /DB_sour/$DB_name \
      > $WORKDIRPATH/${output}/${assembly_name%.*}_results.tsv
    echo "Results saved to $output"
}

binning_execute()
{
  echo "vamb"
  output="vamb_${label}"
  echo "-v $assembly_path:/input"
  echo "-v $bam_path:/bam_files"
  echo "-v $WORKDIRPATH/$output:/output"
}

QC_nanopore()
{
  # tested
  echo "Starting nanopore QC"
  output="nanoporeQC_${label}"
  mkdir -p $output
  docker run --rm -it --cpus="${CPU}" \
    -v ${WORKDIRPATH}/${output}:/output \
    -v ${seqSum_path}:/QCTutorial/RawData \
    replikation/nanopore_qc \
    R --slave -e 'rmarkdown::render("Nanopore_SumStatQC_Tutorial.Rmd", output_format = "html_document", output_file = "QC_reads.html", output_dir = "/output/")'
  echo "Results saved to $output"
}

#############################
###   Start of script    ####
#############################
echo "                                               _____________________________"
echo "______________________________________________/ Created by Christian Brandt \___"
echo " "
echo -e "${YEL}$SCRIPTNAME -h ${NC}for help/usage"
echo " "
# you could add a output flag
while getopts 'a:b:n:s:D:t:l:BLCQkh' flag; do
    case "${flag}" in
      a) assembly_file="${OPTARG}" ;;
      b) bam_file="${OPTARG}" ;;
      n) nano_reads="${OPTARG}" ;;
      s) seqSum="${OPTARG}" ;;
      D) sour_DB="${OPTARG}" ;;
      t) CPU="${OPTARG}" ;;
      l) label="${OPTARG}" ;;
      B) binning='true' ;;
      L) tax_read='true';;
      C) tax_assembly='true';;
      Q) nanoQC='true';;
      k) centrif_DB='true';;
      h) usage;;
      *) usage
         exit 1 ;;
    esac
done

# getting dir names
  assembly_dir=$(dirname "$assembly_file") 2>/dev/null
  bam_dir=$(dirname "$bam_file") 2>/dev/null
  nano_dir=$(dirname "$nano_reads") 2>/dev/null
  seqSum_dir=$(dirname "$seqSum") 2>/dev/null
  sour_dir=$(dirname "$sour_DB") 2>/dev/null
# getting absolute paths
  assembly_path=$(cd $assembly_dir && pwd) 2>/dev/null
  bam_path=$(cd "$bam_dir" && pwd) 2>/dev/null
  nano_path=$(cd "$nano_dir" && pwd) 2>/dev/null
  seqSum_path=$(cd "$seqSum_dir" && pwd) 2>/dev/null
  sour_path=$(cd "$sour_dir" && pwd) 2>/dev/null
# getting filename
  assembly_name=${assembly_file##*/} 2>/dev/null
  nano_name=${nano_reads##*/} 2>/dev/null
  DB_name=${sour_DB##*/} 2>/dev/null
echo " "
#############################
## Choose Executable(s)    ##
#############################

# Taxonomic read classification
  if [ ! -z "${tax_read}" ]; then
      if [ ! -z "${nano_reads}" ]; then centrifuge_execute ; else error=true ; fi
  fi
# Taxonomic assembly classification
  if [ ! -z "${tax_assembly}" ]; then
      if [ ! -z "${assembly_file}" ]; then
        if [ ! -z "${sour_DB}" ]; then sourmash_execute ; else error=true ; fi
      else error=true ; fi
  fi
# Taxonomic assembly classification
  if [ ! -z "${binning}" ]; then
      if [ ! -z "${assembly_file}" ]; then
        if [ ! -z "${bam_file}" ]; then binning_execute ; else error=true ; fi
      else error=true ; fi
  fi
# Nanopore QC
  if [ ! -z "${nanoQC}" ]; then
    if [ ! -z "${seqSum}" ]; then QC_nanopore ; else error=true ; fi
  fi

# error filled
  if [ ! -z "${error}" ]; then usage ; fi
