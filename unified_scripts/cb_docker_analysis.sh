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
    echo "Usage:    $SCRIPTNAME [-a assembly.fa ] [-b file.bam] [-n nanopore reads] [-D path_to_database] [OPTIONS]"
    echo " "
    echo "Inputs:"
    echo -e "          [-a]    ${YEL}Assembly file${NC}; .fa or .fasta"
    echo -e "          [-b]    ${YEL}bam file${NC}; .bam"
    echo -e "          [-n]    ${YEL}Nanopore read file${NC}; .fastq"
    # echo -e "          [-s]    ${YEL}sequencing_summary.txt location${NC};name must be sequencing_summary.txt"
    echo -e "          [-D]    ${YEL}Sourmash database location${NC}; .sbt.json"
    echo "Options:"
    echo -e "          [-t]    Default: ${GRE}-t ${CPU}${NC} - amount of cores"
    echo "Analysis options - add flag(s) to use one or more of these"
    echo -e "          [-B]    vamb, metagenome binning of contig file; needs ${YEL}-a -b${NC}"
    echo -e "          [-L]    centrifuge, taxonomic classification of Long Reads; needs ${YEL}-n${NC}"
    echo -e "          [-C]    sourmash, taxonomic classification on contigs; needs ${YEL}-a -D${NC}"
    exit;
  }

centrifuge_execute()
{
  # untested
  echo "Starting centrifuge"
  output="centrifuge_${assembly_name%.*}"
  mkdir $output
  docker run --rm -i \
    -v $nano_path:/input \
    -v $WORKDIRPATH/${output}:/output \
    replikation/centrifuge \
    centrifuge -p $CPU -x /centrifuge/database/p_compressed \
    -U /input/$nano_name -S /output/centrifuge_out.txt --report-file /output/centrifuge_out.log
    # create report for pavian
    docker run --rm -i -v $WORKDIRPATH/${output}:/output replikation/centrifuge \
    centrifuge-kreport -x /centrifuge/database/p_compressed /output/centrifuge_out.txt \
    > $WORKDIRPATH/${output}/${nano_name%.*}_pavian_report.csv
}

sourmash_execute()
{
  # untested
  echo "Starting sourmash"
  output="sourmash_${assembly_name%.*}"
  mkdir $output
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
}

binning_execute()
{
  echo "vamb"
  output="vamb_${assembly_name%.*}"
  echo "-v $assembly_path:/input"
  echo "-v $bam_path:/bam_files"
  echo "-v $WORKDIRPATH/$output:/output"
}

QC_nanopore()
{
  # mount seq.summary  directly to the /QCTutorial/RawData of the git inside
  config=config
  mkdir $config
  docker run --rm -it \
    -v $PWD:/output \ #change to your output
    -v $PWD/${config}:/QCTutorial/RawData \ #this will be the sequencing summary location via flag
    replikation/nanopore_qc \
    R --slave -e 'rmarkdown::render("Nanopore_SumStatQC_Tutorial.Rmd", output_format = "html_document", output_dir = "/output/")'
    # für render gibts noch einen "outputnamen" - damit das ding nicht tutorial heißt"
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
while getopts 'a:b:n:D:t:BLCh' flag; do
    case "${flag}" in
      a) assembly_file="${OPTARG}" ;;
      b) bam_file="${OPTARG}" ;;
      n) nano_reads="${OPTARG}" ;;
      D) sour_DB="${OPTARG}" ;;
      t) CPU="${OPTARG}" ;;
      B) binning='true' ;;
      L) tax_read='true';;
      C) tax_assembly='true';;
      h) usage;;
      *) usage
         exit 1 ;;
    esac
done

# getting dir names
  assembly_dir=$(dirname "$assembly_file") 2>/dev/null
  bam_dir=$(dirname "$bam_file") 2>/dev/null
  nano_dir=$(dirname "$nano_reads") 2>/dev/null
  sour_dir=$(dirname "$sour_DB") 2>/dev/null
# getting absolute paths
  assembly_path=$(cd $assembly_dir && pwd) 2>/dev/null
  bam_path=$(cd "$bam_dir" && pwd) 2>/dev/null
  nano_path=$(cd "$nano_dir" && pwd) 2>/dev/null
  sour_path=$(cd "$sour_dir" && pwd) 2>/dev/null
# getting filename
  assembly_name=${assembly_file##*/} 2>/dev/null
  nano_name=${nano_reads##*/} 2>/dev/null
  DB_name=${sour_DB##*/} 2>/dev/null
echo " "

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

# error filled
if [ ! -z "${error}" ]; then usage ; fi
