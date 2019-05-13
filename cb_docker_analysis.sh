#!/usr/bin/env bash

## Docker ##
  type docker >/dev/null 2>&1 || { echo -e >&2 "${RED}Docker not found.${NC} Please run: ${GRE}sudo apt install docker.io${NC}"; exit 1; }
## clear all variables - to be sure - paranoid ftw
  unset assembly_file plasflow nano_reads fwd_reads rev_reads input_folder DB_in
  unset CPU binning tax_read sour_cluster nanoQC centrif_DB AB_res recentrifuge deepvir
## colors, needed for echos
  RED='\033[0;31m'
  YEL='\033[0;33m'
  BLU='\033[0;94m'
  NC='\033[0m'
  GRE='\033[0;32m'
  DIM='\033[1;30m'
## Parameters ##
  WORKDIRPATH=$(pwd) # for docker mountpoint (-v)
  SCRIPTNAME=$(basename -- "$0")
  # default for -t $CPU cores
  CPU=$(lscpu -p | egrep -v '^#' | wc -l)
  # default for -l
  label='results'

###############
##  Modules  ##
###############

usage()
  {
    echo "Usage:    $SCRIPTNAME [INPUT(S)]  [OPTION(S)] [ANALYSIS TOOL(S)]"
    echo -e "${DIM}e.g.:     $SCRIPTNAME -n reads.fastq -t 20 -l first_assembly -Lk${NC}"
    echo " "
    echo "Inputs:"
    echo -e "          [-a]    ${YEL}Assembly or contig file${NC}; .fa or .fasta"
    echo -e "          [-n]    ${YEL}Nanopore read file${NC}; .fastq"
    echo -e "          [-1]    ${YEL}Illumina forward read file${NC}; .fastq or .fastq.gz"
    echo -e "          [-2]    ${YEL}Illumina reverse read file${NC}; .fastq or .fastq.gz"
    echo -e "          [-s]    ${YEL}sequencing_summary.txt file${NC}; name must be: sequencing_summary.txt"
    echo -e "          [-f]    ${YEL}directory location${NC}; e.g. centrifuge_results/"
    echo "Options:"
    echo -e "          [-t]    Default: ${GRE}-t ${CPU}${NC} - amount of cores"
    echo -e "          [-l]    Label added to output dir e.g. ${GRE}-l Seq_run_2018${NC}"
    echo "Analysis tools:"
    echo -e "${DIM}Metagenome${NC}"
    echo -e "          [-B]    ${BLU}B${NC}inning via metabat-checkM of contigs; Input: ${YEL}-1 -2 -a${NC}"
    echo -e "          [-C]    ${BLU}C${NC}entrifuge, tax. classif. of reads (bacteria & archaea); Input: ${YEL}-n${NC} or ${YEL}-1 -2${NC}"
    echo -e "              [-c]   custom DB, absolut path to database w/o suffixes like .1.cf (${YEL}-n${NC} only)"
    echo -e "          [-K]    ${BLU}K${NC}rona-recentrifuge; Input: ${YEL}-f${NC}, contains atleast 1 *.out file(s) from [-C]"
    echo -e "${DIM}QC${NC}"
    echo -e "          [-Q]    ${BLU}Q${NC}C for nanopore reads, QC results for reads; Input: ${YEL}-s${NC}"
    echo -e "${DIM}Screening${NC}"
    echo -e "          [-P]    ${BLU}P${NC}lasflow, plasmid binning; Input: ${YEL}-a${NC}"
    echo -e "          [-S]    ${BLU}S${NC}ourmash, sequence clustering; Input: ${YEL}-a${NC}"
    echo -e "          [-R]    ${BLU}R${NC}esistance gene screening via ABRicate; Input: ${YEL}-a${NC} or ${YEL}-n${NC}"
    echo -e "          [-D]    ${BLU}D${NC}eepVirFinder, predicts viral sequences; Input: ${YEL}-a${NC}"
    exit;
  }

centrifuge_nanopore()
{
  # Standard parameters
  dockerimage_centri='centrifuge'
  DB_default='/centrifuge/database/p_compressed'
  tag_centri='bacteria_archeae'
  unset DB_in
  # changing default parameters for custom database
    if [ ! -z "${centrif_DB}" ]; then
      dockerimage_centri='centrifuge_small'
      DB_default="/DB_c/$centrif_DB_file"
      DB_in="-v ${centrif_DB_path}:/DB_c"
      tag_centri='custom_DB'
    fi
  # running centrifuge
  echo "Starting centrifuge for ${tag_centri}"
  output="centrifuge_nanopore_${tag_centri}_${label}"
  mkdir -p $output
  docker run --rm -i --cpus="${CPU}" ${DB_in} \
    -v $nano_path:/input \
    -v $WORKDIRPATH/${output}:/output \
    replikation/$dockerimage_centri \
    centrifuge -p $CPU -x $DB_default -k 5 --min-hitlen 16 \
    -U /input/$nano_name -S /output/centrifuge_results.out --report-file /output/centrifuge_out.log
  # filter reads by $4 = score and $6 = min hit length
  < ${output}/centrifuge_results.out awk '{if(NR < 2 || $4 >= 250) {print}}' | awk '{if(NR < 2 || $6 >= 150) {print}}' \
  > ${output}/centrifuge_filtered.out
  # create report for pavian
  docker run --rm -i --cpus="${CPU}" ${DB_in} -v $WORKDIRPATH/${output}:/output replikation/$dockerimage_centri \
    centrifuge-kreport -x $DB_default --min-score 300 --min-length 500 /output/centrifuge_filtered.out \
    > $WORKDIRPATH/${output}/${nano_name%.*}_pavian_report_filtered.csv
  echo "Results saved to $output"
}

centrifuge_illumina()
{
  echo "Starting centrifuge for bacteria and archaea (illumina)"
  output="centrifuge_illumina_bac.arch_${label}"
  mkdir -p $output
  docker run --rm -i --cpus="${CPU}"\
    -v $fwd_path:/input_fwd \
    -v $rev_path:/input_rev \
    -v $WORKDIRPATH/${output}:/output \
    replikation/centrifuge \
    centrifuge -p $CPU -x /centrifuge/database/p_compressed -k 5 --min-hitlen 22 \
    -1 /input_fwd/${fwd_file} -2 /input_rev/${rev_file} -S /output/centrifuge_results.out --report-file /output/centrifuge_out.log
    # create report for pavian
    docker run --rm -i -v $WORKDIRPATH/${output}:/output replikation/centrifuge \
    centrifuge-kreport -x /centrifuge/database/p_compressed --min-score 300 --min-length 500 /output/centrifuge_results.out \
    > $WORKDIRPATH/${output}/${fwd_file%.*}_pavian_report.csv
  echo "Results saved to $output"
}

recentrifuge()
{
  output="recentrifuge_${label}"
  mkdir -p $output
  echo -e "Searching for *.out files in ${YEL}${infolder_path}${NC} ..."
    test_files=$(ls -1 ${infolder_path}/*.out 2> /dev/null)
    if [ -z "${test_files}" ]; then echo -e "  Can't find .out files in ${YEL}${infolder_path}${NC}, exiting"; exit 1; fi
    # adjusting command for more than 1 file
      filenames=$(find ${infolder_path}/*.out -type f -print0  | xargs -0 -n 1 basename)
      input_for_docker="${filenames//$'\n'/ -f /input/}"
      input_for_docker='-f /input/'${input_for_docker}
      Num_of_samples=$(echo "$filenames" | wc -l)
  echo -e "Found ${YEL}${Num_of_samples}${NC} file(s)"
  docker run --rm -it \
    -v ${infolder_path}:/input \
    -v $WORKDIRPATH/${output}:/output \
    replikation/recentrifuge \
    -n /database/ncbi_node $input_for_docker -o /output/${Num_of_samples}_samples_overview.html
}

plasflow_execute()
{
  echo "Starting plasflow, removing contigs below 2000 bp"
  output="plasflow_${label}"
  docker run --rm -it --cpus="${CPU}"\
    -v $WORKDIRPATH/${output}:/output \
    -v $assembly_path:/input \
    replikation/plasflow \
    filter_sequences_by_length.pl -input /input/${assembly_name} \
    -output /output/${assembly_name} -thresh 2000
  docker run --rm -it --cpus="${CPU}"\
    -v $WORKDIRPATH/${output}:/output \
    replikation/plasflow \
    PlasFlow.py --input /output/${assembly_name} \
    --output /output/plasflow_predictions --threshold 0.7
  echo "Results saved to $output"
}

QC_nanopore()
{
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

sourmash_cluster()
{
  echo "Starting sourmash clustering"
  output="sourmash_cluster_${label}"
  mkdir -p $output
  # remove .fasta from file for better picture?
  docker run --rm -it --cpus="${CPU}"\
    -v $WORKDIRPATH/${output}:/output \
    -v $assembly_path:/input \
    replikation/sourmash \
    /bin/sh -c "cd /input/ && sourmash compute --scaled 1000 -k 31 $assembly_name -o /output/signature.sig --singleton"
  docker run --rm -it -v $WORKDIRPATH/${output}:/output replikation/sourmash \
    sourmash compare /output/signature.sig -o /output/results_sig
  docker run --rm -it -v $WORKDIRPATH/${output}:/output replikation/sourmash \
    /bin/sh -c "sourmash plot --pdf --subsample=250 --labels /output/results_sig && mv *.pdf /output/"
  echo -e "${RED}Hint: Won't plot more than 250 samples${NC}"
  }

resistance_screen()
{
  echo "Starting ABRicate resistance gene screening"
  docker pull replikation/abricate
  output="ABRicate_resistances_${label}"
  mkdir -p $output
  DB_list=$(echo -e "resfinder\nncbi\ncard\nplasmidfinder\nargannot")
  # parallel
  echo "$DB_list" | xargs -I% -P ${CPU} \
    sh -c "docker run --rm \
    -v $assembly_path:/input \
    replikation/abricate \
    /input/$assembly_name --nopath --quiet --mincov 25 --db % > ${WORKDIRPATH}/${output}/%.tab "
}

resistance_read_screen()
{
  echo "Starting ABRicate resistance gene screening against reads greater 1000 bp"
  docker pull replikation/abricate
  output="ABRicate_resistances_vs_reads_${label}"
  mkdir -p $output
  echo "Preparing reads first (takes time)..."
  sed -n '1~4s/^@/>/p;2~4p' $nano_reads | \
  sed ':a;N;/^>/M!s/\n//;ta;P;D' | \
  awk '/^>/ { getline seq } length(seq) >1000 { print $0 "\n" seq }' > $output/filtered_reads.fasta
  echo "Preparing done."
  DB_list=$(echo -e "resfinder\nncbi\ncard\nplasmidfinder\nargannot")
  # parallel
  echo "$DB_list" | xargs -I% -P ${CPU} \
    sh -c "docker run --rm \
    -v ${WORKDIRPATH}/${output}:/input \
    replikation/abricate \
    /input/filtered_reads.fasta --nopath --quiet --mincov 25 --db % > ${WORKDIRPATH}/${output}/%.tab "
  rm -f $output/filtered_reads.fasta
}

binning_execute()
{
  # untested
  echo "Starting binning with metabat"
  output="metabat_${label}"
  mkdir -p $output
  # read mapping
  docker run --rm -it --cpus="${CPU}"\
    -v $fwd_path:/input_fwd \
    -v $rev_path:/input_rev \
    -v $assembly_path:/input \
    -v $WORKDIRPATH/${output}:/output \
    --entrypoint "/bin/sh" \
    replikation/unicycler \
    -c "bowtie2-build /input/${assembly_name} /output/assembly.contigs && \
    bowtie2 -p ${CPU} -x /output/assembly.contigs -1 /input_fwd/$fwd_file -2 /input_rev/$rev_file | samtools view -bS -o binary.bam && \
    samtools sort binary.bam -o /output/${assembly_name}.srt.bam && \
    samtools index /output/${assembly_name}.srt.bam"
  # binning
  docker run --rm -it --cpus="${CPU}"\
    -v $assembly_path:/input \
    -v $WORKDIRPATH/${output}:/output \
    metabat/metabat \
    sh -c "cp /input/${assembly_name} /output/ && cd /output/ && runMetaBat.sh -m 1500 ${assembly_name} ${assembly_name}.srt.bam \
    && mv *metabat-bins1500 metabat_bins"
  # bin integrety & plotting
  docker run --rm -it --cpus="${CPU}"\
     -v $WORKDIRPATH/${output}:/output \
     replikation/checkm \
     lineage_wf -x fa /output/metabat_bins -t ${CPU} /output/checkm/
  docker run --rm -it --cpus="${CPU}"\
     -v $WORKDIRPATH/${output}:/output \
     replikation/checkm \
     bin_qa_plot -x fa /output/checkm/ /output/metabat_bins /output/plots
}

deepvirfinder_excecute()
{
  echo "Starting deepvirfinder predictions"
  output="deepvirfinder_${label}"
  mkdir -p $output
  docker run --rm -it --cpus="${CPU}"\
  -v $WORKDIRPATH/${output}:/output \
  -v $assembly_path:/input \
  replikation/deepvirfinder \
  -c ${CPU} -i /input/${assembly_name} -o /output/
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
while getopts 'a:1:2:n:s:f:t:l:BCc:KSQPDRh' flag; do
    case "${flag}" in
      a) assembly_file="${OPTARG}" ;;
      1) fwd_reads="${OPTARG}" ;;
      2) rev_reads="${OPTARG}" ;;
      n) nano_reads="${OPTARG}" ;;
      s) seqSum="${OPTARG}" ;;
      f) input_folder="${OPTARG}" ;;
        t) CPU="${OPTARG}" ;;
        l) label="${OPTARG}" ;;
      B) binning='true' ;;
      C) tax_read='true';;
        c) centrif_DB="${OPTARG}";;
      K) recentrifuge='true';;
      S) sour_cluster='true';;
      Q) nanoQC='true';;
      P) plasflow='true';;
      D) deepvir='true';;
      R) AB_res='true';;
        h) usage;;
        *) usage
         exit 1 ;;
    esac
done

# Paths for Docker compatibility
# getting dir names from file inputs
  if [ ! -z "${assembly_file}" ]; then assembly_dir=$(dirname "$assembly_file"); fi
  if [ ! -z "${nano_reads}" ]; then nano_dir=$(dirname "$nano_reads"); fi
  if [ ! -z "${fwd_reads}" ]; then fwd_dir=$(dirname "$fwd_reads"); fi
  if [ ! -z "${rev_reads}" ]; then rev_dir=$(dirname "$rev_reads"); fi
  if [ ! -z "${seqSum}" ]; then seqSum_dir=$(dirname "$seqSum"); fi
  if [ ! -z "${centrif_DB}" ]; then centrif_DB_dir=$(dirname "$centrif_DB"); fi
# getting absolute paths
  if [ ! -z "${assembly_file}" ]; then assembly_path=$(cd "$assembly_dir" && pwd); fi
  if [ ! -z "${nano_reads}" ]; then nano_path=$(cd "$nano_dir" && pwd); fi
  if [ ! -z "${fwd_reads}" ]; then fwd_path=$(cd "$fwd_dir" && pwd); fi
  if [ ! -z "${rev_reads}" ]; then rev_path=$(cd "$rev_dir" && pwd); fi
  if [ ! -z "${seqSum}" ]; then  seqSum_path=$(cd "$seqSum_dir" && pwd); fi
  if [ ! -z "${input_folder}" ]; then infolder_path=$(cd "$input_folder" && pwd); fi
  if [ ! -z "${centrif_DB}" ]; then centrif_DB_path=$(cd "$centrif_DB_dir" && pwd); fi
# getting filename w/o path
  assembly_name=${assembly_file##*/}
  nano_name=${nano_reads##*/}
  fwd_file=${fwd_reads##*/}
  rev_file=${rev_reads##*/}
  centrif_DB_file=${centrif_DB##*/}
#############################
## Choose Executable(s)    ##
#############################
# [-C] Taxonomic read classification
  if [ ! -z "${tax_read}" ] && [ ! -z "${nano_reads}" ]; then centrifuge_nanopore ; fi
  if [ ! -z "${tax_read}" ] && [ ! -z "${fwd_reads}" ] && [ ! -z "${rev_reads}" ]; then centrifuge_illumina; fi
  # [-K] Krona summary
  if [ ! -z "${recentrifuge}" ] && [ ! -z "${input_folder}" ]; then recentrifuge ; fi
# Nanopore QC
  if [ ! -z "${nanoQC}" ] && [ ! -z "${seqSum}" ]; then QC_nanopore; fi
# Plasmid binning
  if [ ! -z "${plasflow}" ] && [ ! -z "${assembly_file}" ]; then plasflow_execute; fi
# Sourmash assembly classification
  if [ ! -z "${sour_cluster}" ] && [ ! -z "${assembly_file}" ]; then sourmash_cluster; fi
# [-D] Deepvir prediction
  if [ ! -z "${deepvir}" ] && [ ! -z "${assembly_file}" ]; then deepvirfinder_excecute; fi
# [-R] Resistance gene screening
  if [ ! -z "${AB_res}" ] && [ ! -z "${assembly_file}" ]; then resistance_screen; fi
  if [ ! -z "${AB_res}" ] && [ ! -z "${nano_reads}" ]; then resistance_read_screen; fi
# Binning
  if [ ! -z "${binning}" ] && [ ! -z "${assembly_file}" ] && [ ! -z "${fwd_reads}" ] && [ ! -z "${rev_reads}" ]; then binning_execute ; fi
