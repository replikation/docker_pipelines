#!/usr/bin/env bash

## Docker ##
  type docker >/dev/null 2>&1 || { echo -e >&2 "${RED}Docker not found. Please run the installation skript, Aborting.${NC}"; exit 1; }
  echo "Docker identified"

## OPTIONS ##
  WORKDIRPATH=$(pwd) # for docker mountpoint (-v)
  WORKDIRNAME=${PWD##*/} # for docker mountpoint (-v)

  # Foldernames, you can change your folder names in the quotes
    # Server
    FASTQ_serv="<PATH_TO_SERVER>"         # RAWDATA on Server
    ASSEMBLY_serv="<PATH_TO_SERVER>"      # ASSEMBLY on Server
    CLUSTER_serv="<PATH_TO_SERVER>"       # CLUSTER Results on Server
    IP="cat IP"                           # location of IP address
    # Folders in Workingdir for pipeline
    FASTQ_Wdir="FASTQ"                  # FASTQ location (serv & Wdir)
    ASSEMBLY="ASSEMBLY"                 # tmp unicylcer output in WS
    CLUSTER="CLUSTER"                   # cluster analysis

  # CPU cores
    CPU=$(lscpu -p | egrep -v '^#' | wc -l) # can be changed to e.g. CPU="16"
  # colours, needed for echos
    RED='\033[0;31m'
    GRE='\033[0;32m'
    YEL='\033[0;33m'
    NC='\033[0m'

## PARAMETERS ##
# CONTAINS PARAMETERS FOR PROGRAMS, change at your own risk
  # sourmash

  # unicycler

  # PREFIX FASTQ
  PREFIX="_R1_" #write down current prefix to use

###############
### Modules ###
###############

# Workflow
# 1 search either server or current workdir for folders wir R1 and R2
    #a. compare $FASTQ_serv to $ASSEMBLY_serv, download missing to $FASTQ_Wdir
    #b. check what reads it can find in $FASTQ_Wdir/ (non server)

# 2. fastqc each and put html results into $ASSEMBLY
    #a. remove unused files

# 3. unicycler each folder with name and put assembly into $ASSEMBLY

# 4. clusteranalysis all files in $ASSEMBLY and put results to $CLUSTER
  # a.) cluster analysis should be a two option again.
  # either use all in $ASSEMBLY or feed a list and get it from the server


serv_check_assemblies()
  {
    echo "download all raw data without a assembly from $FASTQ_serv to $FASTQ_Wdir"
  }



## fastqc ##
read_quality()
  {
  echo "fastqc docker, put it into the assembly folder"
  }


unicycler_execute()
  {
    echo "unicycler"
  }

## Assembler ##
direct_read_use()
  {
  echo "filters now sequences greater than 2000 bp into a no-assembly-fasta-file for analysis"
  for fastqfile in ${FASTQ}/*.fastq; do sed -n '1~4s/^@/>/p;2~4p' $fastqfile > ${fastqfile%.fastq}.fasta ; done
  for fastafile in ${FASTQ}/*.fasta; do sed ':a;N;/^>/M!s/\n//;ta;P;D' $fastafile > ${fastafile%.fasta}_oneliner.fasta ; done
  for fastafile in ${FASTQ}/*_oneliner.fasta; do awk '/^>/ { getline seq } length(seq) >2000 { print $0 "\n" seq }' $fastafile > ${fastafile%_oneliner.fasta}_no_assembly_reads-only.fa ; done
  mv ${FASTQ}/*.fa ${ASSEMBLY}/
  rm ${FASTQ}/*.fasta
  echo "No-assembly-fasta-file stored under ${ASSEMBLY}"
  }

wtdbg2_execute()
  {
  echo -e "${RED} Running wtdbg2 with -e ${wtdbg2_edge_depth}; -S ${wtdbg2_subsampling}${NC}"
  for fastqfile in ${FASTQ}/*.fastq; do
    filename=$(basename $fastqfile)
    # L1
    docker run --rm -it -v ${WORKDIRPATH}:/${WORKDIRNAME} replikation/wtdbg2 \
    wtdbg2 -t $CPU -x ont -e ${wtdbg2_edge_depth} -S ${wtdbg2_subsampling} --rescue-low-cov-edges -i /${WORKDIRNAME}/${fastqfile} -o /${WORKDIRNAME}/${FASTA_raw}/${filename%.fastq}_${wtdbg2_L1}  -L ${wtdbg2_L1}
    # L2
    docker run --rm -it -v ${WORKDIRPATH}:/${WORKDIRNAME} replikation/wtdbg2 \
    wtdbg2 -t $CPU -x ont -e ${wtdbg2_edge_depth} -S ${wtdbg2_subsampling} --rescue-low-cov-edges -i /${WORKDIRNAME}/${fastqfile} -o /${WORKDIRNAME}/${FASTA_raw}/${filename%.fastq}_${wtdbg2_L2}  -L ${wtdbg2_L2}
    # L3
    docker run --rm -it -v ${WORKDIRPATH}:/${WORKDIRNAME} replikation/wtdbg2 \
    wtdbg2 -t $CPU -x ont -e ${wtdbg2_edge_depth} -S ${wtdbg2_subsampling} --rescue-low-cov-edges -i /${WORKDIRNAME}/${fastqfile} -o /${WORKDIRNAME}/${FASTA_raw}/${filename%.fastq}_${wtdbg2_L3}  -L ${wtdbg2_L3}
    # create contigs for all 3 runs
    docker run --rm -it -v ${WORKDIRPATH}:/${WORKDIRNAME} replikation/wtdbg2 \
    wtpoa-cns -t $CPU -i /${WORKDIRNAME}/${FASTA_raw}/${filename%.fastq}_${wtdbg2_L1}.ctg.lay.gz -fo /${WORKDIRNAME}/${ASSEMBLY}/${filename%.fastq}_${wtdbg2_L1}.fa
    docker run --rm -it -v ${WORKDIRPATH}:/${WORKDIRNAME} replikation/wtdbg2 \
    wtpoa-cns -t $CPU -i /${WORKDIRNAME}/${FASTA_raw}/${filename%.fastq}_${wtdbg2_L2}.ctg.lay.gz -fo /${WORKDIRNAME}/${ASSEMBLY}/${filename%.fastq}_${wtdbg2_L2}.fa
    docker run --rm -it -v ${WORKDIRPATH}:/${WORKDIRNAME} replikation/wtdbg2 \
    wtpoa-cns -t $CPU -i /${WORKDIRNAME}/${FASTA_raw}/${filename%.fastq}_${wtdbg2_L3}.ctg.lay.gz -fo /${WORKDIRNAME}/${ASSEMBLY}/${filename%.fastq}_${wtdbg2_L3}.fa
  done
  }

## Metagenome classification ##
centrifuge_execute()
  {
  for fastqfile in ${FASTQ}/*.fastq; do
    seqID=$(basename ${fastqfile%.fastq})
    docker run --rm -i -v ${WORKDIRPATH}:/${WORKDIRNAME} replikation/centrifuge \
    centrifuge -p $CPU -x /centrifuge/database/p_compressed -U /${WORKDIRNAME}/${fastqfile} -S /${WORKDIRNAME}/${TAX}/${seqID}_centrifuge_out.txt --report-file /${WORKDIRNAME}/${TAX}/${seqID}_centrifuge_out.log
    # create report for pavian
    docker run --rm -i -v ${WORKDIRPATH}:/${WORKDIRNAME} replikation/centrifuge \
    centrifuge-kreport -x /centrifuge/database/p_compressed /${WORKDIRNAME}/${TAX}/${seqID}_centrifuge_out.txt > ${WORKDIRPATH}/${TAX}/${seqID}_pavian_report.csv
  done
  }

plasflow_execute()
  {
  for fastafile in ${ASSEMBLY}/*.fa; do
    seqID=$(basename ${fastafile%.fa})
    docker run --rm -i -v ${WORKDIRPATH}:/${WORKDIRNAME} replikation/plasflow \
    filter_sequences_by_length.pl -input /${WORKDIRNAME}/${fastafile} -output /${WORKDIRNAME}/${TAX}/${seqID}.fasta -thresh 2000
    docker run --rm -i -v ${WORKDIRPATH}:/${WORKDIRNAME} replikation/plasflow \
    PlasFlow.py --input /${WORKDIRNAME}/${TAX}/${seqID}.fasta --output /${WORKDIRNAME}/${TAX}/${seqID}_plasflow_predictions.tsv --threshold 0.7
  done
  }

############################
###   Start of script    ###
############################
echo "                                               _____________________________"
echo "______________________________________________/ Created by Christian Brandt \___"
echo " "
while true; do
    echo -e "${GRE}What do you want to do? [a] [c] [e]xit${NC}"
    echo -e "[a] ${YEL}Automated server analysis${NC}"
    echo -e "[c] ${YEL}Current Workingdir analysis${NC}"
    read -p "Automated[a] Current Workdir[c] or [e]xit: " acd
    case $acd in
        [Aa]* ) basecaller_ask; porechop_ask; albacore_execute; porechop_execute; break;;
        [Cc]* ) direct_read_use; wtdbg2_execute; break;;
        [Dd]* ) centrifuge_execute; break;;
        [Ee]* ) echo "Exiting script, bye bye"; exit;;
        * ) echo "  Please answer [b] [m] [t] [p] or [e].";;
    esac
done
echo -e "${GRE} Finished all steps ${NC}"
