#!/usr/bin/env bash

## Docker ##
  type docker >/dev/null 2>&1 || { echo -e >&2 "${RED}Docker not found. Please run the installation skript, Aborting.${NC}"; exit 1; }

## OPTIONS ##
  WORKDIRPATH=$(pwd) # for docker mountpoint (-v)
  WORKDIRNAME=${PWD##*/} # for docker mountpoint (-v)

  # FASTQ prefixes for forward and reverse pairs
    FPREFIX="_R1_" #write down current prefix to use
    RPREFIX="_R2_"

  # Foldernames, you can change your folder names in the quotes
    # Server locations (for option a)
    FASTQ_serv="/volume1/Shared_Regensburg/FASTQ_Rohdaten"        # RAWDATA on Server
    ASSEMBLY_serv="/volume1/Shared_Regensburg/ASSEMBLIES_fertig"  # ASSEMBLY on Server
    CLUSTER_serv="/volume1/Shared_Regensburg/ANALYSEN"            # CLUSTER Results on Server
    IP=$(cat ~/cfg)                                               # location of USERNAME@IPADRESS
    # Foldernames in Current Workdir (for all options)
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

###############
### Modules ###
###############
serv_check_downl()
  {
    echo "Download & Assembly files from $FASTQ_serv ?"
    read -p "Proceed? [yes] or [no] " yn
    case $yn in
        [Yy]* ) echo "Starting Pipeline";;
        [Nn]* ) echo "Exiting..."; exit;;
        * ) echo "  Please answer [y] or [n].";;
    esac
    RAW_files=$(ssh $IP ls ${FASTQ_serv})
    while IFS= read -r samplename ; do
       echo "$samplename"
       if ssh < /dev/null $IP stat ${ASSEMBLY_serv}/${samplename} \> /dev/null 2\>\&1; then
            echo "${samplename} exists"
          else
            echo "${samplename} does not exist, starting download"
            mkdir -p $FASTQ_Wdir
            scp -r ${IP}:${FASTQ_serv}/${samplename} ${FASTQ_Wdir}
       fi
    done < <(printf '%s\n' "$RAW_files")
  }

## QC ##
read_quality()
  {
    mkdir -p $ASSEMBLY
    for illuminarun in ${FASTQ_Wdir}/* ; do
      F_file=$(ls ${illuminarun}/*${FPREFIX}*)
      R_file=$(ls ${illuminarun}/*${RPREFIX}*)
      SampleID=${illuminarun##*/}
      mkdir -p ${ASSEMBLY}/${SampleID}
      # fastqc run
      docker run --rm -it -v ${WORKDIRPATH}:/${WORKDIRNAME} replikation/fastqc \
      ${WORKDIRNAME}/${F_file} /${WORKDIRNAME}/${R_file} -t $CPU --outdir=/${WORKDIRNAME}/${ASSEMBLY}/${SampleID}
    done
  }

## Assembler ##
unicycler_execute()
  {
    for illuminarun in ${FASTQ_Wdir}/* ; do
      F_file=$(ls ${illuminarun}/*${FPREFIX}*)
      R_file=$(ls ${illuminarun}/*${RPREFIX}*)
      SampleID=${illuminarun##*/}
      mkdir -p ${ASSEMBLY}/${SampleID}
      # fastqc run
      docker run --rm -it -v ${WORKDIRPATH}:/${WORKDIRNAME} replikation/unicycler \
      -1 /${WORKDIRNAME}/$F_file -2 /${WORKDIRNAME}/$R_file -t $CPU -o /${WORKDIRNAME}/${ASSEMBLY}/${SampleID}
      # rename assembly file
      mv ${illuminarun}/assembly.fasta ${illuminarun}/${SampleID}.fasta
    done
    # remove all gfa files
    rm ${ASSEMBLY}/*/*.gfa
  }

## sourmash ##
serv_get_fasta()
  {
    echo "Getting all Assembly files"
    mkdir -p ${ASSEMBLY}
    scp -r ${IP}:${ASSEMBLY_serv}/* ${ASSEMBLY}/
  }


sourmash_execute()
  {
  #  mkdir -p $CLUSTER/signatures
  #  for illuminarun in ${ASSEMBLY}/*/*.fasta ; do
  #    SampleID=${illuminarun##*/}
  #    docker run --rm -it -v ${WORKDIRPATH}:/${WORKDIRNAME} replikation/sourmash \
  #    sourmash compute -n 5000 -k 31,51 /${WORKDIRNAME}/${illuminarun} -o /${WORKDIRNAME}/${CLUSTER}/signatures/${SampleID%.fasta}.sig
    #done
    # index
    #docker run --rm -it -v ${WORKDIRPATH}:/${WORKDIRNAME} replikation/sourmash \
    #/bin/sh -c "sourmash compare /${WORKDIRNAME}/${CLUSTER}/signatures/* -o /${WORKDIRNAME}/${CLUSTER}/signatures/index_sig"
    # plot
    #docker run --rm -t -v ${WORKDIRPATH}:/${WORKDIRNAME}  replikation/sourmash \
    #/bin/sh -c "sourmash plot --pdf --labels /${WORKDIRNAME}/${CLUSTER}/signatures/index_sig ; mv *.pdf /${WORKDIRNAME}/${CLUSTER} ; mv *.png /${WORKDIRNAME}/${CLUSTER}"

####OPTION 2
### trying to get better lables out of it without paths and .fasta
  mkdir -p $CLUSTER/signatures
  docker run --rm -it -v ${WORKDIRPATH}/${ASSEMBLY}:/${ASSEMBLY} replikation/sourmash \
  /bin/sh -c 'mkdir /signatures_tmp && cd /signatures_tmp && \
              cp /${ASSEMBLY}/*/*.fasta . && \
              for fasta in *.fasta; do mv ${fasta} ${fasta%.fasta}; done && \
              for fasta in *; do sourmash compute -n 5000 -k 51 ${fasta} -o ${fasta}.sig; done && \
              sourmash compare *.sig -o index_sig && \
              sourmash plot --pdf --labels index_sig && mv *.pdf /${WORKDIRNAME}/${CLUSTER}'
  }

serv_upload()
  {
    echo "uploading to server"
    scp -r $ASSEMBLY/*/ $IP:$ASSEMBLY_serv #check how the file structure looks after upload
    scp -r ${CLUSTER}/* $IP:$CLUSTER_serv #check how the file structure looks after upload
  }

############################
###   Start of script    ###
############################
echo "                                               _____________________________"
echo "______________________________________________/ Created by Christian Brandt \___"
echo " "
while true; do
    echo "Illumina Paired end, automated Pipeline"
    echo -e "The following fastq IDs for forward/reverse are used: ${YEL}${FPREFIX}${NC} and ${YEL}${RPREFIX}${NC}"
    echo " "
    echo -e "${GRE}What do you want to do? [a] [c] [s] [e]xit${NC}"
    echo -e "[a] ${YEL}Automated server analysis - assembly, clusteranalysis${NC}"
    echo "    Start in a empty folder, downloads and uploades files automatically"
    echo -e "[c] ${YEL}Current Workingdir Assembly${NC} - unicycler assembly"
    echo "    Expects one folder for each illumina run in ./$FASTQ_Wdir/, give folders a meaningful name or ID"
    echo -e "[s] ${YEL}Current Workingdir Clusteranalysis${NC} - sourmash cluster analysis"
    echo "    Expects one folder for each assembly in ./$ASSEMBLY/, give fasta files a meaningful name or ID"
    echo " "
    read -p "Automated[a] Current Wdir Assembly[c] Sourmash cluster [s] or [e]xit: " acsde
    case $acsde in
        [Aa]* ) serv_check_downl; read_quality; unicycler_execute; serv_get_fasta; sourmash_execute; serv_upload; break;;
        [Cc]* ) read_quality; unicycler_execute; break;;
        [Ss]* ) sourmash_execute; break;;
        [Dd]* ) sourmash_execute; break;;
        [Ee]* ) echo "Exiting script, bye bye"; exit;;
        * ) echo "  Please answer [a] [c] [s] or [e].";;
    esac
done
echo -e "${GRE} Finished all steps ${NC}"
