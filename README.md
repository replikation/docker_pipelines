# Nextflow - Docker pipeline

![](https://img.shields.io/badge/nextflow-20.07.1-brightgreen)
![](https://img.shields.io/badge/uses-docker-blue.svg)
![](https://img.shields.io/badge/licence-GPL--3.0-lightgrey.svg)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4030719.svg)](https://doi.org/10.5281/zenodo.4030719)

![](https://github.com/replikation/docker_pipelines/workflows/Syntax_check/badge.svg)

## Overview

* This is a pipeline overview containing various tools in an easy install free nextflow environment
* you dont need to install any tools with the exeption of docker and nextflow

### brief tool overview

* abricate -> resistance detection (via fasta or fastq files)
* centrifuge -> taxonomic classification of nanopore reads (fastq)
* gtdbtk -> evaluation of bins (dir input)
* guppygpu -> gpu basecalling via guppy of fast5 files (dir input)
* metamaps -> taxonomic classification of nanopore reads (fastq) 
* nanoplot -> QC of nanopore reads (fastq)
* plasflow -> Plasmid prediction of assembly (fasta)
* sourclass -> Determines the species of a assembly via sourmash (fasta)
* sourcluster -> kmer based tree/clustering of contigs (fasta) or multiple assmeblies (via dir/)
* sourmeta -> metagenomic WIMP determination via sourmash
* deepHumanPathogen -> maps illumina reads against human to remove "human reads" (fastqPairs)
* (and more)

## Installation

* intended for linux systems

### Docker

* install `docker` via `sudo apt install docker.io`
	* or `docker-ce` via [this](https://docs.docker.com/install/linux/docker-ce/ubuntu/)
* add the docker group to your user via `sudo usermod -a -G docker $USER`
* reboot after this

### Nextflow

* install java runtime on your system
* install nextflow via ``curl -s https://get.nextflow.io | bash``

## Usage

````bash
nextflow run replikation/docker_pipelines --help
````

* update docker pipelines via 

````bash
nextflow pull replikation/docker_pipelines
````


### Citation

> Christian Brandt. (2020, September 15). Docker pipelines: a collection of nextflow workflows related to nanopore data, taxonomy and antibiotic resistance (Version v.1.0.0). Zenodo. http://doi.org/10.5281/zenodo.4030719

