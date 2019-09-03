#!/usr/bin/env nextflow
nextflow.preview.dsl=2

/*
* Nextflow -- Analysis Pipeline
* Author: christian.jena@gmail.com
*/

// terminal prints for pipeline overview
println " "
println "\u001B[32mProfile: $workflow.profile\033[0m"
println " "
println "\033[2mPipeline-version: $workflow.commitId"
println "Nextflow-version: $nextflow.version"
println "Starting time: $nextflow.timestamp"
println "Workdir location:"
println "  $workflow.workDir\u001B[0m"
println " "
if (workflow.profile == 'standard') {
println "\033[2mCPUs to use: $params.cores"
println "Output dir name: $params.output\u001B[0m"
println " "}

// helpmessages & user input checks
if (params.help) { exit 0, helpMSG() }
if (params.profile) {
    exit 1, "--profile is WRONG use -profile" }
if (params.fasta == '' &&  params.fastq == '' &&  params.dir == '') {
    exit 1, "input missing, use [--fasta] [--fastq] or [--dir]"}
if (params.fasta && params.fastq) {
    exit 1, "please us either: [--fasta] or [--fastq]"}   
if (params.fastq && params.metamaps && params.tax_db == '') {
    exit 1, "taxonomic database location not specified via [--tax_db]"}

// fasta input or via csv file
if (params.fasta && params.list) { fasta_input_ch = Channel
        .fromPath( params.fasta, checkIfExists: true )
        .splitCsv()
        .map { row -> ["${row[0]}", file("${row[1]}")] }
        .view() }
else if (params.fasta) { fasta_input_ch = Channel
        .fromPath( params.fasta, checkIfExists: true)
        .map { file -> tuple(file.baseName, file) }
        .view() }

if (params.fastq && params.list) { fastq_input_ch = Channel
        .fromPath( params.fastq, checkIfExists: true )
        .splitCsv()
        .map { row -> ["${row[0]}", file("${row[1]}")] }
        .view() }
else if (params.fastq) { fastq_input_ch = Channel
        .fromPath( params.fastq, checkIfExists: true)
        .map { file -> tuple(file.baseName, file) }
        .view() }

if (params.dir && params.list) { dir_input_ch = Channel
        .fromPath( params.dir, checkIfExists: true, type: 'dir' )
        .splitCsv()
        .map { row -> ["${row[0]}", file("${row[1]}")] }
        .view() }
if (params.dir) { dir_input_ch = Channel
        .fromPath( params.dir, checkIfExists: true, type: 'dir')
        .map { file -> tuple(file.name, file) }
        .view() }


/************************** 
* DATABASES
**************************/
// sourmash
if (params.sourmeta || params.sourclass) {
    if (params.sour_db) { database_sourmash = file(params.sour_db) }
    else {  include 'modules/sourmashgetdatabase'
            sourmash_download_db() 
            database_sourmash = sourmash_download_db.out } }
// metamaps
if (params.tax_db) {
database_metamaps = file(params.tax_db, checkIfExists: true, type: 'dir') }

// gtdbtk
if (params.gtdbtk) {
    if (params.gtdbtk_db) { database_gtdbtk = file(params.gtdbtk_db) }
    else {  include 'modules/gtdbtkgetdatabase'
            gtdbtk_download_db() 
            database_gtdbtk = gtdbtk_download_db.out } }


/************************** 
* MODULES
**************************/
/*************  
* --gtdbtk | tax classification of fastas
*************/
if (params.gtdbtk && params.dir) { include 'modules/gtdbtk' params(output: params.output) 
    gtdbtk(dir_input_ch,database_gtdbtk) }

/*************  
* --guppy-gpu | basecalling via guppy
*************/
if (params.guppygpu && params.dir) { include 'modules/basecalling' params(output: params.output, flowcell: params.flowcell, barcode: params.barcode, kit: params.kit ) 
    basecalling(dir_input_ch) }

/*************  
* --metamaps | taxonomic read classification
*************/
if (params.metamaps && params.fastq && params.tax_db) { 
    include 'modules/metamaps' params(output: params.output, memory: params.memory) 
    include 'modules/krona' params(output: params.output)
    krona(metamaps(fastq_input_ch, database_metamaps)) }

/*************  
* --nanoplot | read quality via nanoplot
*************/
if (params.nanoplot && params.fastq) { include 'modules/nanoplot' params(output: params.output) 
    nanoplot(fastq_input_ch) }

/*************  
* --plasflow | plasmidprediction
*************/
if (params.plasflow && params.fasta) { include 'modules/plasflow' params(output: params.output) 
    plasflow(fasta_input_ch) }

/*************  
* --sourmeta | Metagenomic classification via sourmash
*************/
if (params.sourmeta && params.fastq) {
    include 'modules/sourmeta' params(output: params.output, fasta: params.fasta, fastq: params.fastq)
    include 'modules/fastqTofasta' params(output: params.output)
    include 'modules/rmetaplot' params(output: params.output)
    rmetaplot(sourmashmeta(fastqTofasta(fastq_input_ch),database_sourmash)) }

if (params.sourmeta && params.fasta) { include 'modules/sourmeta' params(output: params.output, fasta: params.fasta, fastq: params.fastq)
    sourmashmeta(fasta_input_ch,database_sourmash) }

/*************  
* --sourclass | Taxonomic classification via sourmash
*************/
if (params.sourclass && params.fasta) { include 'modules/sourclass' params(output: params.output) 
    sourmashclass(fasta_input_ch,database_sourmash) }

/*************  
* --sourcluster | Sequence clustering via sourmash 
*************/
if (params.sourcluster && params.fasta) { include 'modules/sourclusterfasta' params(output: params.output) 
    sourmashclusterfasta(fasta_input_ch) }
else if (params.sourcluster && params.dir ) { include 'modules/sourclusterdir' params(output: params.output) 
    sourmashclusterdir(dir_input_ch) }

// /*************  
// * --abricate | resistance screening   TODO
// *************/
// if (params.abricate) { include 'modules/sourclass' params(output: params.output, cpus: params.cpus) 
//     sourmashclass(fasta_input_ch,database_sourmash) }

// /*************  
// * --abricate-reads | resistance screening  of reads TODO
// *************/
// // add a fastq input, needs matching file name to work (.join)
// if (params.abricate) { include 'modules/sourclass' params(output: params.output, cpus: params.cpus) 
//     sourmashclass(fasta_input_ch,database_sourmash) }



/*************  
* --help
*************/
def helpMSG() {
    c_green = "\033[0;32m";
    c_reset = "\033[0m";
    c_yellow = "\033[0;33m";
    c_blue = "\033[0;34m";
    c_dim = "\033[2m";
    log.info """
    ____________________________________________________________________________________________
    
    Nextflow Analysis modules for easy use, by Christian Brandt
    
    ${c_yellow}Usage example:${c_reset}
    nextflow run replikation/docker_pipelines --fasta '*/*.fasta' --sourmeta --sourclass

    ${c_yellow}Input:${c_reset}
    ${c_green} --fasta ${c_reset}            '*.fasta'   -> assembly file(s) - uses filename
    ${c_green} --fastq ${c_reset}            '*.fastq'   -> read file(s) in fastq, one sample per file - uses filename
    ${c_green} --dir  ${c_reset}             'foobar*/'  -> folder(s) as input - uses dirname
    ${c_dim}  ..change above input to csv:${c_reset} ${c_green}--list ${c_reset}            

    ${c_yellow}Workflows:${c_reset}
    ${c_blue} --gtdbtk ${c_reset}            tax. class. via marker genes        ${c_green}[--dir]${c_reset}
    ${c_dim}  ..option flags:            [--gtdbtk_db] path to your own DB instead ${c_reset}
    ${c_blue} --sourmeta ${c_reset}          metagenomic sourmash analysis       ${c_green}[--fasta]${c_reset}
    ${c_dim}  ..option flags:            [--sour_db] path to your own DB instead ${c_reset}
    ${c_blue} --sourclass ${c_reset}         taxonomic sourmash classification   ${c_green}[--fasta]  ${c_reset}
    ${c_dim}  ..option flags:            [--sour_db] path to your own DB instead ${c_reset}
    ${c_blue} --sourcluster ${c_reset}       sequence comparision with kmers     ${c_green}[--fasta]${c_reset} or ${c_green}[--dir]${c_reset}
    ${c_dim}  ..inputs:                  multi-fasta: --fasta; multiple files: --dir${c_reset}
    ${c_blue} --nanoplot  ${c_reset}         read quality via nanoplot           ${c_green}[--fastq]${c_reset}
    ${c_blue} --guppygpu ${c_reset}          basecalling via guppy-gpu-nvidia   ${c_green} [--dir]${c_reset}
    ${c_dim}  ..option flags:            [--flowcell] [--kit] [--barcode]
      ..default settings:        [--flowcell $params.flowcell] [--kit $params.kit] ${c_reset}
    ${c_blue} --plasflow ${c_reset}          predicts & seperates plasmid-seqs${c_green}   [--fasta]${c_reset}
    ${c_blue} --metamaps ${c_reset}          metagenomic classification of long reads  ${c_green} [--fastq]${c_reset}
    ${c_dim}  ..mandatory:               [--memory] [--tax_db] e.g. --memory 100 --tax_db /databases/miniSeq+H 

    ${c_yellow}Options:${c_reset}
    --cores             max cores for local use [default: $params.cores]
    --memory            80% of available RAM in GB for --metamaps [default: $params.memory]
    --output            name of the result folder [default: $params.output]

    ${c_yellow}Database(s) behaviour${c_reset}
    The Priority is:
    1. Use your own DB via the flags, e.g. [--sour_db] [--gtdbtk_db]
    2. Without a flag the pre downloaded database will be used, located in ./db_auto-build
    3. If 1. and 2. are not available it will download the database and put it into ./db_auto-build
    ${c_dim}If a auto download is not possible or implemented the workflow will tell you

    ${c_dim}Nextflow options:
    -with-report rep.html    cpu / ram usage (may cause errors)
    -with-dag chart.html     generates a flowchart for the process tree
    -with-timeline time.html timeline (may cause errors)

    Profile:
    -profile                 standard, gcloud [default: standard] ${c_reset}
    """.stripIndent()
}