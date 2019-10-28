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
println "\033[2mCurrent User: $workflow.userName"
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
if (params.fasta == '' &&  params.fastq == '' &&  params.dir == '' && params.fastqPair == '' && !params.dev) {
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

if (params.fastqPair ) { fastqPair_input_ch = Channel
        .fromFilePairs( params.fastqPair , checkIfExists: true )
        .view() }

workflow {
/************************** 
* DATABASES
**************************/
// sourmash
    if (params.sour_db) { database_sourmash = file(params.sour_db) }

    else if (workflow.profile == 'gcloud' && (params.sourmeta || params.sourclass)) {
            sour_db_preload = file("gs://databases-nextflow/databases/sourmash/genbank-k31.lca.json")    
        if (sour_db_preload.exists()) { database_sourmash = sour_db_preload }
        else {  include './modules/sourmashgetdatabase'
                sourmash_download_db() 
                database_sourmash = sourmash_download_db.out } }

    else if (params.sourmeta || params.sourclass) {
                include './modules/sourmashgetdatabase'
                sourmash_download_db() 
                database_sourmash = sourmash_download_db.out }

// metamaps
    if (params.tax_db) {
    database_metamaps = file(params.tax_db, checkIfExists: true, type: 'dir') }

// gtdbtk
    if (params.gtdbtk_db) { database_gtdbtk = file(params.gtdbtk_db) }

    else if (workflow.profile == 'gcloud' && params.gtdbtk) {
            gtdbtk_preload = file("gs://databases-nextflow/databases/gtdbtk/gtdbtk_r89_data.tar.gz")    
        if (gtdbtk_preload.exists()) { database_gtdbtk = gtdbtk_preload }
        else {  include './modules/gtdbtkgetdatabase'
                gtdbtk_download_db() 
                database_gtdbtk = gtdbtk_download_db.out } }

    else if (params.gtdbtk) {
                include './modules/gtdbtkgetdatabase'
                gtdbtk_download_db() 
                database_gtdbtk = gtdbtk_download_db.out }

// centrifuge
    if (params.centrifuge_db) { database_centrifuge = file(params.centrifuge_db) }

    else if (workflow.profile == 'gcloud' && params.centrifuge) {
            centrifuge_preload = file("gs://databases-nextflow/databases/centrifuge/gtdb_r89_54k_centrifuge.tar")    
        if (centrifuge_preload.exists()) { database_centrifuge = centrifuge_preload }
        else {  include './modules/centrifugegetdatabase'
                centrifuge_download_db() 
                database_centrifuge = centrifuge_download_db.out } }

    else if (params.centrifuge) {
                include './modules/centrifugegetdatabase'
                centrifuge_download_db() 
                database_centrifuge = centrifuge_download_db.out } 

/************************** 
* MODULES
**************************/
/*************  
* --abricate | resistance screening
*************/
    if (params.abricate && params.fasta) { 
        include './modules/PARSER/abricateParserFASTA' params(output: params.output)
        include './modules/PLOTS/abricatePlotFASTA' params(output: params.output)
        include './modules/abricate' params(output: params.output)
        method = ['argannot', 'card', 'ncbi', 'plasmidfinder', 'resfinder']
        abricatePlotFASTA(abricateParserFASTA(abricate(fasta_input_ch, method))) }

    if (params.abricate && params.fastq) { 
        include './modules/abricateBatch'
        include './modules/PARSER/abricateParser' params(output: params.output)
        include './modules/PLOTS/abricatePlot' params(output: params.output)
        include './modules/fastqTofasta' params(output: params.output)
        method = ['ncbi', 'plasmidfinder']
        //split fastq batches
        abricateBatch(fastqTofasta(fastq_input_ch.splitFastq(by: 100000, file: true)), method) 
            abricateBatch.out.collectFile(storeDir: "${params.output}/abricate-batch", skip: 1, keepHeader: true)
            newChannel = abricateBatch.out.collectFile(skip: 1, keepHeader: true).map { file -> tuple(file.baseName, file)}
        abricatePlot(abricateParser(newChannel)) }

/*************  
* --centrifuge | tax classification of fastq
*************/
    if (params.centrifuge && params.fastq) { include './modules/centrifuge' params(output: params.output) 
        centrifuge(fastq_input_ch,database_centrifuge) }
/*************  
* --deepHumanPathogen | Deepsequencing of Human pathogens 
*************/
    if (params.deepHumanPathogen && params.fastqPair) { 
            include './modules/downloadHuman' params(output: params.output) 
            include './modules/bwaUnmapped' params(output: params.output) 
            include './modules/removeViaMapping' params(output: params.output) 
        removeViaMapping(bwaUnmapped(fastqPair_input_ch,downloadHuman())) }
/*************  
* --gtdbtk | tax classification of fastas
*************/
    if (params.gtdbtk && params.dir) { include './modules/gtdbtk' params(output: params.output) 
        gtdbtk(dir_input_ch,database_gtdbtk) }
/*************  
* --guppy-gpu | basecalling via guppy
*************/
    if (params.guppygpu && params.dir) { include './modules/basecalling' params(output: params.output, flowcell: params.flowcell, barcode: params.barcode, kit: params.kit ) 
        basecalling(dir_input_ch) }
/*************  
* --metamaps | taxonomic read classification
*************/
    if (params.metamaps && params.fastq && params.tax_db) { 
        include './modules/metamaps' params(output: params.output, memory: params.memory) 
        include './modules/krona' params(output: params.output)
        krona(metamaps(fastq_input_ch, database_metamaps)) }
/*************  
* --nanoplot | read quality via nanoplot
*************/
    if (params.nanoplot && params.fastq) { include './modules/nanoplot' params(output: params.output) 
        nanoplot(fastq_input_ch) }
/*************  
* --plasflow | plasmidprediction
*************/
    if (params.plasflow && params.fasta) { include './modules/plasflow' params(output: params.output) 
        plasflow(fasta_input_ch) }
/*************  
* --sourmeta | Metagenomic classification via sourmash
*************/
    if (params.sourmeta && params.fastq) {
        include './modules/sourmeta' params(output: params.output, fasta: params.fasta, fastq: params.fastq)
        include './modules/fastqTofasta' params(output: params.output)
        include './modules/rmetaplot' params(output: params.output)
        rmetaplot(sourmashmeta(fastqTofasta(fastq_input_ch),database_sourmash)) }

    if (params.sourmeta && params.fasta) { include './modules/sourmeta' params(output: params.output, fasta: params.fasta, fastq: params.fastq)
        sourmashmeta(fasta_input_ch,database_sourmash) }
/*************  
* --sourclass | Taxonomic classification via sourmash
*************/
    if (params.sourclass && params.fasta) { include './modules/sourclass' params(output: params.output) 
        sourmashclass(fasta_input_ch,database_sourmash) }
/*************  
* --sourcluster | Sequence clustering via sourmash 
*************/
    if (params.sourcluster && params.fasta) { include './modules/sourclusterfasta' params(output: params.output) 
        sourmashclusterfasta(fasta_input_ch) }
    else if (params.sourcluster && params.dir ) { include './modules/sourclusterdir' params(output: params.output) 
        sourmashclusterdir(dir_input_ch) }

/*************  
* --dev | Quick test and defs modules
*************/

if (params.dev ) { 
    include './modules/dev' params(output: params.output)
    repeater = ['8', '16', '24', '32', '40', '48']
    databasefile = file("gs://databases-nextflow/databases/thinspace/4centrifuge.tar.gz")
        dev(databasefile, repeater) }

}
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
    ${c_blue} --abricate ${c_reset}          antibiotic and plasmid screening    ${c_green}[--fasta]${c_reset} or ${c_green}[--fastq]${c_reset}
    ${c_blue} --centrifuge ${c_reset}        metagenomic classification of long reads  ${c_green} [--fastq]${c_reset}
    ${c_dim}  ..option flags:            [--centrifuge_db] path to your own DB instead, either .tar or .tar.gz ${c_reset}
    ${c_blue} --deepHumanPathogen ${c_reset} pathogen identification in human  ${c_green} [--fastqPair '*_R{1,2}.fastq.gz']${c_reset}
    ${c_blue} --gtdbtk ${c_reset}            tax. class. via marker genes        ${c_green}[--dir]${c_reset}
    ${c_dim}  ..option flags:            [--gtdbtk_db] path to your own DB instead ${c_reset}
    ${c_blue} --sourmeta ${c_reset}          metagenomic analysis "WIMP"         ${c_green}[--fasta]${c_reset} or ${c_green}[--fastq]${c_reset}
    ${c_dim}  ..option flags:            [--sour_db] path to your own DB instead ${c_reset}
    ${c_blue} --sourclass ${c_reset}         taxonomic classification            ${c_green}[--fasta]  ${c_reset}
    ${c_dim}  ..option flags:            [--sour_db] path to your own DB instead ${c_reset}
    ${c_blue} --sourcluster ${c_reset}       sequence comparision with kmers     ${c_green}[--fasta]${c_reset} or ${c_green}[--dir]${c_reset}
    ${c_dim}  ..inputs:                  multi-fasta: --fasta; multiple files: --dir${c_reset}
    ${c_blue} --nanoplot  ${c_reset}         read quality via nanoplot           ${c_green}[--fastq]${c_reset}
    ${c_blue} --guppygpu ${c_reset}          basecalling via guppy-gpu-nvidia   ${c_green} [--dir]${c_reset}
    ${c_dim}  ..option flags:            [--flowcell] [--kit] [--barcode]
      ..default settings:        [--flowcell $params.flowcell] [--kit $params.kit] ${c_reset}
    ${c_blue} --plasflow ${c_reset}          predicts & seperates plasmid-seqs${c_green}   [--fasta]${c_reset}
    ${c_blue} --metamaps ${c_reset}          metagenomic classification of long reads  ${c_green} [--fastq]${c_reset}
    ${c_dim}  ..mandatory flags:         [--memory] [--tax_db] e.g. --memory 100 --tax_db /databases/miniSeq+H 

    ${c_yellow}Options:${c_reset}
    --cores             max cores for local use [default: $params.cores]
    --memory            80% of available RAM in GB for --metamaps [default: $params.memory]
    --output            name of the result folder [default: $params.output]

    ${c_yellow}Database(s) behaviour${c_reset}
    The Priority is:
    1. Use your own DB via the flags, e.g. [--sour_db] [--gtdbtk_db]
    2. Without a flag it downloads/retrives a database to/from: ./nextflow-autodownload-databases
    ${c_dim}If a auto download is not possible or implemented the workflow will tell you

    ${c_dim}Nextflow options:
    -with-report rep.html    cpu / ram usage (may cause errors)
    -with-dag chart.html     generates a flowchart for the process tree
    -with-timeline time.html timeline (may cause errors)

    Profile:
    -profile                 standard, gcloud, UKJ [default: standard] ${c_reset}
    """.stripIndent()
}