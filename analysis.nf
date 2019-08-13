nextflow.preview.dsl=2

/*
* Nextflow -- Analysis Pipeline
* Author: christian.jena@gmail.com
*/

if (params.help) { exit 0, helpMSG() }
if (params.fasta == '' &&  params.fastq == '' &&  params.dir == '') {
    exit 1, "input missing, use [--fasta] [--fastq] or [--dir]"}

// fasta input or via csv file
if (params.fasta && params.list) { fasta_input_ch = Channel
        .fromPath( params.fasta, checkIfExists: true )
        .splitCsv()
        .map { row -> ["${row[0]}", file("${row[1]}")] }
        .view() }
else if (params.fasta) { fasta_input_ch = Channel
        .fromPath( params.fasta , checkIfExists: true)
        .map { file -> tuple(file.baseName, file) }
        .view() }

if (params.fastq && params.list) { fastq_input_ch = Channel
        .fromPath( params.fastq, checkIfExists: true )
        .splitCsv()
        .map { row -> ["${row[0]}", file("${row[1]}")] }
        .view() }
else if (params.fastq) { fastq_input_ch = Channel
        .fromPath( params.fastq , checkIfExists: true)
        .map { file -> tuple(file.baseName, file) }
        .view() }

if (params.dir) { dir_input_ch = Channel
        .fromPath( params.dir , checkIfExists: true)
        .map { file -> tuple(file.name, file) }
        .view() }


/************* 
* DATABASES
*************/
if (params.sourmeta || params.sourclass) {
    sour_db_preload = file(params.sour_db_present)
    if (params.sour_db) { database_sourmash = file(params.sour_db) }
    else if (sour_db_preload.exists()) { database_sourmash = sour_db_preload }
    else {  include 'modules/sourmashgetdatabase'
            sourmash_download_db() 
            database_sourmash = sourmash_download_db.out } }

/*************  
* --sourmeta | Metagenomic classification via sourmash
*************/
if (params.sourmeta && params.fasta) { include 'modules/sourmeta' params(output: params.output, cpus: params.cpus)
    sourmashmeta(fasta_input_ch,database_sourmash) }

/*************  
* --sourclass | Taxonomic classification via sourmash
*************/
if (params.sourclass && params.fasta) { include 'modules/sourclass' params(output: params.output, cpus: params.cpus) 
    sourmashclass(fasta_input_ch,database_sourmash) }

/*************  
* --nanoplot | read quality via nanoplot
*************/
if (params.nanoplot && params.fastq) { include 'modules/nanoplot' params(output: params.output, cpus: params.cpus) 
    nanoplot(fastq_input_ch) }

/*************  
* --guppy-gpu | basecalling via guppy
*************/
if (params.guppygpu && params.dir) { include 'modules/basecalling' params(output: params.output, flowcell: params.flowcell, barcode: params.barcode, kit: params.kit ) 
    basecalling(dir_input_ch) }

// ***  TODO:

// /*************  
// * --sourcluster | Sequence clustering via sourmash   TODO
// *************/
// if (params.sourclass) { include 'modules/sourclass' params(output: params.output, cpus: params.cpus) 
//     sourmashclass(fasta_input_ch,database_sourmash) }

// /*************  
// * --plasflow | plasmidprediction   TODO
// *************/
// if (params.plasflow) { include 'modules/sourclass' params(output: params.output, cpus: params.cpus) 
//     sourmashclass(fasta_input_ch,database_sourmash) }

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
    log.info """

    Usage example:
    nextflow run replikation/docker_pipelines --fasta '*/*.fasta' --sourmeta --sourclass

    Input:
    --fasta             '*.fasta'   -> assembly file(s) 
    --fastq             '*.fastq'   -> read file(s) in fastq, one sample per file
    --dir               'foobar*/'  -> a folder(s) as input
    --list              activates csv input for --fasta --fastq instead of fasta/q files
 

    Workflows [Input needed]:
    --sourmeta          metagenomic sourmash analysis       [--fasta]
    --sourclass         taxonomic sourmash classification   [--fasta]  
    --nanoplot          read quality via nanoplot           [--fastq]
    --guppygpu          basecalling via guppy-gpu-nvidia    [--dir]
    .. option flags:            [--flowcell] [--kit] [--barcode]
    .. default settings:        [--flowcell $params.flowcell] [--kit $params.kit]


   Options:
    --cores             max cores for local use [default: $params.cores]
    --output            name of the result folder [default: $params.output]


    Nextflow options:
    -with-report rep.html    cpu / ram usage (may cause errors)
    -with-dag chart.html     generates a flowchart for the process tree
    -with-timeline time.html timeline (may cause errors)

    Profile:
    -profile                 standard [default: standard]
    """.stripIndent()
}