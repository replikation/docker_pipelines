nextflow.preview.dsl=2

/*
* Nextflow -- Analysis Pipeline
* Author: christian@jena@gmail.com
*/

if (params.help) { exit 0, helpMSG() }
if (params.fasta == '' ) {
    exit 1, "input missing, use [--fasta]"}

// nanopore reads input
/*
if (params.fasta && params.list) {
  Channel
        .fromPath( params.nano, checkIfExists: true )
        .splitCsv()
        .map { row -> ["${row[0]}", file("${row[1]}")] }
        .into { fastalist_report_ch; fastalist_analyse_ch }
fastalist_report_ch.subscribe { println "Got nano reads: ${it}" }
}
else if (params.fasta) {  }*/

fasta_input_ch = Channel.fromPath( params.fasta , checkIfExists: true)
                        .map { file -> tuple(file.simpleName, file) }
                        .view()
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
if (params.sourmeta) { include 'modules/sourmeta' params(output: params.output, cpus: params.cpus)
    sourmashmeta(fasta_input_ch,database_sourmash) }

/*************  
* --sourclass | Taxonomic classification via sourmash
*************/
if (params.sourclass) { include 'modules/sourclass' params(output: params.output, cpus: params.cpus) 
    sourmashclass(fasta_input_ch,database_sourmash) }



/*************  
* --help
*************/
def helpMSG() {
    log.info """

    Usage:
    nextflow run analysis.nf --fasta '*.fasta' [Workflows] [Options]

    Mandatory 
    --fasta             assembly file or files
 
    Options:
    --cores             max cores for local use [default: $params.cores]
    --output            name of the result folder [default: $params.output]

    Workflows:
    --sourmeta          Metagenomic sourmash analysis
    --sourclass         Taxonomic sourmash classification

    Nextflow options:
    -with-report rep.html    **not working** needs "ps" in all containers. (CPU and RAM usage report in rep.html)
    -with-dag chart.html     generates a flowchart for the process tree
    -with-timeline time.html **untested** generates a timeline

    Profile:
    -profile                 standard, gcloud (wip) [default: standard]
    """.stripIndent()
}