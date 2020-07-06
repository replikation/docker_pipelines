#!/usr/bin/env nextflow
nextflow.preview.dsl=2

/*
* Nextflow -- Analysis Pipeline
* Author: christian.jena@gmail.com
*/

/************************** 
* HELP messages & USER INPUT checks
**************************/
if (params.help) { exit 0, helpMSG() }

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

/************* 
* ERROR HANDLING
*************/
// profiles
if ( workflow.profile == 'standard' ) { "Using default profile [-profile local,docker]" }

// params tests
if (params.profile) {
    exit 1, "--profile is WRONG use -profile" }
if (!params.fasta &&  !params.fastq &&  !params.dir && !params.fastqPair && !params.dev) {
    exit 1, "input missing, use [--fasta] [--fastq] or [--dir]"}
if (params.fasta && params.fastq) {
    exit 1, "please us either: [--fasta] or [--fastq]"}   
if (params.fastq && params.metamaps && params.tax_db == '') {
    exit 1, "taxonomic database location not specified via [--tax_db]"}

if (params.watchFast5 && !params.samplename ) {
    exit 1, "please specify a sample name via [--samplename]"}

// fasta input or via csv file
if (params.fasta && params.list) { fasta_input_ch = Channel
        .fromPath( params.fasta, checkIfExists: true )
        .splitCsv()
        .map { row -> [row[0], file("${row[1]}", checkIfExists: true)]  }
        }
else if (params.fasta) { fasta_input_ch = Channel
        .fromPath( params.fasta, checkIfExists: true)
        .map { file -> tuple(file.baseName, file) }
        }

// fastq input or via csv file
if (params.fastq && params.list) { fastq_input_ch = Channel
        .fromPath( params.fastq, checkIfExists: true )
        .splitCsv()
        .map { row -> tuple("${row[0]}", file("${row[1]}", checkIfExists: true))  }
        }
else if (params.fastq) { fastq_input_ch = Channel
        .fromPath( params.fastq, checkIfExists: true)
        .map { file -> tuple(file.baseName, file) }
        }

// dir input or via csv file
if (params.dir && params.list) { dir_input_ch = Channel
        .fromPath( params.dir, checkIfExists: true )
        .splitCsv()
        .map { row -> [row[0], file("${row[1]}", checkIfExists: true , type: 'dir')] }
        .view() }
else if (params.dir) { dir_input_ch = Channel
        .fromPath( params.dir, checkIfExists: true, type: 'dir')
        .map { file -> tuple(file.name, file) }
        }

// illumina reads input & --list support
if (params.fastqPair && params.list) { fastqPair_input_ch = Channel
        .fromPath( params.fastqPair, checkIfExists: true )
        .splitCsv()
        .map { row -> [row[0], [file("${row[1]}", checkIfExists: true), file("${row[2]}", checkIfExists: true)]] }
        }
else if (params.fastqPair) { fastqPair_input_ch = Channel
        .fromFilePairs( params.fastqPair , checkIfExists: true )
        }

// live folder watching

if (params.watchFast5) { fast5_live_input_ch = Channel
   .watchPath( params.watchFast5 + '*.fast5' )
   .view()
}

if (params.samplename) { sample_name_ch = Channel.of ( params.samplename ) }

/************************** 
* DATABASES
**************************/
workflow sourmash_database_wf {
    main:
        if (params.sour_db) { database_sourmash = file(params.sour_db) }

        else if (workflow.profile == 'gcloud' && (params.sourmeta || params.sourclass)) {
                sour_db_preload = file("${params.database}/sourmash/gtdb.lca.json")    
            if (sour_db_preload.exists()) { database_sourmash = sour_db_preload }
            else {  sourmash_download_db() 
                    database_sourmash = sourmash_download_db.out } }
        else if (params.sourmeta || params.sourclass) {
                    sourmash_download_db() 
                    database_sourmash = sourmash_download_db.out }
    emit: database_sourmash
}

workflow metamaps_database_wf {
    main:
        if (params.tax_db) {
        database_metamaps = file(params.tax_db, checkIfExists: true, type: 'dir') }
    emit: database_metamaps
}

workflow gtdbtk_database_wf {
    main:
        if (params.gtdbtk_db) { database_gtdbtk = file(params.gtdbtk_db) }
        else if (workflow.profile == 'gcloud' && params.gtdbtk) {
                gtdbtk_preload = file("${params.database}/gtdbtk/gtdbtk_r89_data.tar.gz")    
            if (gtdbtk_preload.exists()) { database_gtdbtk = gtdbtk_preload }
            else {  gtdbtk_download_db() 
                    database_gtdbtk = gtdbtk_download_db.out } }
        else if (params.gtdbtk) {  
                    gtdbtk_download_db() 
                    database_gtdbtk = gtdbtk_download_db.out }
    emit: database_gtdbtk
}

workflow centrifuge_database_wf {
    main:
        if (params.centrifuge_db) { database_centrifuge = file( params.centrifuge_db ) }
        else if (!params.cloudProcess) { centrifuge_download_db() ; database_centrifuge = centrifuge_download_db.out}
        else if (params.cloudProcess) { 
            centrifuge_preload = file("${params.database}/centrifuge/gtdb_r89_54k_centrifuge.tar")
            if (centrifuge_preload.exists()) { database_centrifuge = centrifuge_preload }   
            else  { centrifuge_download_db()  ; database_centrifuge = centrifuge_download_db.out }
        }
    emit: database_centrifuge
} 

/************************** 
* MODULES
**************************/
    include abricate from './modules/abricate' 
    include abricateBatch from './modules/abricateBatch'
    include abricateParser from './modules/PARSER/abricateParser' 
    include abricateParserFASTA from './modules/PARSER/abricateParserFASTA'
    include abricatePlot from './modules/PLOTS/abricatePlot'
    include abricatePlotFASTA from './modules/PLOTS/abricatePlotFASTA' 
    include abricate_compare as abricate_chromosomes from './modules/abricate'
    include abricate_compare as abricate_plasmids from './modules/abricate'
    include abricate_compare as abricate_unknown from './modules/abricate'
    include abricate_transposon from './modules/abricate' 
    include baloonplot from './modules/PLOTS/baloonplot'
    include bwaUnmapped from './modules/bwaUnmapped' 
    include centrifuge from './modules/centrifuge' 
    include centrifuge_download_db from './modules/centrifugegetdatabase' 
    include centrifuge_illumina from './modules/centrifuge_illumina' 
    include chromomap from './modules/PLOTS/chromomap'
    include dev from './modules/dev' 
    include downloadHuman from './modules/downloadHuman' 
    include fargene as fargene_chromosomes from './modules/fargene'
    include fargene as fargene_plasmids from './modules/fargene'
    include fargene as fargene_unknown from './modules/fargene'
    include fargene_plasmid_screen from './modules/fargene'
    include fastqTofasta from './modules/fastqTofasta' 
    include fasttree from './modules/fasttree'
    include filter_fasta_by_length from './modules/filter_fasta_by_length'
    include filter_fastq_by_length from './modules/filter_fastq_by_length'
    include flye from './modules/flye'
    include gtdbtk from './modules/gtdbtk' 
    include gtdbtk_download_db from './modules/gtdbtkgetdatabase'
    include guppy_gpu from './modules/guppy_gpu' 
    include gviz from './modules/PLOTS/gviz'
    include krona from './modules/krona' 
    include live_guppy_gpu from './modules/guppy_gpu'
    include mafft from './modules/mafft'
    include mafft_supp from './modules/mafft_supp'
    include medaka from './modules/medaka' 
    include metamaps from './modules/metamaps' 
    include minimap2 from './modules/minimap2'
    include minimap2_polish from './modules/minimap2' 
    include nanoplot from './modules/nanoplot' 
    include overview_parser from './modules/PARSER/overview_parser'
    include parse_plasmidinfo from './modules/PARSER/parse_plasmidinfo' 
    include parse_prokka from './modules/PARSER/parse_prokka'
    include parse_samtools from './modules/PARSER/parse_samtools' 
    include plasflow from './modules/plasflow' 
    include plasflow_compare from './modules/plasflow' 
    include prokka from './modules/prokka' 
    include racon from './modules/racon' 
    include removeViaMapping from './modules/removeViaMapping' 
    include rmetaplot from './modules/rmetaplot' 
    include samtools from './modules/samtools'
    include sourclusterPlot from './modules/PLOTS/sourclusterPlot' 
    include sourmash_download_db from './modules/sourmashgetdatabase'
    include sourmashclassification from './modules/sourclass' 
    include sourmashclusterdir from './modules/sourclusterdir' 
    include sourmashclusterfasta from './modules/sourclusterfasta' 
    include sourmashmeta from './modules/sourmeta' 
    include toytree from './modules/toytree'

/************************** 
* SUB WORKFLOWS
**************************/

workflow centrifuge_wf {
    take:   fastq_input_ch
            centrifuge_DB
    main:   centrifuge(fastq_input_ch,centrifuge_DB) 
}

workflow centrifuge_illumina_wf {
    take:   fastqPair_input_ch
            centrifuge_DB
    main:   centrifuge_illumina(fastqPair_input_ch,centrifuge_DB) 
}

workflow guppy_gpu_wf {
    take:   dir_input_ch
    main:   guppy_gpu(dir_input_ch)
}

workflow deepHumanPathogen_wf {
    take:   fastqPair_input_ch
    main:   removeViaMapping(bwaUnmapped(fastqPair_input_ch,downloadHuman()))
}
   
workflow nanoplot_wf {
    take:   fastq_input_ch
    main:   nanoplot(fastq_input_ch)
}

workflow plasflow_wf {
    take:    fasta_input_ch
    main:   plasflow(fasta_input_ch)
}

workflow abricate_FASTQ_wf {
    take:    fastq_input_ch
    main:   method = ['ncbi', 'plasmidfinder']
            abricateBatch(fastqTofasta(fastq_input_ch.splitFastq(by: 100000, file: true)), method) 
            abricateBatch.out.collectFile(storeDir: "${params.output}/abricate-batch", skip: 1, keepHeader: true)
                collectResults = abricateBatch.out.collectFile(skip: 1, keepHeader: true).map { file -> tuple(file.baseName, file)}
            abricatePlot(abricateParser(collectResults))
}

workflow abricate_FASTA_wf {
    take:   fasta_input_ch
    main:   method = ['argannot', 'card', 'ncbi', 'plasmidfinder', 'resfinder']
            abricatePlotFASTA(abricateParserFASTA(abricate(fasta_input_ch, method)))
}

workflow abricate_FASTA_transposon_wf {
    take:   fasta_input_ch
    main:   mobile_database = Channel.fromPath( workflow.projectDir + "/data/IS.fna", checkIfExists: true )
            abricate_transposon(fasta_input_ch.combine(mobile_database))
    emit:   abricate_transposon.out // used in plasmid_comparision_wf {}
}

workflow gtdbtk_wf {
    take:   dir_input_ch
            gtdbtk_DB
    main:   gtdbtk(dir_input_ch,gtdbtk_DB)
}

workflow metamaps_wf {
    take:   fastq_input_ch
            metamaps_DB
    main:   krona(metamaps(fastq_input_ch, metamaps_DB))
}

workflow sourmash_WIMP_FASTA_wf {
    take:    fasta_input_ch
            sourmash_DB
    main:   sourmashmeta(fasta_input_ch,sourmash_DB)
}

workflow sourmash_WIMP_FASTQ_wf {
    take:   fastq_input_ch
            sourmash_DB
    main:   rmetaplot(sourmashmeta(fastqTofasta(fastq_input_ch),sourmash_DB))
}

workflow sourmash_tax_classification_wf {
    take:   fasta_input_ch
            sourmash_DB
    main:   sourmashclassification(fasta_input_ch,sourmash_DB)
}

workflow dev_build_centrifuge_DB_cloud_wf {
    main:       
    //repeater = ['8', '16', '24', '32', '40', '48']
        databasefile = file("gs://databases-nextflow/databases/thinspace/4centrifuge.tar.gz")
        dev(databasefile)
}

workflow sourmash_CLUSTERING_FASTA_wf {
    take:   fasta_input_ch
    main:   sourmashclusterfasta(fasta_input_ch)
}

workflow sourmash_CLUSTERING_DIR_wf {
    take:   fastq_input_ch
    main:   sourmashclusterdir(dir_input_ch)
            sourclusterPlot(sourmashclusterdir.out[1])
}

workflow amino_acid_tree_wf {
    take:   dir_input_ch
    main:   toytree(fasttree(mafft(dir_input_ch)))
}

workflow amino_acid_tree_supp_wf {
    take:   dir_input_ch
            proteins
    main:   toytree(fasttree(mafft_supp(dir_input_ch, proteins.map{ it -> it[1]})))
}

workflow resistance_comparision_wf {
  take: 
    fastas    //val(name), path(file))
  main:
    input_ch_plasflow = fastas.splitFasta(by: 50000, file: true)
                        .map { it -> tuple ( it[0], file(it[1]).getName(), it[1] ) }

    plasflow_compare( input_ch_plasflow ) 
    
    // abricate
    method = ['ncbi']
    abricate_chromosomes(plasflow_compare.out.genome, method)
    abricate_plasmids(plasflow_compare.out.plasmids, method)
    abricate_unknown(plasflow_compare.out.unclassified, method)

    // fargene
    hmm_method = ['class_a', 'class_b_1_2', 'class_b_3', 'class_c', 'class_d_1', 'class_d_2']
    fargene_chromosomes(plasflow_compare.out.genome, hmm_method)
    fargene_plasmids(plasflow_compare.out.plasmids, hmm_method)
    fargene_unknown(plasflow_compare.out.unclassified, hmm_method) 

    //summarize this visualy
    overview_parser(    abricate_chromosomes.out.map{ it -> it [1]}
                        .mix(abricate_plasmids.out.map{ it -> it [1]})
                        .mix(abricate_unknown.out.map{ it -> it [1]})
                        .mix(fargene_chromosomes.out)
                        .mix(fargene_plasmids.out)
                        .mix(fargene_unknown.out)
                        .collect(),
                        fastas.map { samplenames -> samplenames[0] }.collect() 
                    )

    baloonplot(overview_parser.out)

} 

workflow plasmid_comparision_wf {
  take: 
    fastas    //val(name), path(file)
  main:
    input_ch_plasflow = fastas.splitFasta(by: 50000, file: true)
                        .map { it -> tuple ( it[0], file(it[1]).getName(), it[1] ) }

    plasflow_compare( input_ch_plasflow )

    // abricate
    method = ['ncbi', 'plasmidfinder']
    hmm_method = ['class_a', 'class_b_1_2', 'class_b_3', 'class_c', 'class_d_1', 'class_d_2']

    abricate_plasmids(plasflow_compare.out.plasmids, method)      // *.abricate
    abricate_FASTA_transposon_wf(plasflow_compare.out.plasmids.map { it -> [it[0], it[3]] })  // *.tab
    fargene_plasmid_screen(plasflow_compare.out.plasmids.map { it -> [it[0], it[3]] }, hmm_method) // *.fargene

    // get fargene contigs and annotate them !
    // prokka(samtools(fargene_plasmid_screen.out.groupTuple()))

    group_by_sample = abricate_plasmids.out.groupTuple()
                      .join(abricate_FASTA_transposon_wf.out.groupTuple())
                      .join(fargene_plasmid_screen.out.groupTuple())
                      //.join(prokka.out.groupTuple())
 
    chromomap(parse_samtools(parse_plasmidinfo(group_by_sample).join(fastas)))  
}

workflow assembly_ont_wf {
    take:   fastq
    main:   medaka(racon(minimap2_polish(flye(fastq))))
    emit:   medaka.out
}


/************************** 
* Work in Progress section
**************************/
// Not sure about this one: its mainly implemented in the other workflow
// prokka is hard to parse here
workflow plasmid_annotate_wf {
    take: 
        fastas    //val(name), path(file)
    main:
        filter_fasta_by_length(fastas)
        input_ch_plasflow = filter_fasta_by_length.out
                            .map { it -> tuple ( it[0], file(it[1]).getName(), it[1] ) }

        plasflow_compare( input_ch_plasflow )
        prokka(plasflow_compare.out.plasmids.map { it -> [it[0], it[3]] })      // *.abricate

        group_by_sample = prokka.out.groupTuple()
    
        chromomap(parse_samtools(parse_prokka(group_by_sample).join(fastas)))  
}

// TODO: fastq files are not correctly stored, its just one - some "overwrite" bug i guess??
// could be that my links have the wrong name or so ?
// txt was working so you could add the PWD hast
workflow live_analysis_wf {
    take: 
        sample_name
        fast5_files
        reference_fasta
    main:
        
        live_guppy_gpu(sample_name.combine(fast5_files))
        
        //gviz(
            samtools(
                minimap2(
                    live_guppy_gpu.out.combine(reference_fasta.map { it -> it[1]})
        )) //)

    emit:
    live_guppy_gpu.out
}



/************************** 
* MAIN WORKFLOW
**************************/

workflow {
    if (params.abricate && params.fasta) { abricate_FASTA_wf(fasta_input_ch) }
    if (params.abricate && params.fastq) { abricate_FASTQ_wf(fastq_input_ch) }
    if (params.centrifuge && params.fastq) { centrifuge_wf(fastq_input_ch, centrifuge_database_wf()) }
    if (params.centrifuge && params.fastqPair) { centrifuge_illumina_wf(fastqPair_input_ch, centrifuge_database_wf()) }
    if (params.deepHumanPathogen && params.fastqPair) { deepHumanPathogen_wf(fastqPair_input_ch)}
    if (params.dev ) { dev_build_centrifuge_DB_cloud_wf() }
    if (params.gtdbtk && params.dir) { gtdbtk_wf(dir_input_ch,gtdbtk_database_wf()) }
    if (params.guppygpu && params.dir) { guppy_gpu_wf(dir_input_ch) }
    if (params.metamaps && params.fastq) { metamaps_wf(fastq_input_ch,metamaps_database_wf()) }
    if (params.mobile && params.fasta) { abricate_FASTA_transposon_wf(fasta_input_ch) }
    if (params.nanoplot && params.fastq) { nanoplot_wf(fastq_input_ch) }   
    if (params.plasflow && params.fasta) { plasflow_wf(fasta_input_ch) }
    if (params.plasmid_analysis && params.fasta) { plasmid_comparision_wf(fasta_input_ch) }
    if (params.plasmid_annotate && params.fasta) { plasmid_annotate_wf(fasta_input_ch) }
    if (params.res_compare && params.fasta) { resistance_comparision_wf(fasta_input_ch) }
    if (params.sourclass && params.fasta) { sourmash_tax_classification_wf(fasta_input_ch, sourmash_database_wf()) }
    if (params.sourcluster && params.dir ) { sourmash_CLUSTERING_DIR_wf(dir_input_ch) }
    if (params.sourcluster && params.fasta) { sourmash_CLUSTERING_FASTA_wf(fasta_input_ch) }
    if (params.sourmeta && params.fasta) { sourmash_WIMP_FASTA_wf(fasta_input_ch, sourmash_database_wf()) }
    if (params.sourmeta && params.fastq) { sourmash_WIMP_FASTQ_wf(fastq_input_ch, sourmash_database_wf()) }
    if (params.tree_aa && params.dir && !params.fasta) { amino_acid_tree_wf(dir_input_ch) }
    if (params.tree_aa && params.dir && params.fasta) { amino_acid_tree_supp_wf(dir_input_ch, fasta_input_ch) }
    if (params.assembly_ont && params.fastq) { assembly_ont_wf(fastq_input_ch) }

    // live workflows
    if (params.watchFast5 && params.samplename && params.fasta) { live_analysis_wf(sample_name_ch, fast5_live_input_ch, fasta_input_ch) }
    //if (params.dir_input_ch && params.samplename && params.fasta) { live_analysis_wf(sample_name_ch, fast5_live_input_ch, fasta_input_ch) }
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
    nextflow run replikation/docker_pipelines --fasta '*/*.fasta' --sourmeta --sourclass -profile local,docker 

    ${c_yellow}Inputs:${c_reset}
    ${c_green} --fasta ${c_reset}            '*.fasta'   -> assembly file(s) - uses filename
    ${c_green} --fastq ${c_reset}            '*.fastq'   -> read file(s) in fastq, one sample per file - uses filename
    ${c_green} --fastqPair ${c_reset}        '*_R{1,2}.fastq.gz'   -> fastq file pairs
    ${c_green} --dir  ${c_reset}             'foobar*/'  -> folder(s) as input - uses dirname
    ${c_dim}  ..change above input to csv:${c_reset} ${c_green}--list ${c_reset}${c_dim} e.g. --fasta list_files.csv --list            

    ${c_yellow}Resistance Workflows:${c_reset}
    ${c_blue} --abricate ${c_reset}          antibiotic and plasmid screening    ${c_green}[--fasta]${c_reset} or ${c_green}[--fastq]${c_reset}
    ${c_blue} --mobile ${c_reset}            screens for IS elements             ${c_green}[--fasta]${c_reset}
    ${c_blue} --res_compare ${c_reset}       detailed assembly resistance comparision of 2 or more assemblies ${c_green} [--fasta]${c_reset}
    ${c_dim}  ..option flags:            [--coverage] use coverage info in fasta headers on last position e.g. > name_cov_9.3354 ${c_reset}
    ${c_blue} --plasmid_analysis ${c_reset}  analysis of plasmids with plots     ${c_green}[--fasta]${c_reset}

    ${c_yellow}Cluster and Classifications:${c_reset} 
    ${c_blue} --sourmeta ${c_reset}          metagenomic analysis "WIMP"         ${c_green}[--fasta]${c_reset} or ${c_green}[--fastq]${c_reset}
    ${c_dim}  ..option flags:            [--sour_db] path to your own DB instead ${c_reset}
    ${c_blue} --sourclass ${c_reset}         taxonomic classification            ${c_green}[--fasta]  ${c_reset}
    ${c_dim}  ..option flags:            [--sour_db] path to your own DB instead ${c_reset}
    ${c_blue} --sourcluster ${c_reset}       sequence comparision with kmers     ${c_green}[--fasta]${c_reset} or ${c_green}[--dir]${c_reset}
    ${c_dim}  ..inputs:                  cluster contigs: [--fasta]; cluster fastas: [--dir]${c_reset}
    ${c_dim}  ..option flags:            [--size] figure size; default [--size $params.size]${c_reset}
    ${c_blue} --gtdbtk ${c_reset}            tax. class. via marker genes        ${c_green}[--dir]${c_reset}
    ${c_dim}  ..option flags:            [--gtdbtk_db] path to your own DB instead ${c_reset}

    ${c_yellow}Metagenomic Workflows:${c_reset}
    ${c_blue} --centrifuge ${c_reset}        metagenomic classification of reads ${c_green}[--fastq]${c_reset} or ${c_green}[--fastqPair]${c_reset}
    ${c_dim}  ..option flags:            [--centrifuge_db] path to your own DB instead, either .tar or .tar.gz ${c_reset}
    ${c_blue} --metamaps ${c_reset}          metagenomic class. of long reads    ${c_green}[--fastq]${c_reset}
    ${c_dim}  ..mandatory flags:         [--memory] [--tax_db] e.g. --memory 100 --tax_db /databases/miniSeq+H 

    ${c_yellow}Nanopore specific Workflows:${c_reset}
    ${c_blue} --guppygpu ${c_reset}          basecalling via guppy-gpu-nvidia   ${c_green} [--dir]${c_reset}
    ${c_dim}  ..option flags:            [--flowcell] [--kit] [--barcode] [--modbase]
        ..default settings:        [--flowcell $params.flowcell] [--kit $params.kit] [--modbase FALSE] ${c_reset}
    ${c_dim}  ..config files:            turn on via [--config], modify config type via [--configtype] 
        ..default config type:     [--configtype $params.configtype]  ${c_reset}
    ${c_blue} --nanoplot  ${c_reset}         read quality via nanoplot           ${c_green}[--fastq]${c_reset}
    ${c_blue} --assembly_ont ${c_reset}      simple nanopore assembly            ${c_green}[--fastq]${c_reset}
    ${c_dim}  ..option flags:            [--gsize ${params.gsize}] [--model ${params.model}] [--overlap ${params.overlap}]

    ${c_yellow}Nanopore live analysis Workflows (WIP):${c_reset}
    ${c_blue} --watchFast5 ${c_reset}    watch a dir for fast5 files, basecall them and map against reference
                      Needs: ${c_green}[--samplename]${c_reset} and one multifasta file via ${c_green}[--fasta]${c_reset}

    ${c_yellow}Nanopore Live analysis **WIP**${c_reset}
    [--watchFast5]          directory where fast5 files appear ${c_green} [--watchFast5 fast5/]${c_reset}
    [--samplename]          name of your sample ${c_green} [--samplename "E.coli"]${c_reset}
    [--fasta]               reference multi fastas file to screen your reads against  [--fasta some_genomes.fasta]${c_reset}

    ${c_yellow}Other Workflows:${c_reset}
    ${c_blue} --deepHumanPathogen ${c_reset} pathogen identification in human    ${c_green}[--fastqPair '*_R{1,2}.fastq.gz']${c_reset}    
    ${c_blue} --plasflow ${c_reset}          predicts & seperates plasmid-seqs   ${c_green}[--fasta]${c_reset}
    ${c_blue} --tree_aa ${c_reset}           aminoacid tree of a dir with aa seq ${c_green}[--dir]${c_reset}
    ${c_dim}  ..option flags:            [--filenames] use filenames as labels instead of contig names
    ${c_dim}                             [--fasta] add one multi protein file as "tree enhancer" 
    ${c_dim}                               e.g. [--fasta multipleProteins.aa]${c_reset}

    ${c_reset}Options:
    --cores             max cores for local use [default: $params.cores]
    --memory            80% of available RAM in GB for --metamaps [default: $params.memory]
    --output            name of the result folder [default: $params.output]
    --workdir           location of temporaty files [default: $params.workdir]

    ${c_reset}Database(s) behaviour${c_reset}
    The Priority is:
    1. Use your own DB via the flags, e.g. [--sour_db] [--gtdbtk_db]
    2. Without a flag it downloads/retrives a database to/from: ./nextflow-autodownload-databases
    ${c_dim}If a auto download is not possible or implemented the workflow will tell you

    ${c_dim}Nextflow options:
    -with-report rep.html    cpu / ram usage (may cause errors)
    -with-dag chart.html     generates a flowchart for the process tree
    -with-timeline time.html timeline (may cause errors)

    Profile:
    -profile                 local gcloud docker -> merge profiles e.g. -profile local,docker ${c_reset}
    """.stripIndent()
}
