include { kraken2_illumina_pe } from './process/kraken2.nf' 
include { krona; krona_from_bracken } from './process/krona.nf' 
include { download_database_kraken2 } from './process/download_database_kraken2.nf'
include { bracken } from './process/bracken.nf'
include { krakentools } from './process/krakentools.nf'

workflow read_classification_illumina_pe_wf {
    take:   
        fastq
    main: 

        // database download
        if (params.krakendb) { kraken_db = file("${params.krakendb}") }
        else  { download_database_kraken2(); kraken_db = download_database_kraken2.out } 

        // classification
        kraken2_illumina_pe(fastq, kraken_db)

        // alpha diversity, abundance and korna plots
        krona_from_bracken(krakentools(bracken(kraken2_illumina_pe.out, kraken_db)))

    emit:   
        kraken = kraken2_illumina_pe.out
}





/*
Protocoll here


you might want to add the option to calculate beta diversity here?


https://www.nature.com/articles/s41596-022-00738-y

*/
