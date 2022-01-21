include { spades } from './process/spades.nf'

workflow spades_wf {
    take: 
        fastqPair //tuple val(fasta-basename) path(fastq-files)
    main:
        spades(fastqPair)
}