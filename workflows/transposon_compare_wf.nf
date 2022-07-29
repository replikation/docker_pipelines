include { abricate } from './process/abricate.nf'
include { abricate_samtool_extracter } from './process/abricate_samtool_extracter.nf'
include { prokka_gbk } from './process/prokka_gbk.nf'
include { clinker } from './process/clinker.nf'


workflow transposon_compare_wf {
take: fasta_files

main:

    abricate(fasta_files.map{it -> it[1]}.collect())
    abricate_samtool_extracter(abricate.out)
    prokka_gbk(abricate_samtool_extracter.out.flatten())
    clinker(prokka_gbk.out.collect())

}
