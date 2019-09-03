process fastqTofasta {
      publishDir "${params.output}/${name}/", mode: 'copy', pattern: "filtered_reads.fasta"
      label 'emboss'
    input:
      set val(name), file(fastq)
    output:
      set val(name), file("filtered_reads.fasta")
    script:
      """
      seqret -sequence ${fastq} -outseq filtered_reads.fasta
      """
}