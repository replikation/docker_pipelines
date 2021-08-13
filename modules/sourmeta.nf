process sourmashmeta {
      publishDir "${params.output}/${name}", mode: 'copy', pattern: "metagenomic-composition.txt"
      label 'sourmash'
    input:
      tuple val(name), file(fasta) 
      file(database) 
    output:
      tuple val(name), file("metagenomic-composition.txt")
    script:
      if (params.fasta)
      """
      sourmash sketch dna -p scaled=10000,k=31 ${fasta} -o ${name}.sig 
      sourmash gather  ${name}.sig ${database} --ignore-abundance -o metagenomic-composition.txt
      """
      else if (params.fastq)
      """
      sourmash sketch dna -p scaled=10000,k=31 --track-abundance ${fasta} -o ${name}.sig 
      sourmash gather  ${name}.sig ${database} -o metagenomic-composition.txt
      """
}