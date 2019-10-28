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
      sourmash compute -p ${task.cpus} --scaled 10000 -k 31 ${fasta} -o ${name}.sig 
      sourmash lca gather  ${name}.sig ${database} --ignore-abundance -o metagenomic-composition.txt
      """
      else if (params.fastq)
      """
      sourmash compute -p ${task.cpus} --scaled 1000 --track-abundance -k 31 ${fasta} -o ${name}.sig 
      sourmash lca gather  ${name}.sig ${database} -o metagenomic-composition.txt
      """
}