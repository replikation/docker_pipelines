process sourmashmeta {
      publishDir "${params.output}/${name}", mode: 'copy', pattern: "metagenomic-composition.txt"
      label 'sourmash'
    input:
      set val(name), file(fasta) 
      file(database) 
    output:
      set val(name), file("metagenomic-composition.txt")
    script:
      """
      sourmash compute -p ${task.cpus} --scaled 10000 -k 31 ${fasta} -o ${name}.sig 
      sourmash lca gather  ${name}.sig ${database} --ignore-abundance -o metagenomic-composition.txt
      """
}