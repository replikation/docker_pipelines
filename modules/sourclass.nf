process sourmashclass {
      publishDir "${params.output}/${name}", mode: 'copy', pattern: "taxonomic-classification.txt"
      label 'sourmash'
    input:
      set val(name), file(fasta) 
      file(database) 
    output:
      set val(name), file("taxonomic-classification.txt")
    script:
      """
      sourmash compute -p ${params.cpus} --scaled 10000 -k 31 ${fasta} -o ${name}.sig
      sourmash lca classify --db ${database} --query ${name}.sig -o taxonomic-classification.txt
      """
}