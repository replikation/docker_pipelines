process sourmashclassification {
      publishDir "${params.output}/${name}", mode: 'copy', pattern: "taxonomic-classification.txt"
      label 'sourmash'
    input:
      tuple val(name), file(fasta) 
      file(database) 
    output:
      tuple val(name), file("taxonomic-classification.txt")
    script:
      """
      sourmash sketch dna -p scaled=10000,k=31 ${fasta} -o ${name}.sig
      sourmash lca classify --db ${database} --query ${name}.sig -o taxonomic-classification.txt
      """
}