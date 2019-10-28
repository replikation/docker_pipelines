process sourmashclusterfasta {
      publishDir "${params.output}/${name}/cluster/", mode: 'copy', pattern: "*.pdf"
      label 'sourmash'
    input:
      tuple val(name), file(fasta) 
    output:
      tuple val(name), file("*.pdf")
    script:
      """
      sourmash compute -p ${task.cpus} --scaled 10000 -k 31 ${fasta} -o signature.sig --singleton
      sourmash compare signature.sig -o results_sig
      sourmash plot --pdf --subsample=250 --labels results_sig
      """
}



