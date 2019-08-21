process sourmashclusterdir {
      publishDir "${params.output}/${name}/cluster/", mode: 'copy', pattern: "*.pdf"
      label 'sourmash'
    input:
      set val(name), file(dir) 
    output:
      set val(name), file("*.pdf")
    script:
      """
      cp ${dir}/* .
      sourmash compute -p ${task.cpus} *.fa* --scaled 10000 -k 31
      sourmash compare *.sig -o results_sig
      sourmash plot --pdf --subsample=250 --labels results_sig
      """
}



