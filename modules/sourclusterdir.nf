process sourmashclusterdir {
      publishDir "${params.output}/${name}/cluster/", mode: 'copy', pattern: "*.pdf"
      publishDir "${params.output}/${name}/cluster/", mode: 'copy', pattern: "results.csv"
      label 'sourmash'
    input:
      tuple val(name), file(dir) 
    output:
      tuple val(name), file("*.pdf")
      tuple val(name), file("results.csv")
    script:
      """
      cp ${dir}/*.* .
      sourmash sketch dna -p scaled=10000,k=31 *.fa*
      sourmash compare *.sig -o results_sig --csv results.csv
      sourmash plot --pdf --subsample=250 --labels results_sig
      """
}



