process rmetaplot {
      publishDir "${params.output}/${name}", mode: 'copy', pattern: "Metagenomic-distribution.pdf"
      label 'rmetaplot'
    input:
      tuple val(name), file(sourmeta) 
    output:
      tuple val(name), file("Metagenomic-distribution.pdf")
    script:
      """
      shell-extract.sh
      dot-diagram.R summary.results superkindoms.results 
      """

}