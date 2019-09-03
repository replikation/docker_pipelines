process rmetaplot {
      publishDir "${params.output}/${name}", mode: 'copy', pattern: "Metagenomic-distribution.pdf"
      label 'rmetaplot'
    input:
      set val(name), file(sourmeta) 
    output:
      set val(name), file("Metagenomic-distribution.pdf")
    script:
      """
      shell-extract.sh
      dot-diagram.R summary.results superkindoms.results 
      """

}