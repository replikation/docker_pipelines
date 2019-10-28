process krona {
      publishDir "${params.output}/${name}", mode: 'copy', pattern: "${name}.html"
      label 'krona'
    input:
      tuple val(name), file (krona)
    output:
      tuple val(name), file("${name}.html")
    script:
      """
      ktImportTaxonomy classification_results.EM.reads2Taxon.krona
      mv taxonomy.krona.html ${name}.html
      """
}