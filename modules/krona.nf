process krona {
      publishDir "${params.output}/${name}", mode: 'copy', pattern: "classification_results"
      label 'krona'
    input:
      set val(name), file (krona)
    output:
      set val(name), file("classification_results")
    script:
      """
      ktImportTaxonomy ${krona}
      mv taxonomy.krona.html ${name}.html
      """
}