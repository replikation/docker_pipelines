process clinker {
      publishDir "${params.output}/", mode: 'copy', pattern: "results.html"
      label 'clinker'
    input:
      path(genbankfile) 
    output:
      path("results.html")
    script:
      """
        clinker *.gbk --jobs 20 --plot results.html --identity 0.7
      """
}