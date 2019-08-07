process nanoplot {
    label 'nanoplot'
      publishDir "${params.output}/${name}/", mode: 'copy', pattern: "*.html"
      publishDir "${params.output}/${name}/readquality/", mode: 'copy', pattern: "*_read_quality.txt"
      publishDir "${params.output}/${name}/readquality/", mode: 'copy', pattern: "*.png"
      publishDir "${params.output}/${name}/readquality/", mode: 'copy', pattern: "*.pdf"
    input:
      set val(name), file(reads)
    output:
      set val(name), file("*.html"), file("*.png"), file("*.pdf"), file("${name}_read_quality.txt") 
    script:
      """
      NanoPlot -t ${params.cpus} --fastq ${reads} --title ${name} --color darkslategrey --N50 --plots hex --loglength -f png --store
      NanoPlot -t ${params.cpus} --pickle NanoPlot-data.pickle --title ${name} --color darkslategrey --N50 --plots hex --loglength -f pdf
      mv *.html ${name}_read_quality_report.html
      mv NanoStats.txt ${name}_read_quality.txt
      """
}