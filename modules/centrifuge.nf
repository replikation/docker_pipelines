process centrifuge {
      publishDir "${params.output}/${name}/centrifuge", mode: 'copy', pattern: "${name}_pavian_report_filtered.csv"
      label 'centrifuge'
    input:
      set val(name), file(fastq) 
      file(database) 
    output:
      set val(name), file("${name}_pavian_report_filtered.csv")
    script:
      """
      tar xf ${database}
      centrifuge -p ${task.cpus} -x ex -k 5 --min-hitlen 16 \
      -U ${fastq} -S centrifuge_results.out --report-file centrifuge_out.log

      filter-reads-ont.sh

      centrifuge-kreport -x ex --min-score 300 --min-length 500 centrifuge_filtered.out \
      > ${name}_pavian_report_filtered.csv
      """
}