process centrifuge {
      publishDir "${params.output}/${name}/centrifuge", mode: 'copy', pattern: "${name}_pavian_report_filtered.csv"
      publishDir "${params.output}/${name}/centrifuge/raw", mode: 'copy', pattern: "${name}_centrifuge_filtered.out"
      publishDir "${params.output}/${name}/centrifuge", mode: 'copy', pattern: "${name}_pavian_report_unfiltered.csv"
      publishDir "${params.output}/${name}/centrifuge/raw", mode: 'copy', pattern: "${name}_centrifuge_unfiltered.out"
      label 'centrifuge'
    input:
      tuple val(name), file(fastq) 
      path(database) 
    output:
      tuple val(name), file("${name}_centrifuge_filtered.out"), file("${name}_pavian_report_filtered.csv")
      tuple val(name), file("${name}_centrifuge_unfiltered.out"), file("${name}_pavian_report_unfiltered.csv")
    script:
      """
      case "${database}" in
      *.tar.gz)
        tar xzf ${database}
        ;;
      *.gz | *.tgz ) 
        gzip -d ${database}
        ;;
      *.tar)
        tar xf ${database}
        ;;
      esac
      
      DBname=\$(ls *.[1-9].cf | head -1 | cut -f1 -d".")

      # command
        centrifuge -p ${task.cpus} -x \${DBname} -k 5 --min-hitlen 16 \
        -U ${fastq} -S ${name}_centrifuge_unfiltered.out --report-file centrifuge_out.log


      # filter
        filter_centrifuge_hits.sh ${name}_centrifuge_unfiltered.out 150 250 ${name}_centrifuge_filtered.out

      # reports
        centrifuge-kreport -x \${DBname} --min-score 300 --min-length 500 ${name}_centrifuge_filtered.out \
        > ${name}_pavian_report_filtered.csv

        centrifuge-kreport -x \${DBname} --min-score 300 --min-length 500 ${name}_centrifuge_unfiltered.out \
        > ${name}_pavian_report_unfiltered.csv
      """
}




