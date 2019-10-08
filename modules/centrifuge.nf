process centrifuge {
      publishDir "${params.output}/${name}/centrifuge", mode: 'copy', pattern: "${name}_pavian_report_filtered.csv"
      publishDir "${params.output}/${name}/centrifuge", mode: 'copy', pattern: "${name}.out"
      label 'centrifuge'
    input:
      set val(name), file(fastq) 
      file(database) 
    output:
      set val(name), file("${name}.out"), file("${name}_pavian_report_filtered.csv")
    shell:
      """
      case "!{database}" in
      *.tar.gz)
        tar xzf !{database}
        ;;
      *.gz | *.tgz ) 
        gzip -d !{database}
        ;;
      *.tar)
        tar xf !{database}
        ;;
      esac
      
      DBname=\$(ls *.[1-9].cf | head -1 | cut -f1 -d".")

      centrifuge -p !{task.cpus} -x \${DBname} -k 5 --min-hitlen 16 \
      -U !{fastq} -S centrifuge_results.out --report-file centrifuge_out.log

      filter-reads-ont.sh

      centrifuge-kreport -x \${DBname} --min-score 300 --min-length 500 centrifuge_filtered.out \
      > !{name}_pavian_report_filtered.csv

      mv centrifuge_filtered.out !{name}.out
      """
}




