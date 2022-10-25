process centrifuge_illumina {
      publishDir "${params.output}/${name}/centrifuge", mode: 'copy', pattern: "${name}_pavian_report_filtered.csv"
      publishDir "${params.output}/${name}/centrifuge", mode: 'copy', pattern: "${name}.out"
      label 'centrifuge'


    input:
      tuple val(name), file(fastq) 
      path(database) 
    output:
      tuple val(name), file("${name}.out"), file("${name}_pavian_report_filtered.csv")
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

      centrifuge -p ${task.cpus} -x \${DBname} -k 5 \
        -1 ${fastq[0]} -2 ${fastq[0]} -S centrifuge_results.out --report-file centrifuge_out.log

      # filter based on score
      < centrifuge_results.out awk '{if(NR < 2 || \$4 >= 250) {print}}' | awk '{if(NR < 2 || \$6 >= 150) {print}}' > centrifuge_filtered.out

      centrifuge-kreport -x \${DBname} centrifuge_filtered.out > ${name}_pavian_report_filtered.csv

      mv centrifuge_filtered.out ${name}.out

      """
}




