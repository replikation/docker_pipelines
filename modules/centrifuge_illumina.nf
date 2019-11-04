process centrifuge_illumina {
      publishDir "${params.output}/${name}/centrifuge", mode: 'copy', pattern: "${name}_pavian_report_filtered.csv"
      publishDir "${params.output}/${name}/centrifuge", mode: 'copy', pattern: "${name}.out"
      label 'centrifugescale'

      errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
      errorStrategy { task.exitStatus == 14 ? 1 : task.attempt }
      cpus { 12 * task.attempt }
      memory { 70.GB * task.attempt }
      maxRetries 2
      

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

      centrifuge-kreport -x \${DBname} centrifuge_results.out \
      > ${name}_pavian_report_filtered.csv

      mv centrifuge_results.out ${name}.out

      """
}




