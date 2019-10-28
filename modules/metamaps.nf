process metamaps {
      label 'metamaps'
      publishDir "${params.output}/${name}/metamaps", mode: 'copy', pattern: "classification_results.EM.*"
    input:
      tuple val(name), file(fastq) 
      file(database) 
    output:
      tuple val(name), file ("classification_results.EM.*")
    script:
      if (workflow.profile == 'gcloud') {
      """
      metamaps mapDirectly -t ${task.cpus} --all -m 1000 -r ${database}/DB.fa -q ${fastq} -o classification_results
      metamaps classify -t ${task.cpus} --mappings classification_results --DB ${database}
      """
      }
      else {
      """
      metamaps mapDirectly -t ${task.cpus} --all -m 1000 -r ${database}/DB.fa -q ${fastq} -o classification_results --maxmemory ${params.memory}
      metamaps classify -t ${task.cpus} --mappings classification_results --DB ${database}
      """
      }
}