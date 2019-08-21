process metamaps {
      label 'metamaps'
    input:
      set val(name), file(fastq) 
      file(database) 
    output:
      set val(name), file ("${name}.krona")
    script:
      if (workflow.profile == 'gcloud') {
      """
      metamaps mapDirectly -t ${task.cpus} --all -m 2000 -r ${database}/DB.fa -q ${fastq} -o classification_results
      metamaps classify -t ${task.cpus} --mappings classification_results --DB ${database}
      cp classification_results/*.krona ${name}.krona
      """
      }
      else {
      """
      metamaps mapDirectly -t ${task.cpus} --all -m 2000 -r ${database}/DB.fa -q ${fastq} -o classification_results --maxmemory ${params.memory}
      metamaps classify -t ${task.cpus} --mappings classification_results --DB ${database}
      cp classification_results/*.krona ${name}.krona
      """
      }
}