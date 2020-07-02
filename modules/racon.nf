process racon {
      label 'racon'
   input:
      tuple val(name), path(read), path(assembly), path(mapping) 
   output:
   	tuple val(name), file(read), file("${name}_consensus.fasta") 
   shell:
      """
      racon -t ${task.cpus} ${read} ${mapping} ${assembly} > ${name}_consensus.fasta
      """
  }