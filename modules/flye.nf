process flye {
    label 'flye'  
  input:
    tuple val(name), path(read)
  output:
    tuple val(name), path(read), path("${name}.fasta")
  script:
    """
    flye --plasmids --meta -t ${task.cpus} --nano-raw ${read} -o assembly
    mv assembly/assembly.fasta ${name}.fasta --min-overlap ${params.overlap}
    """
  }