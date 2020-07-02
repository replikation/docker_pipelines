
process minimap2 {
  label 'minimap2'
    input:
      tuple val(name), path(reads), path(fastas)
    output:
      tuple val(name), path("ont_sorted.bam")
    shell:
    """
    minimap2 -ax map-ont ${fastas} ${reads} | samtools view -bS - | samtools sort -@ ${task.cpus} - > ont_sorted.bam
    """
}

process minimap2_polish {
  label 'minimap2'
      input:
  	    tuple val(name), file(read), file(assembly) 
      output:
        tuple val(name), file(read), file(assembly), file("${name}.paf") 
      script:
        """
      	minimap2 -x map-ont -t ${task.cpus} ${assembly} ${read} > ${name}.paf
        """
}