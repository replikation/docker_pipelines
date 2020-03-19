
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



/*
process minimap2_to_bam introduces a Math.random() function to add a random long number to the output.
This is introduced as a quick fix, because otherweise we create to much files of the same name which will be a issue in
the process later on when we use all of them together (differential binning)
*/