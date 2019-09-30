process bwaUnmapped {
      publishDir "${params.output}/${name}/", mode: 'copy', pattern: "${name}_non_mappers1_ids.lst"
      label 'bwa' 
    input:
      set val(name), file(shortreads)
      file(assembly) 
    output:
      set val(name), file("${name}_non_mappers1_ids.lst"), file(shortreads)
    script:
      """
      bwa index ${assembly}
      bwa mem -t ${task.cpus} ${assembly} ${shortreads[0]} ${shortreads[1]} > mapped-genome.sam
      samtools view -S -f4 mapped-genome.sam > sample.unmapped.sam
      cut -f1 sample.unmapped.sam | sort | uniq > ${name}_non_mappers1_ids.lst
      """
}

/*
Comment:
This is a modified bwa, which serves to get the "unmapped" reads !! its not a standard bwa module


*/