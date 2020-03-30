process samtools {
  label 'samtools'
  publishDir "${params.output}/", mode: 'copy'
    input:
      tuple val(name) , file(bamfile)
    output:
      tuple val(name), file("results_*.txt")
    shell:
    """
    samtools index ${bamfile}
    samtools idxstats ${bamfile} > results_\${PWD##*/}.txt
    """
}




