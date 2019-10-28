process basecalling {
      echo true
      maxForks 1
      container = 'nanozoo/guppy_gpu:3.2.4-1--195590e'
      containerOptions '--gpus all'
      publishDir "${params.output}/${name}/", mode: 'copy', pattern: "fastq_${params.output}"
    input:
      tuple val(name), file(dir)
    output:
      tuple val(name), file("fastq_${params.output}")
    script:
      if(params.barcode == '')
      """
      guppy_basecaller -r -i ${dir} -s fastq_${params.output} \
      --flowcell ${params.flowcell} --kit ${params.kit} --device auto --trim_strategy dna -q 0
      """
      else
      """
      guppy_basecaller -r -i ${dir} -s fastq_${params.output} \
      --flowcell ${params.flowcell} --kit ${params.kit} --barcode_kits ${params.barcode} --device auto --trim_strategy dna -q 0 \
      --trim_barcodes
      """
}