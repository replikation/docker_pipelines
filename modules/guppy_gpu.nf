process guppy_gpu {
      echo true
      maxForks 1
      container = 'nanozoo/guppy_gpu:3.4.4-1--3dd30a4'
      containerOptions '--gpus all'
      publishDir "${params.output}/${name}/", mode: 'copy', pattern: "fastq_${params.output}"
    input:
      tuple val(name), file(dir)
    output:
      tuple val(name), file("fastq_${params.output}")
    script:
      if (params.config)
      """
      guppy_basecaller -r -i ${dir} -s fastq_${params.output} \
      -c ${params.configtype} --device auto --trim_strategy dna -q 0
      """
      else if(!params.barcode)
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

// dna_r9.4.1_450bps_modbases_dam-dcm-cpg_hac.cfg
