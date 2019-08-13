process basecalling {
    label 'basecalling'
      maxForks 1
      publishDir "${params.output}/${name}/", mode: 'copy', pattern: "fastq_${params.output}"
    input:
      set val(name), file(dir)
    output:
      set val(name), file("fastq_${params.output}")
    script:
      if(params.barcode == '')
      """
      guppy_basecaller -r -i ${dir} -s fastq_${params.output} \
      --flowcell ${flowcell} --kit ${kit} --device auto --trim_strategy dna -q 0
      """
      else
      """
      guppy_basecaller -r -i ${dir} -s fastq_${params.output} \
      --flowcell ${flowcell} --kit ${kit} --barcode_kits ${barcode} --device auto --trim_strategy dna -q 0 \
      --trim_barcodes
      """
}