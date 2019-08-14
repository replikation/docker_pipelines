process plasflow {
      publishDir "${params.output}/${name}/chromosome", mode: 'copy', pattern: "${name}_chromosomes.fasta"
      publishDir "${params.output}/${name}/plasmids", mode: 'copy', pattern: "${name}_plasmids.fasta"
      publishDir "${params.output}/${name}/unclassified", mode: 'copy', pattern: "${name}_unclassified.fasta"
      label 'plasflow'
    input:
      set val(name), file(fasta) 
    output:
      set val(name), file("${name}_chromosomes.fasta"), file("${name}_plasmids.fasta"), file("${name}_unclassified.fasta")
    script:
      """
      PlasFlow.py --input ${fasta} --output ${name} --threshold 0.7
      """
}