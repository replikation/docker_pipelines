process plasflow {
      publishDir "${params.output}/${name}/chromosome", mode: 'copy', pattern: "${name}_chromosomes.fasta"
      publishDir "${params.output}/${name}/plasmids", mode: 'copy', pattern: "${name}_plasmids.fasta"
      publishDir "${params.output}/${name}/unclassified", mode: 'copy', pattern: "${name}_unclassified.fasta"
      label 'plasflow'
    input:
      tuple val(name), file(fasta) 
    output:
      tuple val(name), file("${name}_chromosomes.fasta"), file("${name}_plasmids.fasta"), file("${name}_unclassified.fasta")
    script:
      """
      PlasFlow.py --input ${fasta} --output ${name} --threshold 0.7
      """
}


process plasflow_compare {
      publishDir "${params.output}/plasflow/${name}/chromosome", mode: 'copy', pattern: "${name}_${splitname}_chromosomes.fasta"
      publishDir "${params.output}/plasflow/${name}/plasmids", mode: 'copy', pattern: "${name}_${splitname}_plasmids.fasta"
      publishDir "${params.output}/plasflow/${name}/unclassified", mode: 'copy', pattern: "${name}_${splitname}_unclassified.fasta"
      label 'plasflow'
      //maxForks 4
    input:
      tuple val(name), val(splitname), path(fasta) 
    output:
      tuple val(name), val(splitname), val("chromosome"), path("${name}_${splitname}_chromosomes.fasta"), emit: genome, optional: true
      tuple val(name), val(splitname), val("plasmid"), path("${name}_${splitname}_plasmids.fasta"), emit: plasmids, optional: true
      tuple val(name), val(splitname), val("unclassified"), path("${name}_${splitname}_unclassified.fasta"), emit: unclassified, optional: true
    script:
      """
      PlasFlow.py --input ${fasta} --output ${name}_${splitname} --threshold 0.7
      
      for fastafile in ${name}_*.fasta; do
        lines=\$(cat \${fastafile} | wc -l)
        if [ \${lines} -lt 2 ]; then rm -f \${fastafile}; fi
      done
      """
}