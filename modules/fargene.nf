process fargene {
    label 'fargene'
    errorStrategy 'ignore'
  input:
    tuple val(name), val(splitname), val(type), path(fasta) 
    each method
  output:
    path("${method}_${type}_${name}_*.fargene") 
  script:
    """
    touch "${method}_${type}_${name}_\${PWD##*/}.fargene"
  	fargene -i ${fasta} --hmm-model ${method} -o ${name}_results 
    cp ${name}_results/results_summary.txt  ${method}_${type}_${name}_\${PWD##*/}.fargene
    """
}

process fargene_plasmid_screen {
    label 'fargene'
    validExitStatus 0,1
    //publishDir "${params.output}/fargene/", mode: 'copy'
  input:
    tuple val(name), path(fasta) 
    each method
  output:
    tuple val(name), path("${method}_${name}_*.fargene") optional true
  script:
    """
    filename=\$(ls ${fasta})
  	fargene -i ${fasta} --hmm-model ${method} -o ${name}_results 
    cat ${name}_results/hmmsearchresults/retrieved-genes-*-hmmsearched.out |  grep -v "^#" | sed -e s/"\${filename%.*}_"//g > ${method}_${name}_\${PWD##*/}.fargene
    """
}