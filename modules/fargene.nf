process fargene {
    label 'fargene'
  input:
    tuple val(name), val(type), path(fasta) 
    each method
  output:
    path("${method}_${type}_${name}_summary.fargene")
  script:
    """
  	fargene -i ${fasta} --hmm-model ${method} -o ${name}_results
    cp ${name}_results/results_summary.txt  ${method}_${type}_${name}_summary.fargene
    """
}

