process fargene {
    label 'fargene'
    def random = (Math.random() + Math.random()).toString().md5().toString()
    errorStrategy 'ignore'
  input:
    tuple val(name), val(type), path(fasta) 
    each method
  output:
    path("${method}_${type}_${name}.fargene_${random}") 
  script:
    """
    touch "${method}_${type}_${name}.fargene_${random}"
  	fargene -i ${fasta} --hmm-model ${method} -o ${name}_results 
    cp ${name}_results/results_summary.txt  ${method}_${type}_${name}.fargene_${random} 
    """
}

