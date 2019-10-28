process abricateBatch {
    label 'abricate'   
  input:
    tuple val(name), file(fasta) 
    each method
  output:
	  tuple val(name), file("*.tab") 
  script:
    """
  	abricate ${fasta} --nopath --quiet --mincov 80 --db ${method} > ${method}.tab
    """
}

// tempDir may be neccessary