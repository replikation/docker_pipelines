process abricateBatch {
    label 'abricate'   
  input:
    set val(name), file(fasta) 
    each method
  output:
	  set val(name), file("*.tab") 
  script:
    """
  	abricate ${fasta} --nopath --quiet --mincov 80 --db ${method} > ${method}.tab
    """
}

// tempDir may be neccessary