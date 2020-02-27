process baloonplot {
    publishDir "${params.output}/", mode: 'copy'
    label 'baloonplot'   
  input:
    path(csv) 
  output:
	  path("overview.svg") 
  script:
    """
    baloonrplot.R
    """
}

// tempDir may be neccessary