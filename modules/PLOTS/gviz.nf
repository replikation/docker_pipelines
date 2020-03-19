process gviz {
    publishDir "${params.output}/", mode: 'copy'
    label 'gviz' 
  input:
    path(csv) 
  output:
	  path("overview.svg") 
  script:
    """
    gviz.R
    """

}

