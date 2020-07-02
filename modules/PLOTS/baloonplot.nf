process baloonplot {
    publishDir "${params.output}/", mode: 'copy'
    label 'baloonplot' 
  input:
    path(csv) 
  output:
	  path("overview.svg")
    path("overview.pdf") 
  script:
  if (params.coverage)
    """
    baloonrplot_coverage_count.R
    """
  else if (!params.coverage)
    """
    baloonrplot_gene_count.R
    """
}

