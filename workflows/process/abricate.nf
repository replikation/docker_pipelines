process abricate {
    label 'abricate'
  input:
        path(fasta) 
  output:
	    tuple path("allABR.tab"), path("allfastas.fasta")
  script:
    if (!params.update)
    """
    cat ${fasta} > allfastas.fasta
  	abricate allfastas.fasta --nopath --quiet --mincov 85 --db ncbi > allABR.tab
    abricate allfastas.fasta --nopath --quiet --mincov 85 --db plasmidfinder > allinc.tab
    tail -n+2 allinc.tab >> allabr.tab
    """
}