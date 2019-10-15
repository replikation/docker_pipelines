process abricatePlot {
    publishDir "${params.output}/${name}/ABR-Screening", mode: 'copy', pattern: "*.pdf"
    label 'ggplot2'
    echo true
  input:
    set val(name), file(results) 
  output:
	  set val(name), file("*.pdf") 
  script:
    """
    #!/usr/bin/Rscript

    library(ggplot2)
    
    inputdata <- read.table("${results}", header = TRUE, sep = ";")
    
    pdf("classification-overview.pdf", height = 6, width = 10)
      ggplot(data=inputdata, aes(x=type, y=amount, fill=group)) +
      geom_bar(stat="identity") +
      theme(legend.position = "none") +
      facet_wrap(~ group, scales = "free_y") + coord_flip()
    dev.off()
    """
}

// tempDir may be neccessary
