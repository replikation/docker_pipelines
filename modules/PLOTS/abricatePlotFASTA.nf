process abricatePlotFASTA {
    publishDir "${params.output}/${name}/ABR-Screening", mode: 'copy', pattern: "*.pdf"
    label 'ggplot2'
    errorStrategy { task.exitStatus in 1..1 ? 'ignore' : 'terminate'}
  input:
    tuple val(name), val(method), file(results) 
  output:
	  tuple val(name), val(method), file("*.pdf") optional true
  script:
    """
    #!/usr/bin/Rscript

    library(ggplot2)
    
    inputdata <- read.table("${results}", header = TRUE, sep = ";")
    
    pdf("classification-${method}.pdf", height = 6, width = 10)
      ggplot(data=inputdata, aes(x=type, y=amount, fill=group)) +
      geom_bar(stat="identity") +
      theme(legend.position = "none") +
      facet_wrap(~ group, scales = "free_y") + coord_flip()
    dev.off()
    """
}

// tempDir may be neccessary
