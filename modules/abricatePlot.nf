process abricatePlot {
    publishDir "${params.output}/${name}/ABR-Screening", mode: 'copy', pattern: "*.pdf"
    label 'ggplot2'
  input:
    set val(name), file(results) 
  output:
	  set val(name), file("*.pdf") 
  script:
    """
    #!/usr/bin/Rscript

    library(ggplot2)
    
    inputdata <- read.table("${results}", header = TRUE, sep = ";")
    
    pdf("classification-overview.pdf", height = 5, width = 15)
      ggplot(data=inputdata, aes(x=type, y=amount)) +
      geom_bar(stat="identity")
    dev.off()
    """
}

// tempDir may be neccessary

/*
    #install.packages("ggplot2")
    #library(ggplot2)
svg("classification-overview.svg", height = 5, width = 15)
ggplot(data=inputdata, aes(x=database.taxonomic.level, y=amount.reads, fill=classification)) + 
  geom_bar(position = position_fill(reverse = TRUE),stat="identity") +
  facet_grid(~reactor) + 
  labs(x="", y="Proportion of classified reads", fill="status") + 
  theme(plot.title = element_text(size=25, margin=margin(t=20, b=20))) +
  theme(axis.text.x = element_text(angle = 90))  +
  theme(panel.background = element_blank(), panel.grid.major = element_blank()) +
  scale_fill_manual(values = c("#91AEC1", "#423E37")) +
  theme(legend.title = element_blank())
dev.off()
*/