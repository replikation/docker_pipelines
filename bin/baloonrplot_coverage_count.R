#!/usr/bin/env Rscript

#table structure like this:

	#type,reactor1,reactor2,reactor3,reactor4,type
	#classA,1,2,3,4,genome
	#classB,NA,2,6,2,genome
	#classC,1,1,1,1,plasmid
	#inc,NA,NA,NA,7,genome

	#type,reactor1,reactor2,reactor3,reactor4,type
	#classA,10,25,33.4,22,genome
	#classB,NA,4,16,22,genome
	#classC,1,10,11,5,plasmid
	#inc,NA,NA,10,70,genome


# to install
#install.packages("reshape2")
#install.packages("viridis")
#install.packages("ggplot2")

#libs
library(reshape2)
library(ggplot2)
library(viridis)

dt <- read.csv(file = 'input.csv')
dt2 <- read.csv(file = 'input_coverage.csv')

dm1 <- melt(dt)
dm2 <- melt(dt2)  

names(dm1)[names(dm1) == "value"] <- "Genes"

balloon_melted <- merge(dm1, dm2)

sizew <- ceiling(( ncol(dt) * 0.4 ) + 4 )
sizeh <- ceiling(( nrow(dt) * 0.1 ) + 4 )

options(ggplot2.continuous.colour="viridis")

plot <- ggplot(balloon_melted, aes(x = variable, y = type)) +
	facet_grid(method ~ type.1, scales = "free", space = "free") +
	theme_light() +
	theme(panel.border = element_blank()) + 
	labs(y = element_blank(), x="Sample") +
	geom_point( aes(size=value, colour=Genes) ) +
	scale_size_area(max_size=15) +
	labs(size="Abundance", colour="Gene count") + 
	theme(axis.text.x = element_text(angle = 45, hjust = 1))


svg("overview.svg", height = sizeh, width = sizew)
print(plot)
dev.off()

pdf("overview.pdf", height = sizeh, width = sizew)
print(plot)
dev.off()