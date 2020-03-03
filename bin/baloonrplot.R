#!/usr/bin/env Rscript
# via docker run --rm -it -v $PWD:/input r-base
#install.packages("cowplot")
#install.packages("ggpubr")
#install.packages("viridis")

library(ggpubr)
library(viridis)


# my data

input  <- read.delim("input.csv", row.names = 1, sep = ",")

sizew <- ceiling(( ncol(input) * 0.3 ) + 4 )
sizeh <- ceiling(( nrow(input) * 0.3 ) + 4 )
#sizew <- 8
#sizeh <- 8

#sorting data
data <- input[,sort(names(input))] 


svg("overview.svg", height = sizeh, width = sizew)
#png("phylum.png", height = 1000, width = 800, units = "px", pointsize = 12 )
plot <-	ggballoonplot(data, fill = "value", size.range = c(1, 10)) +
  	scale_fill_viridis_c(option = "inferno")
	ggpar(plot, legend.title = "Number of \n genes",
	xlab = "Samples", ylab = "Amount per gene")
dev.off()


print("Done")

# Docs
# https://www.rdocumentation.org/packages/ggpubr/versions/0.2.2/topics/ggballoonplot
# https://rpkgs.datanovia.com/ggpubr/reference/ggpar.html


