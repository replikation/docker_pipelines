#!/usr/bin/env Rscript
# via docker run --rm -it -v $PWD:/input r-base
#install.packages("cowplot")
#install.packages("ggpubr")
#install.packages("viridis")

#library(ggpubr)
#library(viridis)


# my data
st <- 41.196e6L   # start
en <- 41.202e6L   # end
gtrack <- GenomeAxisTrack(cex = 1)  # set the font size larger
altrack <- AlignmentsTrack(
    "myseq.bam", isPaired = TRUE, col.mates = "deeppink"
)
plotTracks(
    list(gtrack, altrack, grtrack),
    from = st, to = en
)



svg("overview.svg", height = sizeh, width = sizew)
plot <-	
dev.off()


print("Done")

# Docs
# https://blog.liang2.tw/posts/2016/01/plot-seq-depth-gviz/


