# load the libraries
library(GenomicRanges)
library(rtracklayer)
library(ggplot2)
library(lattice)
library(reshape)
library(reshape2)
library(latticeExtra)
library(DESeq)
library(gplots)
library(grDevices)
library(BSgenome.Hsapiens.UCSC.hg19)
library(ShortRead)
library(Rsamtools)
library(biomaRt)
library(GMD)
library(Vennerable)
library(matrixStats)
library(data.table)
library(directlabels)
library(grid)

# to forget about stringAsFactors problem
options(stringsAsFactors = FALSE)

### Functions
# The arrange function below puts more than one plot in the same image window.
multiplot <-function(..., plotlist=NULL, cols){
    require(grid)# Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)

    numPlots = length(plots)# Make the panel
    plotCols = cols                          # Number of columns of plots
    plotRows = ceiling(numPlots/plotCols)# Number of rows needed, calculated from # of cols# Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(plotRows, plotCols)))
    vplayout <-function(x, y)
        viewport(layout.pos.row = x, layout.pos.col = y)
    # Make each plot, in the correct location
    for(i in 1:numPlots){
        curRow = ceiling(i/plotCols)
        curCol =(i-1)%% plotCols +1
        print(plots[[i]], vp = vplayout(curRow, curCol ))
    }
}
