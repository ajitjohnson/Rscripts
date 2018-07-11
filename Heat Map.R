#Heat Map

library(gplots)
library("RColorBrewer")

x=read.table("collapsed correlation.txt", header=T, row.names=1, sep="\t")

x=as.matrix(x)



x=round(x, digits = 2)


colors = c(seq(0,0.5,length=100),seq(0.51,0.64,length=100),seq(0.65,1,length=100))
 my_palette2 <- colorRampPalette(c("#ffffff", "#feb24c", "#f03b20"))(n = 299)
heatmap.2(x, Rowv=FALSE, Colv=FALSE, cellnote=x, notecol="black", notecex=3.8, dendrogram = "none", na.color="black", col=my_palette2, breaks= colors, symm=F,symkey=F,symbreaks=T, scale="none", offsetRow = 0.01, offsetCol = 0.05, margins = c(5,15), trace = "none", cexCol = 1.8)

# colors = c(seq(3,4,length=100), seq(5,6,length=100))
# my_palette2 <- colorRampPalette(c("#ef8a62"))(n = 299)
# heatmap.2(log(x), Rowv=TRUE, Colv=FALSE, dendrogram = "none", na.color="blue",col=my_palette2, breaks= colors, symm=F,symkey=F,symbreaks=T, scale="none", offsetRow = 0.01, offsetCol = 0.05, margins = c(5,15), trace = "none", cexCol = 0.60)