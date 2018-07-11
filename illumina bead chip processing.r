
library(lumi)
fileName <- 'GSE49454_non-normalized.txt'

x.lumi <- lumiR(fileName)
z= lumiN(x.lumi)
write.exprs(z, file="normalised.txt")

library("AnnotationDbi")
library("illuminaHumanv4.db")

x=read.table("normalised.txt", header=T, sep="\t")
y=as.character(x[,1])
annot=select(illuminaHumanv4.db, y, c("SYMBOL"), "PROBEID")
comb=merge(x,annot, by.y="PROBEID", by.x="X")
comb <- comb[,-1]
comb$SYMBOL <- as.factor(comb$SYMBOL)
y=aggregate(.~SYMBOL, data=comb, max)
rownames(y)=y[,1]
y=y[,-1]
write.table(y,file="annotated.txt", sep="\t")


#OR

library(GEOquery)
eset = getGEO('GSE49454')[[1]]
z= lumiN(eset)