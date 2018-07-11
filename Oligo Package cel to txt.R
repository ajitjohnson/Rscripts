library(oligo)
rawData = read.celfiles(as.character(list.celfiles()))
Eset = rma(rawData)
exp = exprs(Eset)

library("AnnotationDbi")
library("hgu133plus2.db") 
library("hugene20sttranscriptcluster.db")
library("hugene10sttranscriptcluster.db")
library("hgu133a.db")

annot=select(hgu133plus2.db, rownames(exp), c("SYMBOL"), "PROBEID")
annot=select(hugene10sttranscriptcluster.db, rownames(exp), c("SYMBOL"), "PROBEID")
annot=select(hgu133a.db, rownames(exp), c("SYMBOL"), "PROBEID")

comb=merge(exp,annot, by.y="PROBEID", by.x="row.names")

#Collapse to one gene

comb <- comb[,-1]
comb$SYMBOL <- as.factor(comb$SYMBOL)
class(comb$SYMBOL)
y=aggregate(.~SYMBOL, data=comb, max)
rownames(y)=y[,1]
y=y[,-1]
y=2^y
write.table(y,file="annotated.txt", sep="\t")


 pdf("Boxraw.pdf",width=6,height=6)
par(cex.lab=1.2)
par(cex.axis=0.5)
par(mar = c(15,5, 4, 2))
boxplot(rawData,las=2,main="Samples Expression Boxplot",ylab = "Expression", which="all")
dev.off()
