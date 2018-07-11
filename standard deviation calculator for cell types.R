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
comb <- comb[!duplicated(comb$Row.names),]
#Collapse to one gene
rownames(comb) <- comb$Row.names
comb <- comb[,-1]
comb$SYMBOL <- as.factor(comb$SYMBOL)
class(comb$SYMBOL)
y=aggregate(.~SYMBOL, data=comb, max)
rownames(y)=y[,1]
y=y[,-1]
write.table(y,file="log-annotated.txt", sep="\t")
y=read.table("log-annotated.txt", header=T, row.names=1, sep="\t")

setwd("C:/Users/v1anirma/Desktop/algorithm formulation/Datasets") 
pdata=read.csv("pdata_blood.csv", header=T)


macrophages <- y[rownames(y) %in% as.character(pdata[pdata$cell == "Macrophages/Dendritic cells", 1]),]
monocytes <- y[rownames(y) %in% as.character(pdata[pdata$cell == "Monocytes", 1]),]
bcell <- y[rownames(y) %in% as.character(pdata[pdata$cell == "B cells", 1]),]
neutrophils <- y[rownames(y) %in% as.character(pdata[pdata$cell == "Neutrophils", 1]),]
nkcells <- y[rownames(y) %in% as.character(pdata[pdata$cell == "NK cells", 1]),]
plasmacells <- y[rownames(y) %in% as.character(pdata[pdata$cell == "Plasma cells", 1]),]
platelets <- y[rownames(y) %in% as.character(pdata[pdata$cell == "Platelets", 1]),]
Tcells <- y[rownames(y) %in% as.character(pdata[pdata$cell == "T cells", 1]),]
Ribosomes <- y[rownames(y) %in% as.character(pdata[pdata$cell == "Ribosomes", 1]),]
IFN <- y[rownames(y) %in% as.character(pdata[pdata$cell == "IFN", 1]),]
CellCycle <- y[rownames(y) %in% as.character(pdata[pdata$cell == "Cell Cycle", 1]),]

macrophagesE=sd(colMeans(macrophages))
monocytesE=sd(colMeans(monocytes))
bcellE=sd(colMeans(bcell))
neutrophilsE=sd(colMeans(neutrophils))
nkcellsE=sd(colMeans(nkcells))
plasmacellsE=sd(colMeans(plasmacells))
plateletsE=sd(colMeans(platelets))
TcellsE=sd(colMeans(Tcells))
RibosomesE=sd(colMeans(Ribosomes))
IFNE=sd(colMeans(IFN))
CellCycleE=sd(colMeans(CellCycle))
expression= rbind(macrophagesE,monocytesE,bcellE,neutrophilsE,nkcellsE,plasmacellsE,plateletsE,TcellsE,RibosomesE,IFNE,CellCycleE)
write.table(expression,file="standard deviation.txt", sep="\t")
