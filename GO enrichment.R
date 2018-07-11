library(org.Hs.eg.db)
library(topGO)

y <- read.table("tissueclusters.txt", header=T, sep="\t")

x <- org.Hs.egSYMBOL
x <- toTable(x)
listEG= merge(x,y, by.x="symbol", by.y="gene")
 
GOenrich <- function(x) 
	{	
	universeEG <- Lkeys(org.Hs.egGENENAME)
	intgenes <- x
	geneList <- factor(as.integer(universeEG %in% intgenes))
	names(geneList)<- universeEG
	GOdata <- new("topGOdata", ontology ="BP", allGenes= geneList, nodeSize= 5, annot = annFUN.org, mapping = "org.Hs.eg.db", ID = "entrez")
	resultFisher <- runTest(GOdata, algorithm="classic", statistic ="fisher")
	allRes <- GenTable(GOdata, classic = resultFisher, ranksOf="classic", topNodes = 5)
	return(allRes)
	}

GOlists<- NULL
for(i in 1:nlevels(as.factor(listEG$clusters))) assign(paste("Cluster", i, sep="_"), GOenrich(listEG[listEG$clusters == i,2]))

cluster_annotation <- Cluster_1
for(i in 2:nlevels(as.factor(listEG$clusters)))	cluster_annotation=rbind(cluster_annotation, get(paste("Cluster", i, sep="_")))

cluster_annotation$cluster <- rep(1:nlevels(as.factor(listEG$clusters)), each=5)

write.table(cluster_annotation,file="cluster_annotated.txt", sep="\t")
