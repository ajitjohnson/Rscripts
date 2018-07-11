
library(limma)
exprs <- read.table("normalised_expression.txt", header=TRUE, row.names=1)
design <- read.table("contrast.matrix.txt", header=TRUE)
fit <- lmFit(exprs,design)
contrast.matrix <- makeContrasts(KO-WT, levels=design)
fit2 <- contrasts.fit(fit,contrast.matrix)
fit2 <- eBayes(fit2)
topTable(fit2, number=20, sort.by="p")
fullList <- topTable(fit2, number=nrow(fit2), sort.by="logFC", resort.by="logFC")
write.table(fullList, file="full_results.txt", row.names=TRUE, sep="\t", quote=FALSE)

