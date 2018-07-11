#Random Forest
library(randomForest)
data <- read.delim("Z:/Ajit/esophagealcancer/TCGA/allsamples/Random forest/top 10 genes.txt", row.names = 1)
fit <- randomForest(Feature ~ ., data, ntree=500)
important= data.frame(importance(fit))
important$dummy = "NA"
colnames(important) = c("Influence", "Dummy")
important= important[order(-important$Influence),]
