
train  <- read.delim("Z:/Ajit/ovariancancer/GSE9891/train/train.txt")
test  <- read.delim("Z:/Ajit/ovariancancer/GSE9891/test/test.txt")
testing = subset( test, select= -group)
sig.list  <- read.delim("Z:/Ajit/ovariancancer/TCGA/3. training dataset/sig list all combinations- Round 3.txt")
sig.list = levels(sig.list[1,])
sig.list = "DPY19L1P1+DNM1P35+LINC00664"

library(ranger)

#Housekeeping
row.names(train) = train[,1]
train = train[,-1]
row.names(testing) = testing[,1]
testing = testing[,-1]
row.names(test) = test[,1]
test = test[,-1]

output = list()
for (i in 1:length(sig.list)){
  formula <- formula(paste('group ~ ', sig.list[i]))
  
  survival_model <- ranger(formula,
                           data = train,
                           seed=1234,
                           verbose = TRUE,
                           num.trees = 10000,
                           write.forest=TRUE )
  
  
  p <- predict(survival_model, testing)
  output = c(output, mean(data.frame(p[3]) == data.frame(test[,4])))
  
}

final = cbind(data.frame(sig.list), data.frame(unlist(output)))
write.table(final,file="Specificity HR.txt", sep="\t")
