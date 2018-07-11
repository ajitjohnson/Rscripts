#Random Sampling
expression  <- read.delim("Z:/Ajit/ovariancancer/TCGA/2. remove dupilicate genes and samples/clinical.txt")
a = expression[,1]
b= sample(a, 100, replace = FALSE, prob = NULL)
write.table(b,file="test.txt", sep="\t")

#Ovarian cancer building the 9 gene signature
#Load data
library(ranger)
library(survival)


exp  <- read.delim("Z:/Ajit/ovariancancer/TCGA/3. training dataset/train_exp_1.txt")
cli  <- read.delim("Z:/Ajit/ovariancancer/TCGA/3. training dataset/train_cli.txt")

#Housekeeping
row.names(exp) = exp[,1]
exp = exp[,-1]
colnames(exp) <- gsub("\\.", "-", colnames(exp))
exp = data.frame(t(exp))

row.names(cli) = cli[,1]
cli = cli[,-1]

#Merge expression and clinical file
ordered = merge (exp, cli, by = "row.names")
row.names(ordered) = ordered[,1]
ordered = ordered[,-1]
ordered = subset(ordered, select = -c(Time, Event, age))

#Building model

survival_model <- ranger(group ~ .,
                         data = ordered,
                         seed = 1234,
                         importance = 'permutation',
                         verbose = TRUE,
                         num.trees = 10000)

#Picking important genes
important= data.frame(importance(survival_model))
important$dummy = "NA"
colnames(important) = c("Influence", "Dummy")
important= important[order(-important$Influence),]
dim(important)
head(important)

#Removing bottom 1/3rd of genes
remove = round((1/3)*(nrow(important)))
rownumber= nrow(important) - remove

imp = important[-(rownumber:nrow(important)),]
imp = rbind(imp, "group"= NA)


ordered <- ordered[row.names(imp)]
dim(ordered)

#Write out the top 9 genes
write.table(important,file="Round 2- top 7 genes.txt", sep="\t")

