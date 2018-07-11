
train  <- read.delim("Z:/Ajit/ovariancancer/TCGA/3. training dataset/train_exp_svm.txt")
test  <- read.delim("Z:/Ajit/ovariancancer/TCGA/4. test dataset/test_exp_svm.txt")
testing = subset( test, select= -group)
sig.list  <- read.delim("Z:/Ajit/ovariancancer/TCGA/3. training dataset/sig list all combinations.txt")
sig.list = levels(sig.list[1,])

#Housekeeping
row.names(train) = train[,1]
train = train[,-1]
row.names(testing) = testing[,1]
testing = testing[,-1]
row.names(test) = test[,1]
test = test[,-1]

library(e1071)
tuned = tune(svm, group~. , data= train, kernel = "polynomial", ranges = list(cost = c(0.001, 0.01, .1, 1, 10, 100)))
summary (tuned)

output = list()
for (i in 1:length(sig.list)){
formula <- formula(paste('group ~ ', sig.list[i]))
svmfit = svm(formula , data= train, kernel = "polynomial", cost = 0.001, scale = FALSE)
#print(svmfit)
#plot(svmfit, train[,col])
p = predict (svmfit, testing, type ="class")
#table (p, test[,9])
output = c(output, mean (p == test[,9]))
}

final = cbind(data.frame(sig.list), data.frame(unlist(output)))
