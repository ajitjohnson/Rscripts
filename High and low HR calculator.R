expression  <- read.delim("Z:/Ajit/ovariancancer/TCGA/2. remove dupilicate genes and samples/expression_trimmed.txt")
survival  <- read.delim("Z:/Ajit/ovariancancer/TCGA/2. remove dupilicate genes and samples/clinical.txt")
tissue.sig <- read.delim("Z:/Ajit/esophagealcancer/TCGA/allsamples/gene-sig.txt")

#Preprocessing for this dataset
#Expression file
row.names(expression) = expression[,1]
expression = expression[,-1]
colnames(expression) <- gsub("\\.", "-", colnames(expression))
#Survival file
row.names(survival)= survival[,1]
survival= survival[,-1]

## Subset both expression and survival objects by samples which are common to both
shared.ids <- intersect(colnames(expression), row.names(survival))
surv <- survival[shared.ids,]
colnames(surv) <- c("Time", "Event", "age")
expr.base <- expression[shared.ids]
## Make sure sample ids are common and in the same ordering
expr <- expr.base[row.names(surv)]
#Remove genes that are not expressed
expr = expr[!(apply(expr, 1, function(y) all( y == 0))),]

#Sig.list
#sig.list <- tissue.sig[,1]
#sig.list = levels(sig.list)

sig.list = row.names(expr)

HR = function(expr, surv, geneset= sig.list){
  #Load Lib
  library("survival")
  #Frame data
  data.base <- data.frame(colMeans(expr[geneset,],  na.rm = TRUE))
  surv.base <- subset(surv, select = c(Time, Event))
  data_surv <- merge(data.base, surv.base, by = "row.names")
  row.names(data_surv) = data_surv[,1]
  data_surv = data_surv[,-1]
  colnames(data_surv) = c("Mean", "Time", "Event")
  data_surv$Time = as.numeric(data_surv$Time)
  data_surv$Event = as.numeric(data_surv$Event)
  data_surv$Mean = as.numeric(data_surv$Mean)
  #Order data and split into chunks
  data_surv <- data_surv[order(data_surv$Mean),]
  #Median
  data_surv$group <- ifelse(data_surv$Mean<= median(data_surv$Mean), "low", "high")
  #Three groups
  #n = 3
  #split(data_surv, nrow(data_surv) )
  #Survival analysis
  fit = coxph(Surv(Time, Event) ~ group, data_surv)
  HR = summary(fit)
  Hazard_Ratio = data.frame(HR["coefficients"])
  row.names(Hazard_Ratio)= geneset[1]
  colnames(Hazard_Ratio)= c("coef", "exp.coef(Hazard L/H)", "se.coef", "z", "Hazard P_Value")
  return(Hazard_Ratio)
}



l = lapply(sig.list, function(x) HR(expr, surv, geneset= x))
output = do.call(rbind.data.frame, l)
write.table(output,file="Median HR.txt", sep="\t")
