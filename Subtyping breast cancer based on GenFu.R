#To identify molecular subtypes of breast cancer using the PAM50 model in genefu package

library(genefu)

#Choosing samples with clinical data and loading expression file
survival=read.table("cli.txt", header=T, row.names=1, sep="\t")
exp1=read.table("exp.txt", header=T, row.names=1, sep="\t")

#cli1= cli[,c("X_EVENT","X_OS")]
#cli2= cli1[complete.cases(cli1),]
#cli3= gsub("-",".", row.names(cli2))
#exp= exp1[,colnames(exp1) %in% cli3]

xxx= t(exp1)


#Loading Annotation
annot= data.frame(exp1$Gene.Symbol <- row.names(exp1))
names(annot)= "Gene.Symbol"

#Loading Annotation
forsubseting= t(exp1)

#Subtype Analysis
y= molecular.subtyping(sbt.model = "pam50", data = xxx, annot = annot, do.mapping = FALSE)

#Summary
split= table(y$subtype)
write.table(split,file="Table.txt", sep="\t")

xxx$PAM50 = y$subtype
LumA= names(y$subtype)[which(y$subtype == "LumA")]
LumB= names(y$subtype)[which(y$subtype == "LumB")]
Basal= names(y$subtype)[which(y$subtype == "Basal")]
Her2= names(y$subtype)[which(y$subtype == "Her2")]
Normal= names(y$subtype)[which(y$subtype == "Normal")]

LumA1= t(forsubseting[row.names(forsubseting) %in% LumA,])
Lumb1= t(forsubseting[row.names(forsubseting) %in% LumB,])
Basal1= t(forsubseting[row.names(forsubseting) %in% Basal,])
Her21= t(forsubseting[row.names(forsubseting) %in% Her2,])
Normal1= t(forsubseting[row.names(forsubseting) %in% Normal,])

#setwd("L:/Extended M drive/Pan Cancer Analysis/2. HZ Cutpoint analysis as per dominics script/3. Breast cancer/TCGA Breast subtype based on GeneFu/expression files")

write.table(LumA1,file="Exp_LumA.txt", sep="\t")
write.table(Lumb1,file="Exp_LumB.txt", sep="\t")
write.table(Basal1,file="Exp_Basal.txt", sep="\t")
write.table(Her21,file="Exp_Her2.txt", sep="\t")
write.table(Normal1,file="Exp_Normal.txt", sep="\t")

#setwd("L:/Extended M drive/Pan Cancer Analysis/2. HZ Cutpoint analysis as per dominics script/3. Breast cancer/TCGA Breast subtype based on GeneFu/clinical data")

#Subsetting Clinical data to be sent for bootstapping
#cli4= survival
#cli4$names= row.names(cli4)
#cli4= as.data.frame(lapply(cli4, function(y) gsub("-",".",y)))
#row.names(cli4)= cli4$names
#cli4= cli4[,c(2,1)]
#colnames(cli4) = c("t.dmfs", "e.dmfs")

cli= survival

LumA2= cli[row.names(cli) %in% LumA,]
Lumb2= cli[row.names(cli) %in% LumB,]
Basal2= cli[row.names(cli) %in% Basal,]
Her22= cli[row.names(cli) %in% Her2,]
Normal2= cli[row.names(cli) %in% Normal,]

write.table(LumA2,file="Cli_LumA.txt", sep="\t")
write.table(Lumb2,file="Cli_LumB.txt", sep="\t")
write.table(Basal2,file="Cli_Basal.txt", sep="\t")
write.table(Her22,file="Cli_Her2.txt", sep="\t")
write.table(Normal2,file="Cli_Normal.txt", sep="\t")

