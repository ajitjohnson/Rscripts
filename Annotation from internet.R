#load required packages

library(oligo)
library(hugene10sttranscriptcluster.db)
library(annotate)

#set wd to myworkingdirectory

setwd("myworkingdirectory")

#read in the raw data from the files and the pDatat

rawData <- read.celfiles(list.celfiles())

#rma normalization

rmaCore <- rma(rawData, target = 'core')

#annotation

ID <- featureNames(rmaCore)
Symbol <- getSYMBOL(ID, "hugene10sttranscriptcluster.db")
Name <- as.character(lookUp(ID, "hugene10sttranscriptcluster.db",
"GENENAME"))

#make a temporary data frame with all the identifiers...

tmpframe <-data.frame(ID=ID, Symbol=Symbol,
Name=Name,stringsAsFactors=F)
tmpframe[tmpframe=="NA"] <- NA

#assign data frame to rma-results

fData(rmaCore) <- tmpframe

#expression table with gene name and annotation info, processed with
sed after export to get the quotations in the right spot and remove NA
lines

write.table(cbind(pData(featureData(rmaCore))[,"Symbol"],exprs(rmaCore
)),file="better_annotation.csv", quote = FALSE, sep = ",")

----------