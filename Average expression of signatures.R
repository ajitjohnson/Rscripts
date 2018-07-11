#Average correlation calculator within a signature

exp=read.table("exp.txt", header=T, sep="\t")
sig = read.table("ImSIg feature selected genes.txt", header=TRUE,sep="\t")

sig.list <-
  lapply(unique(sig$cell), function(x)
    as.character(sig[sig$cell == x, 1]))

names(sig.list) <- unique(sig$cell)

#function to calculate average correlation

average = function(exp, geneset= sig.list){
  
  average_data = exp[exp$Symbols %in% geneset,]
  row.names(average_data) = average_data[,1]
  average_data = average_data[,-1]
  average_value = data.frame(colMeans(average_data, na.rm = TRUE))
  colnames(average_value) = names(geneset)
  return(average_value)
}

l = lapply(sig.list, function(x) average(exp, geneset= x))
output = do.call(cbind.data.frame, l)
write.table(output,file="Average expression of ImSig.txt", sep="\t")
