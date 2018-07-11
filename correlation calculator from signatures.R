#Average correlation calculator within a signature

exp=read.table("exp.txt", header=T, row.names=1, sep="\t")
sig = read.table("sig9-ciber.txt", header=TRUE, sep="\t")

sig.list <-
  lapply(unique(sig$cell), function(x)
    as.character(sig[sig$cell == x, 1]))

names(sig.list) <- unique(sig$cell)

#function to calculate average correlation

correlation = function(exp, geneset= sig.list){
  
  correlation_data = cor(t(exp[geneset,]),method="pearson")
  correlation_data [lower.tri(correlation_data,diag=TRUE)] <-NA
  m =as.data.frame(as.table(correlation_data))
  m1=subset(m,!is.na(m $Freq))
  m2=m1[,3]
  m2=summary(m2)
  correlation_final =m2["Median"]
  return(correlation_final)
}

l = lapply(sig.list, function(x) correlation(exp, geneset= x))
output = do.call(rbind.data.frame, l)
write.table(output,file="preliminaryimsig results.txt", sep="\t")
