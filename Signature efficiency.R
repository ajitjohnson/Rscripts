# Signature efficiency based on median correaltion
# Import Lib's ------------------------------------------------------------
library(plyr)
library(ggplot2)
# Import dataset and signature -------------------------------------------
setwd("L:/Extended M drive/ImSIg/ImSig Part 2/8. Validation datasets/3- Trachoma dataset")

exp=read.table("exp.txt", header=T, row.names=1, sep="\t")
sig = read.table("tissue-sig.txt", header=TRUE, sep="\t") ## Add signature to analyse
sig.list <-
  lapply(unique(sig$cell), function(x)
    as.character(sig[sig$cell == x, 1]))
names(sig.list) <- unique(sig$cell)


#function to calculate median correlation for each gene
correlation = function(exp, geneset= sig.list){
  correlation_data = cor(t(exp[geneset,]),method="pearson")
  correlation_data [lower.tri(correlation_data,diag=TRUE)] <-NA
  m =as.data.frame(as.table(correlation_data))
  m1=subset(m,!is.na(m $Freq))
  Median_value<-data.frame(with(m1, tapply(m1$Freq, m1$Var1, median), na.rm = TRUE))
  return(Median_value)
}

l = lapply(sig.list, function(x) correlation(exp, geneset= x))
output = do.call(rbind.data.frame, l)
output$cell_name<-sapply(strsplit(rownames(output), "\\."), '[',1)
colnames(output) = c("correlation", "celltype")
#write.table(output,file="preliminaryimsig results.txt", sep="\t")


# Overall median correlation ----------------------------------------------

o_correlation = function(exp, geneset= sig.list){
  
  correlation_data = cor(t(exp[geneset,]),method="pearson")
  correlation_data [lower.tri(correlation_data,diag=TRUE)] <-NA
  m =as.data.frame(as.table(correlation_data))
  m1=subset(m,!is.na(m $Freq))
  m2=m1[,3]
  m2=summary(m2)
  correlation_final =m2["Median"]
  return(correlation_final)
}

M = lapply(sig.list, function(x) o_correlation(exp, geneset= x))
output_meds = do.call(rbind.data.frame, M)
output_meds$celltype = row.names(output_meds)
colnames(output_meds) = c("correlation", "celltype")
#write.table(output_meds,file="tissuesig_median_correlation.txt", sep="\t")

# Plot the correlation values ---------------------------------------------
output$celltype = as.factor(output$celltype)

plot_obj <- ggplot(data=output, aes(x=celltype, y=correlation))+ theme_classic()+
  geom_boxplot(outlier.shape=NA) + #avoid plotting outliers twice
  geom_jitter(position=position_jitter(width=.1, height=0))+
  theme(axis.title.x = element_text(face="bold", size=12),
  axis.text.x  = element_text(angle=90, size=12, colour = "black", face="bold", hjust=0.95,vjust=0.2),
  axis.text.y  = element_text(size=12, colour = "black", face="bold"))+
  labs(y="Correlation", x="Signature")+ scale_y_continuous(limits = c(-0.5, 1.2))+
  geom_hline (yintercept= c(0, 0.5), colour = "darkgrey")

plot_obj + geom_text(data=output_meds, aes(x=celltype, y=1.2, label = round(correlation, digits = 2)))

