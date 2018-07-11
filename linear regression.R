#Linear regression

library(ggplot2)

setwd("L:\Extended M drive\ImSIg\Single cell validation\melanoama\patient level validation\linear regression")

dat = read.table("P9-r.txt", header=T, row.names= 1, sep="\t")

linearMod <- lm(actual ~ predicted, data=dat)
linearMod
summary(linearMod)

modelSummary <- summary(linearMod)  # capture model summary as an object
modelCoeffs <- modelSummary$coefficients  # model coefficients
beta.estimate <- modelCoeffs["predicted", "Estimate"]  # get beta estimate for speed
std.error <- modelCoeffs["predicted", "Std. Error"]  # get std.error for speed
t_value <- beta.estimate/std.error  # calc t statistic
p_value <- 2*pt(-abs(t_value), df=nrow(dat)-ncol(dat))  # calc p Value
f_statistic <- linearMod$fstatistic[1]  # fstatistic
f <- summary(linearMod)$fstatistic  # parameters for model p-value calc
f <- summary(linearMod)
model_p <- pf(f[1], f[2], f[3], lower=FALSE)




p_val = paste ("P", round(p_value,3), sep = " = ", collapse = NULL)
r_val = paste ("R", round(f$r.squared, 3), sep = " = ", collapse = NULL)
corr = paste ("r", round(cor(dat$actual,dat$predicted),3), sep = " = ", collapse = NULL)

#Scatter plot
#scatter.smooth(x= dat$actual, y=dat$predicted, main="P1")  # scatterplot


ggplot(dat, aes(x=actual, y=predicted)) + 
  geom_point(aes(fill = row.names(dat)),colour = "black", size = 20, stroke = 0.2, shape = 21) + theme_bw() + 
  geom_smooth(method=lm, se=FALSE, size = 4)+
  theme(plot.title = element_text(face="bold", size=25, hjust = 0.5),
        axis.text.x  = element_text(angle=90, size=30, colour = "black", face="bold", hjust=0.95,vjust=0.2),
        axis.text.y  = element_text(size=30, colour = "black", face="bold"),
        legend.position="none")+
  labs(y="", x="")+ scale_y_continuous(limits = c(0, 5.5))+
  annotate("text", x = 42, y = 5, label = p_val, size = 20, fontface =2) +
  #annotate("text", x = 20, y = 4, label = r_val, size = 20, fontface =2) +
  annotate("text", x = 42, y = 4, label = corr, size = 20, fontface =2)


# Dummy data 
dt <- read.csv("mel_magic_zero_removed.csv", header = T)
genes = dt[,1]
bcell = data.frame(sample(genes, 32, replace = FALSE))
mac = data.frame(sample(genes, 74, replace = FALSE))
tcell = data.frame(sample(genes, 44, replace = FALSE))
nkcell = data.frame(sample(genes, 20, replace = FALSE))

write.table(nkcell, file = "bcell.txt", sep= "\t")
