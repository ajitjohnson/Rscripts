# Diffusion map
library(destiny)
library(Rtsne)
library(ggplot2)
#library(data.table)
#option(stringsAsFactors = F)

setwd("L:/Extended M drive/ImSIg/Single cell validation/HNSC/imsig")
setwd("L:/Extended M drive/ImSIg/Single cell validation/melanoama/imsig validation")

#dt <- read.csv("hnsc_magic_zero_undefined_cells_removed.csv", header = T)
dt <- read.csv("only imsig immune_non-malignant tsne.csv", header = T)
dt <- dt[which(!duplicated(dt$X)),]
rownames(dt) <- dt$X
dt <-  data.frame(dt[,-1])
exp =  data.frame(t(dt), stringsAsFactors = F)

# Convert into numeric and factor
exp$cell <- as.factor(exp$cell)
exp[,2:ncol(exp)] = data.frame(lapply(exp[,2:ncol(exp)], as.numeric))

dat <- as.ExpressionSet(exp)
dm <- DiffusionMap(dat)
plot(dm, col = exp$cells)

#Test
head(exp1[,1:5])
class(exp1$cells)
class(exp$A1BG)

#Only immune cells
exp1 <- exp[which(exp$cells == "T cell" | 
                    exp$cells == "Mast" |
                    exp$cells == "Dendritic" |
                    exp$cells == "Macrophage" |
                    exp$cells == "B cell" ), ]
exp1$cells  <- factor(exp1$cells)

#Only immune melanoma
exp1 <- exp[which(exp$cell == "B" | 
                    exp$cell == "Macrophage" |
                    exp$cell == "NK" |
                    exp$cell == "T" ), ]
exp1$cell  <- factor(exp1$cell)


#Diffusion Plot
dat <- as.ExpressionSet(exp1)
dm <- DiffusionMap(dat)
plot(dm, 1:2, pch = 20, col = exp1$cells)

#TSNE
tsne <- Rtsne(as.matrix(exp1[,2:ncol(exp1)]), perplexity = 100)
plot(tsne$Y, col = factor(exp1$cell))

#Plot in GGplots

p = data.frame(tsne$Y)
p$cells = exp1$cell

plot = ggplot(p, aes(X1, X2))
plot + geom_point(aes( fill = factor(cells)),colour = "black", size = 2, stroke = 0.2, shape = 21) + theme_classic()+
        theme(plot.title = element_text(face="bold", size=25, hjust = 0.5),
        axis.text.x  = element_text(angle=90, size=12, colour = "black", face="bold", hjust=0.95,vjust=0.2),
        axis.text.y  = element_text(size=12, colour = "black", face="bold"))
