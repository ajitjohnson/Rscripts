# Step 1: Order the genes the based on the number of thimes they appear in the signatures
# Step 2: For each signature, rank the genes (1- whatever): Ranked genes
# Step 3: Melt the data in R
# Step 4: Create a table in excel, with total number of genes in each sig: sig number.txt


# load library ------------------------------------------------------------
library(reshape2)
library(circlize)
# load data ---------------------------------------------------------------
setwd("L:/Extended M drive/ImSIg/prev sig comparision")
sig=read.table("ranked genes1.txt", header=T, sep="\t")
number = read.table("sig number.txt", sep="\t")
number$V1 <- factor(number$V1, levels = c("Abbaset.al.2009",
                                          "Bechtet.al.",
                                          "ImSig",
                                          "Bindeaet.al.",
                                          "Newmanet.al.",
                                          "Angelovaet.al.",
                                          "Abbaset.al.2005",
                                          "Watkinset.al." ))
# Melt data ---------------------------------------------------------------

# sig_reshaped <- melt(sig, id = "genes")
# sig_reshaped = sig_reshaped[complete.cases(sig_reshaped), ]
# sig_reshaped = sig_reshaped[,-2]
# df = sig_reshaped[which(sig_reshaped$value!=""),]

# Looop melting

final <- data.frame()
sig_tmp <- sig[,-1]
for(i in 1:ncol(sig_tmp)){
  tmp = melt(sig_tmp, id =colnames(sig_tmp)[i])
  tmp =tmp[complete.cases(tmp), ]
  tmp$study = rep(names(tmp)[1],nrow(tmp))
  colnames(tmp) <- c("fromval","tostudy","toval","fromstudy")
  final <- rbind(final,tmp) 
  print(i)
}

# clear duplicated entries

final = final[duplicated(t(apply(final, 1, sort))),]


# Crap to plot the shit

circos.clear()

#pdf("circos plot1.pdf", width = 10, height = 10)
par(mar=c(1,1,1,1))

circos.par(cell.padding=c(0,0,0,0), 
           start.degree = -65, 
           gap.degree =2
)


circos.initialize(number[,1], x = number[,2])

# a = data.frame(c("#386cb0",
#       "#386cb0",
#       "#000000",
#       "#386cb0",
#       "#386cb0",
#       "#386cb0",
#       "#386cb0",
#       "#386cb0"), stringsAsFactors = F)

a = data.frame(c("#e41a1c", 
                 "#377eb8", 
                 "#000000", 
                 "#4daf4a", 
                 "#984ea3", 
                 "#ff7f00", 
                 "#ffff33", 
                 "#a65628"), stringsAsFactors = F)


colnames(a)[1] = "clr"

a$sig <- c("Abbaset.al.2009",
           "Bechtet.al.",
           "ImSig",
           "Bindeaet.al.",
           "Newmanet.al.",
           "Angelovaet.al.",
           "Abbaset.al.2005",
           "Watkinset.al." )

#Publications
circos.par(track.margin = c(0,0))

circos.trackPlotRegion(ylim = c(0, 1), track.height = 0.02,
                       bg.col = "white", 
                       bg.border = "black",
                       
                       panel.fun = function(x, y) 
                       {
                         #select details of current sector
                         name = get.cell.meta.data("sector.index")
                         i = get.cell.meta.data("sector.numeric.index")
                         xlim = get.cell.meta.data("xlim")
                         ylim = get.cell.meta.data("ylim")
                         #blank in part of main sector
                         circos.rect(xleft=xlim[1], ybottom=-0.1, xright=xlim[2], ytop=0.1, 
                                     col = "white", border = "white")
                         
                         #white line all the way around
                         #circos.rect(xleft=xlim[1], ybottom=0.3, xright=xlim[2], ytop=0.32, col = "white", border = "white")
                         
                       }
                       )

# #tracks
circos.trackPlotRegion(ylim = c(1, 2), track.height = 0.02,
                       bg.col = a$clr, bg.border = NA
)

for(i in 1:nrow(final)) 
{
    circos.link(final$fromstudy[i],final$fromval[i], 
                final$tostudy[i], final$toval[i], col = adjustcolor( a$clr[which(a$sig==final$fromstudy[i])], alpha.f = 0.5)) 

}


dev.off()
circos.clear()



         