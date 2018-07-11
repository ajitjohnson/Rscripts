library(HiClimR)
library(data.table)

setwd("L:/Extended M drive/ImSIg/Single cell validation/melanoama/intre")
dat <- fread("macrophages_nomagic.csv", header = TRUE)
setDF(dat)
rownames(dat) <- paste(rownames(dat),dat$GENE)
dat <- dat[,-1]
#dat = 2^dat


# Remove bottom endf ------------------------------------------------------

Keep_genes = rowSums(dat > 0) > 5
dat = dat[Keep_genes,]

# Correlation -------------------------------------------------------------
dat_cor <- fastCor((dat))
#dat_cor = as.matrix(dist(t(dat)))
#dat_cor = cor(dat, method = "spearman")
diag(dat_cor) <- 0

# Weights -----------------------------------------------------------------
dat_corp <- dat_cor
dat_corp[dat_corp<0] = 0
con_mat <- data.frame(apply(dat_cor,1,max))
colnames(con_mat)[1] <- "max"
con_mat$wts <- colSums(dat_cor, na.rm = T)

# Keep edges --------------------------------------------------------------
mn =10
mx =10
con_mat$keep <- round((mx-mn)*(( con_mat$wts - min(con_mat$wts) )/ (max(con_mat$wts) - min(con_mat$wts))) + mn)
sub_cor <- dat_cor
sub_cor[!is.na(sub_cor)] <- NA
for(i in 1:nrow(dat_cor)){
  tmp <- sort(dat_cor[i,], decreasing = T)
  tmp <- names(tmp[1:con_mat$keep[i]])
  sub_cor[i,tmp] <- dat_cor[i,tmp]
}
tsub_cor <- t(sub_cor)
sub_cor[is.na(sub_cor)] = tsub_cor[is.na(sub_cor)]

sub_cor[sub_cor==0] <- NA
clusters <- data.frame()
joints <- data.frame()
count <- 1
while (ncol(clusters)!=1) {
  clusters <- data.frame(matrix(NA, nrow = nrow(sub_cor))) # All clusters
  nd_clus <- rep(FALSE,nrow(sub_cor)) # checking all nodes for clusters
  nd_clus[which(apply(sub_cor,1,function(x) all(is.na(x))))] <- TRUE # nodes unconnected are given their own cluster
  # Find clusters
  while(any(nd_clus==FALSE)){
    cluster_nd_found <- FALSE # Cluster found(True) or not(False) in iteration
    clus_tmp <- rep(FALSE, nrow(sub_cor)) # representes the cluster being identified, True values for nodes in cluster
    clus_index <- min(which(nd_clus == FALSE)) # find next node cluster
    clus_tmp[clus_index] <- TRUE # initialize search
    while(cluster_nd_found == FALSE){
      prev_clus <- clus_tmp
      for(j in which(clus_tmp)){ #
        clus_tmp[which(!is.na(sub_cor[j,]))] <- TRUE
        nd_clus[which(!is.na(sub_cor[j,]))] <- TRUE
      }
      if(all(prev_clus == clus_tmp)){ # if all nodes of clusters found
        cluster_nd_found <- TRUE
      }
    }
    clusters <- cbind(clusters,clus_tmp)
  }
  
  clusters <- data.frame(clusters[,-1])
  rownames(clusters) <- colnames(sub_cor)
  clus_tmp <- as.numeric()
  
  if(ncol(clusters)!=1){
    # Join clusters
    for(i in 1:ncol(clusters)){
      for(j in 1:ncol(clusters)){
        if(j==i){
          clus_tmp[j] <- NA
        }
        else{
          ed_wt <- max(dat_cor[which(clusters[,i]),which(clusters[,j])], na.rm = T)
          clus_tmp[j] <- ed_wt
        }
      }
      closest_clus <- which.max(clus_tmp)
      tmp <- dat_cor[which(clusters[,i]),which(clusters[,closest_clus])]
      ed_wt <- max(dat_cor[which(clusters[,i]),which(clusters[,closest_clus])], na.rm = T)
      ed <- which(tmp==ed_wt,arr.ind = T)
      x <- rownames(tmp)[ed[1,1]]
      y <- colnames(tmp)[ed[1,2]]
      sub_cor[x,y] <- ed_wt
      sub_cor[y,x] <- ed_wt
      
      # Save joints
      joints[count,1] <- x
      joints[count,2] <- y
      joints[count,3] <- ed_wt
      
      count <- count+1
    }
  }
  x <- paste("Graph components: ",ncol(clusters))
  print(x)
  #cat("\r",paste("Graph components: ",ncol(clusters)))
}


# Plot KNN ----------------------------------------------------------------

sub_cor[is.na(sub_cor)]=0
library(igraph)

g = graph_from_adjacency_matrix(sub_cor, mode = "undirected", weighted = T)
plot(g, vertex.label = NA, vertex.size = 1, layout = layout_with_fr)

# Save for graphia --------------------------------------------------------
sub_cor[sub_cor==0] = NA
sum(!is.na(sub_cor))
k <- 1
out_cor <- data.frame()
for(i in 1:(nrow(sub_cor)-1)){
  for(j in (i+1):ncol(sub_cor)){
    if(!is.na(sub_cor[i,j])){
      out_cor[k,1] <- rownames(sub_cor)[i]
      out_cor[k,2] <- colnames(sub_cor)[j]
      out_cor[k,3] <- sub_cor[i,j]
      k <- k+1
    }
  }
  cat("\r", paste(round(i/(nrow(sub_cor)-1)*100,2), "% Complete"))
}

#setwd("M:/My work/Projects/Network pruning/weight dist keep")
nm <- paste("monocyte_KNN_alog.txt")
write.table(out_cor, nm, sep = "\t", row.names = F, col.names = F)
