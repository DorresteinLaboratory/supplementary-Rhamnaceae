#######################################################
######   Biosynthetically informed Distances   ########
#######################################################
###############   Source Functions   ##################
#######################################################
############################# Robert R. Junker ########
#######################################################


library(vegan)
library(cluster)
library(GUniFrac)

BioSynDist <- function(biosynthesisInfo, SecMetComp){
  
  bioSynDistComp <- vegdist(biosynthesisInfo, method="bray", binary=T,na.rm = TRUE)
  clus_Comp <- hclust(bioSynDistComp, method="average") 
  unifracs <- GUniFrac(SecMetComp, as.phylo(clus_Comp), alpha=c(0, 0.5, 1))$unifracs   
  BSD <- unifracs[, , "d_1"] # Weighted UniFrac 
  #BioSynDist <- unifracs[, , "d_UW"] # Unweighted UniFrac
  #BioSynDist <- unifracs[, , "d_VAW"] # Variance adjusted weighted UniFrac
  #BioSynDist <- unifracs[, , "d_0"]      # GUniFrac with alpha 0
  #BioSynDist <- unifracs[, , "d_0.5"]    # GUniFrac with alpha 0.5
  
  r <- list("BSDist" = BSD,
            "dend_Comp" = clus_Comp)
  
  #return(r)
  invisible(r)
  
}

MergeDist <- function(BioSynDist, Dist, w = 0.878){
  
  BioSynDistStand <- as.dist(BioSynDist / max(BioSynDist))
  DistStand <- as.dist(Dist / max(Dist))
  
  mDist <- w * BioSynDistStand + (1 - w) * DistStand
  
  r <- list("mergedDist" = mDist)
  
  return(r)
  invisible(r)
}



