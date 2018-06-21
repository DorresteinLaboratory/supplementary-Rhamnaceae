##########################################################################################
#                                                                                        #
#   Calculate consensus chemical classes per molecular families using ClassyFire         #
#                                                                                        #
##########################################################################################

sc <- read.table("c1e1c4beb4f6426b9c9eeed05f4e7cf7.tsv",sep="\t",header = T,comment.char = "",stringsAsFactors = F,quote="") # downloaded from https://gnps.ucsd.edu/ProteoSAFe/status.jsp?task=e9e02c0ba3db473a9b1ddd36da72859b (Download Cytoscape data, clusterinfosummarygroup_attributes_withIDs_withcomponentID)

# subdivide molecular family
subclusters <- read.csv("Subclusters.csv",stringsAsFactors = F,header=F)
sc$componentindex[which(sc$cluster.index %in% subclusters$V1[which(subclusters$V2=="Subcluster2")])] <- "Subcluster2"
# subdivide molecular family

ci <- unique(sc$componentindex)
ci <- ci[-which(ci==-1)] # remove single nodes

# downloaded from https://proteomics2.ucsd.edu/ProteoSAFe/status.jsp?task=6b515b235e0e4c76ba539524c8b4c6d8
naph <- read.table("node_attributes_table.tsv",sep="\t",quote = "",comment.char = "",header=T)
naph <- naph[,c(which(colnames(naph) %in% c("cluster.index","Smiles","MetFragSMILES", "ConsensusSMILES", "FusionSMILES")))]

clin <- c()
sm <- c()
for (i in 1:length(ci)){
  
  sc_smiles <- sc[,c(which(colnames(sc)== "cluster.index"), which(colnames(sc)== "componentindex"))]
  sc_smiles <- sc_smiles[sc_smiles$componentindex==ci[i],]
  
  NAPH <- naph[which(naph$cluster.index %in% sc_smiles$cluster.index),]
  
  matches <- length(unique(c(which(NAPH$Smiles!=""),which(NAPH$MetFragSMILES!=""),which(NAPH$FusionSMILES!=""),which(NAPH$ConsensusSMILES!=""))))
  matches_per <- matches/nrow(sc_smiles)
  
  HGNPS <- unname(unlist(sapply(as.character(NAPH$Smiles),strsplit,split=",")))
  HMetFrag <- unname(unlist(sapply(as.character(NAPH$MetFragSMILES),strsplit,split=",")))
  HFusion <- unname(unlist(sapply(as.character(NAPH$FusionSMILES),strsplit,split=",")))
  HConsensus <- unname(unlist(sapply(as.character(NAPH$ConsensusSMILES),strsplit,split=",")))
  
  allsmiles <- unique(c(HGNPS,HMetFrag,HFusion,HConsensus))
  if (length(allsmiles[which(is.na(allsmiles))]) !=0){
    allsmiles <- allsmiles[-which(is.na(allsmiles))]
  }
  
  if(length(allsmiles) != 0){
    sm <- c(sm,allsmiles)
    ident <- paste(as.character(rep(ci[i],length(allsmiles))),round(matches_per,2),as.character(1:length(allsmiles)),sep = "_")
    clin <- c(clin,ident)
  }
}

smnet <- cbind(clin,sm)

write.table(smnet,"ClassyFire_InputSMILES.tsv",sep="\t",row.names = F,quote = F,col.names = F)
# ClassyFire_InputSMILES.tsv was submited to ClassyFire (http://classyfire.wishartlab.com/) [Djoumbou Feunang et al., 2016] for automated chemical classification

###############################
# retrieve ClassyFire results #
###############################

library(rjson)
library(jsonlite)

sc <- read.table("c1e1c4beb4f6426b9c9eeed05f4e7cf7.tsv",sep="\t",header = T,comment.char = "",stringsAsFactors = F,quote="") # downloaded from https://gnps.ucsd.edu/ProteoSAFe/status.jsp?task=e9e02c0ba3db473a9b1ddd36da72859b (Download Cytoscape data, clusterinfosummarygroup_attributes_withIDs_withcomponentID)

# subdivide molecular family
subclusters <- read.csv("Subclusters.csv",stringsAsFactors = F,header=F)
sc$componentindex[which(sc$cluster.index %in% subclusters$V1[which(subclusters$V2=="Subcluster2")])] <- "Subcluster2"
# subdivide molecular family

part1 <- "http://classyfire.wishartlab.com/queries/3011041.json"

files_p1 <- paste(part1, "?page=", 1:769,sep="")
files <- files_p1

jsonl <- lapply(files, function(f) tryCatch({fromJSON(f)},error=identity))

if (length(files[which(vapply(jsonl, is, logical(1), "error")==T)])!=0){
  jsonl <- jsonl[-which(vapply(jsonl, is, logical(1), "error")==T)]
}

kingdom <- list()
superclass <- list()
class <- list()
subclass <- list()
direct_parent <- list()
substituents <- list()
names_sum <- list()

for (i in 1:length(jsonl)){
  kingdom[[i]] <- jsonl[[i]]$entities$kingdom$name
  superclass[[i]] <- jsonl[[i]]$entities$superclass$name
  class[[i]] <- jsonl[[i]]$entities$class$name
  if (length(is.na(jsonl[[i]]$entities$subclass))==10){
    subclass[[i]] <- rep("none",10)
  }
  else {
    subclass[[i]] <- jsonl[[i]]$entities$subclass$name
  }
  direct_parent[[i]] <- jsonl[[i]]$entities$direct_parent$name
  substituents[[i]] <- jsonl[[i]]$entities$substituents
  names_sum[[i]] <- jsonl[[i]]$entities$identifier
}

if(length(which(is.na(unlist(names_sum))))!=0){
  kingdom <- unlist(kingdom)[-which(is.na(unlist(names_sum)))]
  kingdom[which(is.na(kingdom))] <- "none"
  superclass <- unlist(superclass)[-which(is.na(unlist(names_sum)))]
  superclass[which(is.na(superclass))] <- "none"
  class <- unlist(class)[-which(is.na(unlist(names_sum)))]
  class[which(is.na(class))] <- "none"
  subclass <- unlist(subclass)[-which(is.na(unlist(names_sum)))]
  subclass[which(is.na(subclass))] <- "none"
  direct_parent <- unlist(direct_parent)[-which(is.na(unlist(names_sum)))]
  direct_parent[which(is.na(direct_parent))] <- "none"
  substituents <- unlist(lapply(substituents, function(s) lapply(s,paste,collapse="#")))
  substituents <- substituents[-which(is.na(unlist(names_sum)))]
  names_sum <- unlist(names_sum)[-which(is.na(unlist(names_sum)))]
} else {
  kingdom <- unlist(kingdom)
  kingdom[which(is.na(kingdom))] <- "none"
  superclass <- unlist(superclass)
  superclass[which(is.na(superclass))] <- "none"
  class <- unlist(class)
  class[which(is.na(class))] <- "none"
  subclass <- unlist(subclass)
  subclass[which(is.na(subclass))] <- "none"
  direct_parent <- unlist(direct_parent)
  direct_parent[which(is.na(direct_parent))] <- "none"
  substituents <- unlist(lapply(substituents, function(s) lapply(s,paste,collapse="#")))
  substituents <- substituents
  names_sum <- unlist(names_sum)
}

networkgr <- sub("_[^_]+$", "", names_sum)
networks <- unique(networkgr)

clustersummary <- list()
clustersummary_scores <- list()

for (i in 1:length(networks)){
  
  w <- which(networkgr == networks[i])
  
  kingdom_num <- table(kingdom[w])/length(w)
  superclass_num <-table(superclass[w])/length(w)
  class_num <-table(class[w])/length(w)
  subclass_num <-table(subclass[w])/length(w)
  direct_parent_num <-table(direct_parent[w])/length(w)
  substituents_num <- table(unlist(strsplit(substituents[w],split="#")))/length(unlist(strsplit(substituents[w],split="#")))
  
  l <- list()
  l_scores <- list()
  
  l[[1]] <- names(which(kingdom_num==max(kingdom_num)))
  l[[2]] <- names(which(superclass_num==max(superclass_num)))
  l[[3]] <- names(which(class_num==max(class_num)))
  l[[4]] <- names(which(subclass_num==max(subclass_num)))
  l[[5]] <- names(which(direct_parent_num==max(direct_parent_num)))
  
  for (j in 1:length(l)){
    if(length(l[[j]]) !=1){
      l[[j]] <- paste(sort(l[[j]]),collapse = ";")
    }
    else{
      l[[j]] <- l[[j]]
    }
  }
  
  l[[6]] <- paste(sort(names(tail(sort(substituents_num),10))),collapse = ";")
  
  clustersummary[[i]] <- unlist(l)
  clustersummary_scores[[i]] <- c(round(max(kingdom_num),3),
                                  round(max(superclass_num),3),
                                  round(max(class_num),3),
                                  round(max(subclass_num),3),
                                  round(max(direct_parent_num),3),
                                  paste(round(substituents_num[which(names(substituents_num) %in% sort(names(tail(sort(substituents_num),10))))],3),collapse = ";"))
  
}

df <- data.frame(matrix(unlist(clustersummary), nrow=length(clustersummary), byrow=T))
df$CF_score <- gsub(".*_","",networks)
df$componentindex <- sub("_[^_]+$", "",networks)
colnames(df) <- c("CF_kingdom","CF_superclass","CF_class","CF_subclass","CF_directparent","CF_substituents","CF_score","componentindex")
df_scores <- data.frame(matrix(unlist(clustersummary_scores), nrow=length(clustersummary_scores), byrow=T))
df_scores$componentindex <- sub("_[^_]+$", "",networks)
colnames(df_scores) <- c("CF_kingdom_scores","CF_superclass_scores","CF_class_scores","CF_subclass_scores","CF_directparent_scores","CF_substituents_scores","componentindex")

finalout <- merge(df,df_scores,sort=FALSE)
finalout <- merge(sc[,which(colnames(sc) %in% c("cluster.index","componentindex"))],finalout,sort=FALSE,all = T)
finalout <- finalout[,c(which(colnames(finalout)=="cluster.index"),which(colnames(finalout)!="cluster.index"))]
colnames(finalout)[which(colnames(finalout)=="componentindex")] <- "CF_componentindex"

write.table(finalout,"ClassyFire_Output_forCytoscape.tsv",sep="\t",row.names = F,quote = F)

# import ClassyFire_Output_forCytoscape.tsv as table into Cytoscape version 3.4.0 [Shannon et al., 2003]. Import all data columns 
# as list of strings (data type), and semicolon as list delimiter.