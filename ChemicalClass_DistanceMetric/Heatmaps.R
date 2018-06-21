############################################
#                                          #
#    Heatmap with customized dendrograms   #
#                                          #
############################################

library(dendextend)
library(vegan)
library(gplots)
library(RColorBrewer)

############################
# Subclasses bray curtis   #
############################

## import species dendrogram
dm <- read.csv("BioSynDist_cutoff1000_subclasses.tsv",check.names=F,row.names=1,sep="\t")

metadata <- read.csv("MetaData_Rhamnaceae.txt",sep="\t")

dend <-  dm %>% as.dist %>%
  hclust(method = "complete") %>% as.dendrogram

genera <- metadata$Genus[match(rownames(dm)[order.dendrogram(dend)],metadata$filename)]
  
gen_cols <- c("purple3", "orangered2","purple","slateblue1","black","steelblue",
                "tomato","red","royalblue4","rosybrown2","yellow","tan1","darkorange","green","blue")
  
true_species_cols <- gen_cols[as.numeric(genera)]
  
dend <-  dm %>% as.dist %>%
  hclust(method = "complete") %>% as.dendrogram %>% 
  set("branches_lwd", 2) %>%
  set("labels_colors","black") %>% 
  set("labels_cex", c(0.3)) 
  #color_branches(k=71,col=true_species_cols)
plot(dend)

## import chemical class dendrogram
subcl <- read.table("classlist_subclass.tsv", header=T, row.names=1,sep="\t") 

ft <- read.csv("featuretable_subclasses_cutoff1000.tsv",sep="\t",row.names = 1)
ft <- ft[which(rownames(ft) %in% rownames(subcl)),]
ft <- t(ft)
subcl <- subcl[which(rownames(subcl) %in% colnames(ft)),]
ft <- ft[,match(rownames(subcl),colnames(ft))]
ft_norm <- ft/rowSums(ft) 

bioSynDistComp <- vegdist(subcl, method="bray", binary=T,na.rm = TRUE)
clus_Comp <- hclust(bioSynDistComp, method="average") 

chemdend <-  clus_Comp %>% as.dendrogram %>% set("branches_lwd", 2) 

my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 15)

generarow <- metadata$Genus[match(gsub("[.]"," ",rownames(ft)),gsub('[^[:alnum:]]', ' ', metadata$filename))]
Rowcols <- gen_cols[as.numeric(generarow)]

colpal <- "Dark2"
chemcols <- c(rep(brewer.pal(8, colpal)[1],5),
              rep(brewer.pal(8, colpal)[2],13), 
              brewer.pal(8, colpal)[3:6], 
              rep(brewer.pal(8, colpal)[7],4),
              rep("black",1),
              rep(brewer.pal(8, colpal)[8],17))[match(colnames(ft),colnames(ft)[order.dendrogram(chemdend)])]

pdf(file="SubclassHeatMapDendrogram_ScaleRow.pdf", width=11, height=9)
heatmap.2(ft, Rowv=dend, Colv=chemdend, cexRow = 0.5,scale="row",col = my_palette,RowSideColors=Rowcols,ColSideColors=chemcols,tracecol=NA) 
dev.off()