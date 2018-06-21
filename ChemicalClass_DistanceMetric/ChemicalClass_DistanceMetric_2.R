########################################################################################
#                                                                                      #
# Calculate chemically informed distance metric based on ClassyFire chemical ontology  #
#                                                                                      #
########################################################################################

# R script adapted from Junker (2018), Chemoecology 28, 29-37. (https://link.springer.com/article/10.1007/s00049-017-0250-4)

# Load source functions
source("SI1_sourceFunctions_BioSynDist.R") # retrieved from from Junker (2018), Chemoecology 28, 29-37. (https://link.springer.com/article/10.1007/s00049-017-0250-4)

#-----------------------------------------------------------------------------------------
# Load DATA - hypothetical example

BioSyn <- read.table("classlist_directparents.tsv", header=T, row.names=1,sep="\t") # or alternatively, load "classlist_subclass.tsv"
BioSyn <- BioSyn[-63,]

SampleCompoundMatrix <- read.table("featuretable_directparent_cutoff1000.tsv", header=T, row.names=1,sep="\t") # or alternatively, load "featuretable_subclasses_cutoff1000.tsv"
SampleCompoundMatrix <- SampleCompoundMatrix[which(rownames(SampleCompoundMatrix) %in% rownames(BioSyn)),]
SampleCompoundMatrix <- t(SampleCompoundMatrix)

# Calculate chemically informed distances d(A,B) between samples
r <- BioSynDist(BioSyn, SampleCompoundMatrix)
bsDist <- as.dist(r$BSDist)
bsynmat <- as.matrix(bsDist)
write.table(bsynmat,"BioSynDist_cutoff1000_directparents.tsv",sep = "\t",quote = F,row.names = T)

# Plot dendrogram of chemical classes based on the hierarchical levels of the ClassyFire chemical ontology 
pdf(file="DirectParentTree.pdf", width=8, height=8)
plot(as.phylo(r$dend_Comp),tip.color = "black", label.offset=0.02, no.margin=TRUE, cex=0.3)
dev.off()