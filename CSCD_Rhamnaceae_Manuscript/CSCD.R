
pair.cscd <- function(featureTable, dcosMatrix, idxPair) {
  ab <- matrix(featureTable[idxPair[1],], ncol=1) %*% matrix(featureTable[idxPair[2],], nrow=1)
  aa <- matrix(featureTable[idxPair[1],], ncol=1) %*% matrix(featureTable[idxPair[1],], nrow=1)
  bb <- matrix(featureTable[idxPair[2],], ncol=1) %*% matrix(featureTable[idxPair[2],], nrow=1)
  
  cosXaa = dcosMatrix * aa
  cosXbb = dcosMatrix * bb
  cosXab = dcosMatrix * ab
  
  sum(cosXab)/max(c(sum(cosXaa), sum(cosXbb)))
}

cscd <- function(featureTable, cosineEdgeList, norm=TRUE, cosineThreshold=0.6) {
	if(norm) featureTable <- t(apply(featureTable, 1, function(x) x/sum(as.numeric(x))))

	samp_pairs <- combn((1:nrow(featureTable)),2)

	cosineEdgeList[cosineEdgeList[,5] < cosineThreshold,5] <- 0
	cosineMatrix <- matrix(0, nrow=ncol(featureTable), ncol=ncol(featureTable))
	  
	cosineMatrix[cbind(match(cosineEdgeList[,1], colnames(featureTable)), match(cosineEdgeList[,2], colnames(featureTable)))] <- cosineEdgeList[,5]
	cosineMatrix[cbind(match(cosineEdgeList[,2], colnames(featureTable)), match(cosineEdgeList[,1], colnames(featureTable)))] <- cosineEdgeList[,5]
	diag(cosineMatrix) <- 1

	cscs <- apply(samp_pairs, 2, function(x) pair.cscd(featureTable, cosineMatrix, x))  
	cscd <- 1 - cscs
	cscdMatrix <- matrix(0, nrow=nrow(featureTable), ncol=nrow(featureTable))
	cscdMatrix[t(samp_pairs)] <- cscd
	cscdMatrix[t(samp_pairs[2:1,])] <- cscd
	rownames(cscdMatrix) <- colnames(cscdMatrix) <- rownames(featureTable)
	return(cscdMatrix)
}
