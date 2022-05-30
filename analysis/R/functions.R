#------------
#  Packages
#------------

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(yaGST))
suppressPackageStartupMessages(library(doParallel))

#--------------
#  Utilities
#--------------

case_barcode <- function(aliquot) {

    return( paste0(unlist(strsplit(aliquot,"-"))[1:3], collapse = "-") )
  
}

sample_barcode <- function(aliquot) {

    return( paste0(unlist(strsplit(aliquot,"-"))[1:4], collapse = "-") )
  
}

`%notlike%` <- Negate(`%like%`)

#--------------------
# mtDNA copy number
#--------------------

calc_mtDNA_copy <- function(gDNA,mtDNA,ploidy,purity) {

    copy_number <- (mtDNA / gDNA) * ((purity * ploidy) + 2*(1-purity))
    return(copy_number)
}

#--------------------------
# single sample MMW test
#--------------------------

ssMwwGST <- function(geData, geneSet){
  
  means <- rowMeans(geData)
  sds <- apply(geData, 1, sd)
  
  ans <- foreach(ss = 1:ncol(geData)) %do% {
    currentSample <- (geData[, ss] - means)/sds
    rankedList <- sort(currentSample, decreasing = T)
    
    aMwwGST <- lapply(geneSet, function(x) mwwGST(rankedList, geneSet = x, minLenGeneSet = 5, alternative = "two.sided", verbose = F))
    aMwwGST <- aMwwGST[sapply(aMwwGST, length) != 0]
    tmp_NES <- sapply(aMwwGST, function(x) x$log.pu)
    tmp_pValue <- sapply(aMwwGST, function(x) x$p.value)
    
    ans <- list(tmp_NES = tmp_NES, tmp_pValue = tmp_pValue)
    print(ss)
    return(ans)
  }
  
  if(ncol(geneSet) == 1) {
    NES <- matrix(sapply(ans, function(x) x$tmp_NES), nrow = 1, ncol = ncol(geData))
    pValue <- matrix(sapply(ans, function(x) x$tmp_pValue), nrow = 1, ncol = ncol(geData))
  }else{
    NES <- sapply(ans, function(x) x$tmp_NES)
    pValue <- sapply(ans, function(x) x$tmp_pValue)
  }

  colnames(NES) <- colnames(pValue) <- colnames(geData)
  FDR <- t(apply(pValue, 1, function(x) p.adjust(x, method = "fdr")))
  res <- list(NES = NES, pValue = pValue, FDR = FDR)
  return(res)
}
