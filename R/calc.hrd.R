#' Determnine number of loss-off heterosigosity
#'
#' @param seg segmentation data
#' @param nA column number of copy number of A allele
#' @return number of LOH
calc.hrd<-function(seg, nA=7){
  #seg: segmentation file.
  #2nd column: chromosome name
  #nA: is the column where copy number of A allele is found, default:7.
  nB <- nA+1

  output <- rep(0, length(unique(seg[,1])))
  names(output) <- unique(seg[,1])
  for(i in unique(seg[,1])){
    segSamp <- seg[seg[,1] %in% i,]
    chrDel <-vector()

    for(j in unique(segSamp[,2])){
      if(all(segSamp[segSamp[,2] == j,nB] == 0)) {
        chrDel <- c(chrDel, j)
      }
    }

    ## Added 2015-06-02, joins adjacent regions of LOH
    segSamp[segSamp[,nA] > 1,nA] <- 1
    segSamp <- shrink.seg.ai.wrapper(segSamp)
    ##
    segLOH <- segSamp[segSamp[,nB] == 0 & segSamp[,nA] != 0,,drop=F]
    segLOH <- segLOH[segLOH[,4]-segLOH[,3] > 15e6,,drop=F]
    segLOH <- segLOH[!segLOH[,2] %in% chrDel,,drop=F]
    output[i] <- nrow(segLOH)
  }
  return(output)
}
