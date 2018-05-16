#' Determine number of loss-off heterosigosity
#'
#' @param seg segmentation data
#' @param nA column number of copy number of A allele
#' @param sizelimit lower limit of the size of LOHs
#' @return number of LOH
calc.hrd<-function(seg, nA=7, return.loc=FALSE,sizelimit1){
  nB <- nA+1
  output <- rep(0, length(unique(seg[,1])))
  names(output) <- unique(seg[,1])
  if(return.loc) {
    out.seg <- matrix(0,0,9)
    colnames(out.seg) <- c(colnames(seg)[1:8],'HRD breakpoint')
  }
  #For multiple patients
  for(i in unique(seg[,1])){
    segSamp <- seg[seg[,1] %in% i,]
    chrDel <-vector()
    for(j in unique(segSamp[,2])){
      if(all(segSamp[segSamp[,2] == j,nB] == 0)) {
        chrDel <- c(chrDel, j)
      }
    }
    segSamp[segSamp[,nA] > 1,nA] <- 1
    segSamp <- shrink.seg.ai.wrapper(segSamp)
    segLOH <- segSamp[segSamp[,nB] == 0 & segSamp[,nA] != 0,,drop=F]
    segLOH <- segLOH[segLOH[,4]-segLOH[,3] > sizelimit1,,drop=F]
    segLOH <- segLOH[!segLOH[,2] %in% chrDel,,drop=F]
    output[i] <- nrow(segLOH)
    if(return.loc){
      if(nrow(segLOH) < 1){next}
      segLOH <- cbind(segLOH[,1:8], 1)
      colnames(segLOH)[9] <- 'HRD breakpoint'
      out.seg <- rbind(out.seg, segLOH)
    }
  }
  if(return.loc){
    return(out.seg)
  } else {
    return(output)
  }
}

