#' Preprocessing for further analysis
#'
#' @param seg segmentation data
#' @return preprocessed data
prerocess.hrd<-function(seg){
  #remove sex chromosomes, these are unreliable in males
  seg <- seg[!seg[,2] %in% c(paste('chr',c('X','Y','x','y',23,24),sep=''),c('X','Y','x','y',23,24)),]

  #sample names
  seg[,1] <- as.character(seg[,1])

  if(! all(seg[,8] <= seg[,7]) ){
    cat("Warning!! nB  not always <= nA!!  -- Correcting for internal use (only!)\n") # In case ASCAT people change the algorithm
    tmp <- seg
    seg[tmp[,8] > tmp[,7],7]  <- tmp[tmp[,8] > tmp[,7],8]
    seg[tmp[,8] > tmp[,7],8]  <- tmp[tmp[,8] > tmp[,7],7]
  }
  #Join segments with same copy number, after moving to ASCAT 2.3, this becomes necessary, as Raw segments gives multiple with same cn.
  seg <- shrink.seg.ai.wrapper(seg)
  return(seg)

}
