#' unite neighbouring segments if possible
#'
#' @param seg segmentation data
shrink.seg.ai.wrapper<-function(seg){
  new.seg <- seg[1,]
  for(j in unique(seg[,1])){
    sample.seg <- seg[seg[,1] %in% j,]
    new.sample.seg <- seg[1,]
    for(i in unique(sample.seg[,2])){
      sample.chrom.seg <- sample.seg[sample.seg[,2] %in% i,,drop=F]
      if(nrow(sample.chrom.seg) > 1){
        sample.chrom.seg <- shrink.seg.ai(sample.chrom.seg)
      }
      new.sample.seg <- rbind(new.sample.seg, sample.chrom.seg)
    }
    new.seg <- rbind(new.seg, new.sample.seg[-1,])
  }
  seg <- new.seg[-1,]
  return(seg)
}
