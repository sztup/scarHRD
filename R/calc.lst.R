#' Determine large-scale transitions
#'
#' @param seg segmentation data
#' @param chrominfo reference genome version
#' @param nA column number of copy number of A allele
#' @param chr.arm option to use chromosome arms defined during segmentation
#' @return number of LSTs
calc.lst<-function(seg, chrominfo=chrominfo,nA=7,chr.arm='no'){
  nB <- nA+1
  samples <- unique(seg[,1])
  output <- setNames(rep(0,length(samples)), samples)
  for(j in samples){
    sample.seg <- seg[seg[,1] %in% j,]
    sample.lst <- c()
    chroms <- unique(sample.seg[,2])
    chroms <- chroms[!chroms %in% c(23,24,'chr23','chr24','chrX','chrx','chrY','chry')]
    for(i in chroms){
      sample.chrom.seg <- sample.seg[sample.seg[,2] %in% i,,drop=F]
      if(chr.arm !='no'){
        p.max <- if(any(sample.chrom.seg[,chr.arm] == 'p')){max(sample.chrom.seg[sample.chrom.seg[,chr.arm] == 'p',4])}
        q.min <- min(sample.chrom.seg[sample.chrom.seg[,chr.arm] == 'q',3])
      }
      if(nrow(sample.chrom.seg) < 2) {next}
      sample.chrom.seg.new <- sample.chrom.seg
      if(chr.arm == 'no'){
        p.arm <- sample.chrom.seg.new[sample.chrom.seg.new[,3] <= chrominfo[i,2],,drop=F] # split into chromosome arms
        q.arm <- sample.chrom.seg.new[sample.chrom.seg.new[,4] >= chrominfo[i,3],,drop=F]
        q.arm<- shrink.seg.ai(q.arm)
        q.arm[1,3] <- chrominfo[i,3]
        if(nrow(p.arm) > 0){
          p.arm<- shrink.seg.ai(p.arm)
          p.arm[nrow(p.arm),4] <- chrominfo[i,2]
        }
      }
      if(chr.arm != 'no'){
        q.arm <- sample.chrom.seg.new[sample.chrom.seg.new[,chr.arm] == 'q',,drop=F]
        q.arm<- shrink.seg.ai(q.arm)
        q.arm[1,3] <- q.min
        if(any(sample.chrom.seg.new[,chr.arm] == 'p')){
          p.arm <- sample.chrom.seg.new[sample.chrom.seg.new[,chr.arm] == 'p',,drop=F] # split into chromosome arms
          p.arm<- shrink.seg.ai(p.arm)
          p.arm[nrow(p.arm),4] <- p.max
        }
      }
      n.3mb <- which((p.arm[,4] - p.arm[,3]) < 3e6)
      while(length(n.3mb) > 0){
        p.arm <- p.arm[-(n.3mb[1]),]
        p.arm <- shrink.seg.ai(p.arm)
        n.3mb <- which((p.arm[,4] - p.arm[,3]) < 3e6)
      }
      if(nrow(p.arm) >= 2){
        p.arm <- cbind(p.arm[,1:8], c(0,1)[match((p.arm[,4]-p.arm[,3]) >= 10e6, c('FALSE','TRUE'))])
        for(k in 2:nrow(p.arm)){
          if(p.arm[k,9] == 1 & p.arm[(k-1),9]==1 & (p.arm[k,3]-p.arm[(k-1),4]) < 3e6){
            sample.lst <- c(sample.lst, 1)
          }
        }
      }
      n.3mb <- which((q.arm[,4] - q.arm[,3]) < 3e6)
      while(length(n.3mb) > 0){
        q.arm <- q.arm[-(n.3mb[1]),]
        q.arm <- shrink.seg.ai(q.arm)
        n.3mb <- which((q.arm[,4] - q.arm[,3]) < 3e6)
      }
      if(nrow(q.arm) >= 2){
        q.arm <- cbind(q.arm[,1:8], c(0,1)[match((q.arm[,4]-q.arm[,3]) >= 10e6, c('FALSE','TRUE'))])
        for(k in 2:nrow(q.arm)){
          if(q.arm[k,9] == 1 & q.arm[(k-1),9]==1 & (q.arm[k,3]-q.arm[(k-1),4]) < 3e6){
            sample.lst <- c(sample.lst, 1)

          }
        }
      }
    }
    output[j] <- sum(sample.lst)
  }
  return(output)
}
