#' Preprocessing seqz files
#'
#' Preprocesses seqz files
#' @param seg A segmentation file
#' @param ploidy0 ploidy, optional
#' @return Output is
preprocess.seqz<-function(seg, ploidy0=NULL,outputdir=NULL){
  if (is.null(ploidy0)){
    ploidy01 = seq(1, 5.5, 0.1)
  } else {
    ploidy01= seq(ploidy0-0.5,ploidy0+0.5,0.1)
  }
  
  if (is.null(outputdir)){
    outputdir = getwd()
  }
  
  run_name<-gsub(".*/","",gsub("_small.seqz","",gsub("gz","",seg)))
  extract<-sequenza.extract(seg, chromosome.list=paste('chr',c(1:24),sep=''),gamma = 60, kmin = 50)
  extract.fit<-sequenza::sequenza.fit(extract, N.ratio.filter = 10, N.BAF.filter = 1, segment.filter = 3e6, mufreq.treshold = 0.10, ratio.priority = FALSE,ploidy=ploidy01, mc.cores = 1)
  #  sequenza.results(extract, extract.fit, out.dir = getwd(),sample.id =run_name)

  seg.tab <- do.call(rbind, extract$segments[extract$chromosomes])
  seg.len <- (seg.tab$end.pos - seg.tab$start.pos)/1e+06
  cint <- get.ci(extract.fit)
  cellularity <- cint$max.cellularity
  ploidy <- cint$max.ploidy
  avg.depth.ratio <- extract$avg.depth.ratio #avg.depth.ratio <- mean(extract$gc$adj[, 2])
  info_seg<-c(cellularity,ploidy,avg.depth.ratio)
  names(info_seg)<-c("cellularity","ploidy","avg.depth.ratio")
  write.table(t(info_seg),paste0(outputdir,"/",run_name,"_info_seg.txt"),sep="\t",row.names=F)
  allele.cn <- sequenza:::baf.bayes(Bf = seg.tab$Bf, CNt.max = 20, depth.ratio = seg.tab$depth.ratio, avg.depth.ratio = 1,
                                   cellularity = cint$max.cellularity, ploidy = cint$max.ploidy,
                                   sd.ratio = seg.tab$sd.ratio, weight.ratio = seg.len, sd.Bf = seg.tab$sd.BAF,
                                   weight.Bf = 1, ratio.priority = FALSE, CNn = 2)
  seg.tab$CN <- allele.cn[,1]
  allele.cn <- as.data.table(allele.cn)
  #Making imput file
  seg <- data.frame(SampleID = as.character(run_name), Chromosome = seg.tab$chromosome, Start_position = seg.tab$start.pos,
                    End_position = seg.tab$end.pos, Nprobes = 1, total_cn = allele.cn$CNt, A_cn = allele.cn$B,
                    B_cn = allele.cn$A, ploidy = ploidy)
  seg$contamination <- 1
  seg<-seg[!is.na(seg$A_cn),]
  seg<-seg[!is.na(seg$B_cn),]
  return(seg)
}
