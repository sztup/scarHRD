#' Determine Allelic imbalance
#'
#' @param seg segmentation data
#' @param chrominfo reference genome version
#' @param min.size min.size
#' @param ploidyByChromosome option to determine ploidy per chromosome
#' @param shrink shrink
#' @return number of Allelic Imbalances
calc.ai_new<-function(seg, chrominfo, min.size=1e6, cont = 0,ploidyByChromosome=TRUE, shrink=TRUE){
  seg <- seg[seg[,4]- seg[,3] >= min.size,]
  seg <- seg[seg[,10] >= cont,]
  if(shrink){
    seg <- shrink.seg.ai.wrapper(seg)
  }
  AI <- rep(NA, nrow(seg))
  seg <- cbind(seg, AI)
  samples <- as.character(unique(seg[,1]))

  ascat.ploidy <- setNames(seg[!duplicated(seg[,1]),9], seg[!duplicated(seg[,1]),1])
  for(j in samples){
    sample.seg <- seg[seg[,1] %in% j,]
    if(!ploidyByChromosome){
      ploidy <- vector()
      for(k in unique(sample.seg[,6])){
        tmp <- sample.seg[sample.seg[,6] %in% k,]
        ploidy <- c(ploidy, setNames(sum(tmp[,4]-tmp[,3]), k))
      }
      ploidy <- as.numeric(names(ploidy[order(ploidy,decreasing=T)]))[1]
      sample.seg[,9] <- ploidy
      if(ploidy %in% c(1,seq(2, 200,by=2))){
        sample.seg[,'AI'] <- c(0,2)[match(sample.seg[,7] == sample.seg[,8], c('TRUE', 'FALSE'))]
      }
      if(!ploidy %in%  c(1,seq(2, 200,by=2))){
        sample.seg[,'AI'] <- c(0,2)[match(sample.seg[,7] + sample.seg[,8] == ploidy & sample.seg[,7] != ploidy, c('TRUE', 'FALSE'))]
      }
    }
    new.sample.seg<- sample.seg[1,]
    for(i in unique(sample.seg[,2])){
      sample.chrom.seg <- sample.seg[sample.seg[,2] %in% i,,drop=F]
      if(nrow(sample.chrom.seg) == 0){ next}
      if(ploidyByChromosome){
        ploidy <- vector()
        for(k in unique(sample.seg[,6])){
          tmp <- sample.chrom.seg[sample.chrom.seg[,6] %in% k,]
          ploidy <- c(ploidy, setNames(sum(tmp[,4]-tmp[,3]), k))
          ploidy <- ploidy[!names(ploidy) %in% 0]
        }
        ploidy <- as.numeric(names(ploidy[order(ploidy,decreasing=T)]))[1]
        sample.chrom.seg[,9] <- ploidy # update "ploidy" column, so the new calculated value can be returned
        if(ploidy %in% c(1,seq(2, 200,by=2))){
          sample.chrom.seg[,'AI'] <- c(0,2)[match(sample.chrom.seg[,7] == sample.chrom.seg[,8], c('TRUE', 'FALSE'))]
        }
        if(!ploidy %in%  c(1,seq(2, 200,by=2))){
          sample.chrom.seg[,'AI'] <- c(0,2)[match(sample.chrom.seg[,7] + sample.chrom.seg[,8] == ploidy & sample.chrom.seg[,8] != 0, c('TRUE', 'FALSE'))]
        }
        sample.seg[sample.seg[,2] %in% i,9] <-ploidy
        sample.seg[sample.seg[,2] %in% i,'AI'] <-sample.chrom.seg[,'AI']
      }

      if(class(chrominfo) != 'logical'){# Here we consider the centromere
        if(sample.chrom.seg[1,'AI'] == 2 & nrow(sample.chrom.seg) != 1 & sample.chrom.seg[1,4] < (chrominfo[i,2])){
          sample.seg[sample.seg[,2]==i,'AI'][1] <- 1
        }
        if(sample.chrom.seg[nrow(sample.chrom.seg),'AI'] == 2 & nrow(sample.chrom.seg) != 1 & sample.chrom.seg[nrow(sample.chrom.seg),3] > (chrominfo[i,3])){
          sample.seg[sample.seg[,2]==i,'AI'][nrow(sample.seg[sample.seg[,2]==i,])] <- 1
        }
      }
      if(nrow(sample.seg[sample.seg[,2]==i,]) == 1 & sample.seg[sample.seg[,2]==i,'AI'][1] != 0){
        sample.seg[sample.seg[,2]==i,'AI'][1] <- 3
      }
    }
    seg[seg[,1] %in% j,] <- sample.seg
  }
  samples <- as.character(unique(seg[,1]))
  #0 = no AI, 1=telomeric AI, 2=interstitial AI, 3= whole chromosome AI
  no.events <- matrix(0, nrow=length(samples), ncol=12)
  rownames(no.events) <- samples
  colnames(no.events) <- c("Telomeric AI", "Mean size", "Interstitial AI", "Mean Size", "Whole chr AI", "Telomeric LOH",  "Mean size", "Interstitial LOH", "Mean Size", "Whole chr LOH", "Ploidy", "Aberrant cell fraction")
  for(j in samples){
    sample.seg <- seg[seg[,1] %in% j,]
    no.events[j,1] <- nrow(sample.seg[sample.seg[,'AI'] == 1,])
    no.events[j,2] <- mean(sample.seg[sample.seg[,'AI'] == 1,4] - sample.seg[sample.seg[,'AI'] == 1,3])
    no.events[j,3] <- nrow(sample.seg[sample.seg[,'AI'] == 2,])
    no.events[j,4] <- mean(sample.seg[sample.seg[,'AI'] == 2,4] - sample.seg[sample.seg[,'AI'] == 2,3])
    no.events[j,5] <- nrow(sample.seg[sample.seg[,'AI'] == 3,])
    no.events[j,11] <- ascat.ploidy[j]
    no.events[j,12] <- unique(sample.seg[,10]) # aberrant cell fraction
    #Here we restrict ourselves to real LOH
    sample.seg <- sample.seg[sample.seg[,8] == 0,]
    no.events[j,6] <- nrow(sample.seg[sample.seg[,'AI'] == 1,])
    no.events[j,7] <- mean(sample.seg[sample.seg[,'AI'] == 1,4] - sample.seg[sample.seg[,'AI'] == 1,3])
    no.events[j,8] <- nrow(sample.seg[sample.seg[,'AI'] == 2,])
    no.events[j,9] <- mean(sample.seg[sample.seg[,'AI'] == 2,4] - sample.seg[sample.seg[,'AI'] == 2,3])
    no.events[j,10] <- nrow(sample.seg[sample.seg[,'AI'] == 3,])
  }
  return(no.events)
}
