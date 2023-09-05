#' Scar score
#'
#' Determining genomic scar score (telomeric allelic imbalance, loss-off heterozigosity, large-scle transitions), signs of homologous recombination deficiency
#' @param seg Imput file: either a binned sequenza out file (or any other segmentation file with the following columns: chromosome, position, base.ref, depth.normal, depth.tumor, depth.ratio, Af, Bf, zygosity.normal, GC.percent, good.reads, AB.normal, AB.tumor, tumor.strand) or an allele-specific segmentation file with the following columns: 1st column: sample name, 2nd column: chromosome, 3rd column: segmentation start, 4th column: segmentation end, 5th column: total copynumber, 6th column: copy number of A allele, 7th column: copy number of B allele
#' @param reference: reference genome: either grch38 or grch37. Default is grch38.
#' @param seqz Optional parameter, set to TRUE, if input file is a binned sequenza file.
#' @param ploidy Optional parameter, may be used if the ploidy of the sample is known.
#' @param chr.in.names Optional parameter, default: TRUE, set to FALSE if input file does not contain 'chr' in chromosome names.
#' @return Output is, with the following columns: HRD	Telomeric AI	Mean size	Interstitial AI	Mean Size	Whole chr AI	Telomeric LOH	Mean size	Interstitial LOH	Mean Size	Whole chr LOH	Ploidy	Aberrant cell fraction	LST	HRDscore	adjustedHRDscore
#' @export
#' @import sequenza
#' @import data.table
#' 
args <- commandArgs(T)
source("/mnt/GenePlus002/genecloud/Org_terminal/org_52/terminal/licq_18810756054/Research/HRD/07_we7h/013_HRD_CNV/code/scarHRD/R/calc.lst.bak.R")
source("/mnt/GenePlus002/genecloud/Org_terminal/org_52/terminal/licq_18810756054/Research/HRD/07_we7h/013_HRD_CNV/code/scarHRD/R/calc.ai_new.R")
source("/mnt/GenePlus002/genecloud/Org_terminal/org_52/terminal/licq_18810756054/Research/HRD/07_we7h/013_HRD_CNV/code/scarHRD/R/preprocess.hrd.R")
source("/mnt/GenePlus002/genecloud/Org_terminal/org_52/terminal/licq_18810756054/Research/HRD/07_we7h/013_HRD_CNV/code/scarHRD/R/preprocess.seqz.R")
source("/mnt/GenePlus002/genecloud/Org_terminal/org_52/terminal/licq_18810756054/Research/HRD/07_we7h/013_HRD_CNV/code/scarHRD/R/shrink.seg.ai.wrapper.R")
source("/mnt/GenePlus002/genecloud/Org_terminal/org_52/terminal/licq_18810756054/Research/HRD/07_we7h/013_HRD_CNV/code/scarHRD/R/shrink.seg.ai.R")
source("/mnt/GenePlus002/genecloud/Org_terminal/org_52/terminal/licq_18810756054/Research/HRD/07_we7h/013_HRD_CNV/code/scarHRD/R/calc.hrd.R")
chrominfo_grch37 = read.table("/mnt/GenePlus002/genecloud/Org_terminal/org_52/terminal/licq_18810756054/Research/HRD/07_we7h/013_HRD_CNV/code/scarHRD/R/chrominfo_grch37.2.list", header = T)
chrominfo_grch38 = read.table("/mnt/GenePlus002/genecloud/Org_terminal/org_52/terminal/licq_18810756054/Research/HRD/07_we7h/013_HRD_CNV/code/scarHRD/R/chrominfo_grch38.2.list", header = T)


scar_score<-function(seg,reference = "grch38", seqz=FALSE, facets=FALSE, facets.stat=NULL, ploidy=NULL, sizelimitLOH=15e6, outputdir=NULL){

  if (is.null(outputdir)){
  outputdir=getwd()
  }
  if (reference == "grch38"){
    chrominfo = chrominfo_grch38
  } else if(reference == "grch37"){
    chrominfo = chrominfo_grch37
  } else {
    stop()
  }

  purity <- NULL
  emflags <- NULL
  if (seqz==TRUE){
    cat('Preprocessing started...\n')
    seg<-preprocess.seqz(seg,ploidy0=ploidy)
    cat('Preprocessing finished \n')
  } else if(facets==TRUE){
    cat('Using Facets result, Preprocessing started...\n')
    if (! is.null(facets.stat)){
        tmp <- read.table(facets.stat,header=F,stringsAsFactors = F, sep="\t" )
        if (tmp$V1[3] == "ploidy" && tmp$V2[3] != "NA"){
            cat ("I:ploidy is ", tmp$V2[3], "\n")
            ploidy <- as.numeric(tmp$V2[3])
        }
        if (is.na(tmp$V2[5])){
            cat ("W:possible unreliable facets result, ", tmp$V2[5], "\n")
        }
        purity <- tmp$V2[2]
        emflags <- tmp$V2[5]
    }
    seg<-preprocess.facets(seg, ploidy=ploidy)
    cat('Preprocessing finished \n')
  }else {
    seg<-read.table(seg,header=T, check.names = F, stringsAsFactors = F, sep="\t")
    seg[,9]<-seg[,8]
    seg[,8]<-seg[,7]
    seg[,7]<-seg[,6]
    seg[,6]<-seg[,5]
    seg[,10]<-rep(1,dim(seg)[1])

  }
  #prep
  cat('Determining HRD-LOH, LST, TAI \n')
  seg<-preprocess.hrd(seg)
  #Calculating the hrd score:
  out_loh <- calc.hrd(seg,sizelimit1=sizelimitLOH)
  res_hrd = out_loh[[1]]
  segLOH = out_loh[[2]]
  #Calculating the telomeric allelic imbalance score:
  out_tai<- calc.ai_new(seg = seg, chrominfo = chrominfo, min.size = 11e6) #<-- the first column is what I need
  res_ai <- out_tai[[1]]
  seg_ai <- out_tai[[2]]
  #Calculating the large scale transition scores:
  out_lst <- calc.lst(seg = seg, chrominfo = chrominfo) #<-- You need to use the chrominfo.snp6 file! Nicolai sent it to you!
  res_lst <- out_lst[[1]]
  seg_lst <- out_lst[[2]]
  sum_HRD0<-res_lst+res_hrd+res_ai[1]
  if (is.null(ploidy)){
    sum_HRDc<-NA
  } else {
    sum_HRDc<-res_lst-15.5*ploidy+res_hrd+res_ai[1]
  }

  #HRDresulst<-c(res_hrd,res_ai,res_lst,sum_HRD0,sum_HRDc)
  #names(HRDresulst)<-c("HRD",colnames(res_ai),"LST", "HRD-sum","adjusted-HRDsum")
  HRDresulst<-as.character(c(res_hrd,res_ai[1],res_lst,sum_HRD0, toString(purity), toString(ploidy), toString(emflags)))
  names(HRDresulst)<-c("HRD-LOH",colnames(res_ai)[1],"LST", "HRD-sum", "cnv-purity", "cnv-ploidy", "cnv-emflags")
  run_name<-names(sum_HRD0)
  write.table(t(HRDresulst),paste0(outputdir,"/",run_name,"_HRDresults.txt"),col.names=NA,sep="\t",row.names=unique(seg[,1]),quote = FALSE)
  write.table(segLOH,paste0(outputdir,"/",run_name,"_cnv.LOH.txt"),sep="\t",quote = FALSE,row.names=FALSE)
  write.table(seg_ai,paste0(outputdir,"/",run_name,"_cnv.TAI.txt"),sep="\t",quote = FALSE,row.names = FALSE)
  write.table(seg_lst,paste0(outputdir,"/",run_name,"_cnv.LST.txt"),sep="\t",quote = FALSE,row.names = FALSE)
  return(t(HRDresulst))
}

scar_score(args[1],reference = "grch37",seqz=TRUE,outputdir=args[2])
