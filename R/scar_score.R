#' Scar score
#'
#' Determining genomic scar score (telomeric allelic imbalance, loss of heterozygosity, large-scale transitions), signs of homologous recombination deficiency
#' @param seg Imput file: either a binned sequenza out file (or any other segmentation file with the following columns: chromosome, position, base.ref, depth.normal, depth.tumor, depth.ratio, Af, Bf, zygosity.normal, GC.percent, good.reads, AB.normal, AB.tumor, tumor.strand) or an allele-specific segmentation file with the following columns: 1st column: sample name, 2nd column: chromosome, 3rd column: segmentation start, 4th column: segmentation end, 5th column: total copynumber, 6th column: copy number of A allele, 7th column: copy number of B allele
#' @param reference: reference genome: either grch38 or grch37. Default is grch38.
#' @param seqz Optional parameter, set to TRUE, if input file is a binned sequenza file.
#' @param ploidy Optional parameter, may be used if the ploidy of the sample is known.
#' @return Output is, with the following columns: HRD	Telomeric AI	Mean size	Interstitial AI	Mean Size	Whole chr AI	Telomeric LOH	Mean size	Interstitial LOH	Mean Size	Whole chr LOH	Ploidy	Aberrant cell fraction	LST	HRDscore	adjustedHRDscore
#' @export
#' @import sequenza
#' @import data.table
scar_score<-function(seg,reference = "grch38", seqz=FALSE, ploidy=NULL, sizelimitLOH=15e6, outputdir=NULL){

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

  if (seqz==TRUE){
    cat('Preprocessing started...\n')
    seg<-preprocess.seqz(seg,ploidy0=ploidy)
    cat('Preprocessing finished \n')
  } else {
    seg<-read.table(seg,header=T, check.names = F, stringsAsFactors = F, sep="\t")
    seg[,9]<-seg[,8]
    seg[,8]<-seg[,7]
    seg[,7]<-seg[,6]
    seg[,10]<-rep(1,dim(seg)[1])

  }
  #prep
  cat('Determining HRD-LOH, LST, TAI \n')
  seg<-preprocess.hrd(seg)
  #Calculating the hrd score:
  res_hrd <- calc.hrd(seg,sizelimit1=sizelimitLOH)
  #Calculating the telomeric allelic imbalance score:
  res_ai<- calc.ai_new(seg = seg, chrominfo = chrominfo) #<-- the first column is what I need
  #Calculating the large scale transition scores:
  res_lst <- calc.lst(seg = seg, chrominfo = chrominfo) #<-- You need to use the chrominfo.snp6 file! Nicolai sent it to you!
  sum_HRD0<-res_lst+res_hrd+res_ai[1]

  if (is.null(ploidy)){
    sum_HRDc<-NA
  } else {
    sum_HRDc<-res_lst-15.5*ploidy+res_hrd+res_ai[1]
  }

  #HRDresulst<-c(res_hrd,res_ai,res_lst,sum_HRD0,sum_HRDc)
  #names(HRDresulst)<-c("HRD",colnames(res_ai),"LST", "HRD-sum","adjusted-HRDsum")
  HRDresulst<-c(res_hrd,res_ai[1],res_lst,sum_HRD0)
  names(HRDresulst)<-c("HRD",colnames(res_ai)[1],"LST", "HRD-sum")
  run_name<-names(sum_HRD0)
  write.table(t(HRDresulst),paste0(outputdir,"/",run_name,"_HRDresults.txt"),col.names=NA,sep="\t",row.names=unique(seg[,1]))
  return(t(HRDresulst))
}
