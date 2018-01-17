## source function 

### vcf2bed
# a function to transfer input vcf to output bed format
## add CNVnator-DUP-overlapped INV into DUP for breakdancer and delly

vcf2bed <-function(inputVCF,outputBEDpre,Type=c("DEL","DUP"),filterChr=TRUE,keepChr=1:22,filterLong=FALSE,longCutoff=500000){
  
  invcf=read.table(inputVCF,sep="\t",as.is=T,comment.char="#",header=F)
  colnames(invcf)[1:8]=c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO")
  
  tmp=strsplit(invcf$INFO,split=";")
  
  # SVtype
  SVTYPEtmp=sapply(tmp,function(x) return(x[startsWith(x,"SVTYPE=")]))
  SVtype=gsub("SVTYPE=","",as.character(SVTYPEtmp))
  
  # END
  SVENDtmp=sapply(tmp,function(x) return(x[startsWith(x,"END=")]))
  SVend=as.numeric(gsub("END=","",as.character(SVENDtmp)))
  
  # SVlen
  SVLENtmp=sapply(tmp,function(x) return(x[startsWith(x,"SVLEN=")]))
  SVlen=as.numeric(gsub("SVLEN=","",as.character(SVLENtmp)))
  
  # if SVlen is NA, use END as endPos; otherwise use SVlen to calculate
  startPos=invcf$POS
  endPos=rep(0,times=nrow(invcf))
  for(i in 1:nrow(invcf)){
    if(is.na(SVlen[i])){
      endPos[i]=SVend[i]
    }else{
      endPos[i]=startPos[i]+abs(SVlen[i])-1
    }
  }

  df=data.frame(chr=invcf$CHROM,startPos=startPos,endPos=endPos,name=SVtype)
  
  ## check endPos>startPos
  switchInd=which(endPos<startPos)
  tmp=df$startPos[switchInd]
  df$startPos[switchInd]=df$endPos[switchInd]
  df$endPos[switchInd]=tmp
  
  # filter out END-POS > cutoff
  if(filterLong==TRUE){
    longInd=(df$endPos-df$startPos)>longCutoff
    df=df[which(longInd==FALSE),]
  }

  # split the data by chr
  if(filterChr==TRUE){
    keepInd=which(df$chr%in%keepChr)
    df=df[keepInd,]
  }
  
  # split the data by type
  for(type in Type){
    dfSub=df[df$name==type,]
    write.table(dfSub, file=paste(outputBEDpre,"_",type,".bed",sep=""),sep="\t",col.names=F, row.names=F, quote=F)
  }
  
}


