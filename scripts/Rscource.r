
########## prepare function ############

#### prepare_truth
# to prepare truth data: filtering, split by SV type
prepare_true <- function(para,sampleDir,Sample){
  for(sample in Sample){
    print(sample)
    
    trueVCFtmp=list.files(path=paste(sampleDir,"/",sample,sep=""),pattern="true.vcf")
    trueVCF=paste(sampleDir,"/",sample,"/",trueVCFtmp,sep="")
    trueBEDallPre=paste(para$outDir,"/",para$trueDir,"/",sample,sep="")
    
    # if the file has been processed, no need to processed again
    if(!file.exists(paste(trueBEDallPre,"_",para$Type[1],".bed",sep=""))){
      # vcf to bed, with filtering
      
      # read in vcf file
      invcf=read.table(trueVCF,sep="\t",as.is=T,comment.char="#",header=F)
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
      
      # filter the data by chr
      keepInd=which(as.character(df$chr)%in%para$Chr)
      df=df[keepInd,]
      
      # split the data by type
      for(type in para$Type){
        dfSub=df[df$name==type,]
        write.table(dfSub, file=paste(trueBEDallPre,"_",type,".bed",sep=""),sep="\t",col.names=F, row.names=F, quote=F)
      }
      
    }
  }
}


#### prepare_breakdancer
# to prepare breakdancer data: filtering, split by SV length bin and SV type
prepare_breakdancer <- function(para,sampleDir,Sample){
  
  ## for test, parameter
  #sampleDir=para$trainDataDir
  #Sample=para$trainSample
  
  tool="breakdancer"
  toolDir=paste(para$outDir,"/",para$breakdancerDir,sep="")
  trueDir=paste(para$outDir,"/",para$trueDir,sep="")
  
  for(sample in Sample){
    print(sample)
    
    toolVCFtmp=list.files(path=paste(sampleDir,"/",sample,sep=""),pattern=paste(tool,".vcf",sep=""))
    toolVCF=paste(sampleDir,"/",sample,"/",toolVCFtmp,sep="")
    toolBED1=paste(toolDir,"/",sample,".1.bed",sep="")
    toolBEDpre=paste(toolDir,"/",sample,sep="")
    
    # if the file has been processed, no need to processed again
    if(!file.exists(paste(toolBEDpre,"_",para$Type[1],".bed",sep=""))){
      
      # read in vcf file
      invcf=read.table(toolVCF,sep="\t",as.is=T,comment.char="#",header=F)
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
      
      df=data.frame(chr=invcf$CHROM,startPos=startPos,endPos=endPos,name=SVtype,score=invcf$QUAL)
      
      ## check endPos>startPos
      switchInd=which(endPos<startPos)
      tmp=df$startPos[switchInd]
      df$startPos[switchInd]=df$endPos[switchInd]
      df$endPos[switchInd]=tmp
      
      # filter the data by chr
      keepInd=which(as.character(df$chr)%in%para$Chr)
      df=df[keepInd,]
      
      # filter the data by type
      keepInd=which(df$name%in%para$Type)
      df=df[keepInd,]
      
      # filter the data with start=end
      keepInd=which(df$startPos < df$endPos)
      df=df[keepInd,]
      
      # filter the N region
      if(para$filterNregion==TRUE){
        write.table(df[,1:3],file=toolBED1,sep="\t",col.names=F, row.names=F, quote=F)
        Ncount=as.numeric(system(paste(para$bedtools," getfasta -tab -fi ",para$reference," -bed ",toolBED1," | awk -F \"\\t\" '{print $2}' | awk -F \"N|n\" '{print NF-1}'",sep=""),intern=T))
        Npercent=Ncount/(df$endPos-df$startPos+1)
        keepInd=which(Npercent<para$Nrate)
        df=df[keepInd,]
      }
      
      # split the data by type
      for(type in para$Type){
        
        toolBEDtype=paste(toolBEDpre,"_",type,".bed",sep="")
        trueBEDtype=paste(trueDir,"/",sample,"_",type,".bed",sep="")
        toolTrueMarktype=paste(toolDir,"/",sample,"_",type,"_mark.txt",sep="")
        
        dfSub=df[df$name==type,]
        write.table(dfSub, file=toolBEDtype,sep="\t",col.names=F, row.names=F, quote=F)
        
        # intersect with truth
        s=paste(para$bedtools," intersect -a ", toolBEDtype, " -b ",trueBEDtype," -f 0.5 -r -c > ",toolTrueMarktype,sep="")
        system(s)
      }
    }else{
      print(paste("Sample has been processed, skip",sep=""))
    }
  }  # end for(sample in Sample)
} # end prepare_breakdancer function



