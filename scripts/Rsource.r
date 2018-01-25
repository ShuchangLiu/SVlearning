## R source functions

#### count_N_percent
# a functiont to calculate N percentage given a string
count_N_percent <- function(S){
  Ssplit=unlist(strsplit(S,split=""))
  if(length(Ssplit)==0){
    return(1)
  }else{
    return(sum(Ssplit%in%c("N","n"))/length(Ssplit))
  }
}



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
      startPos=as.integer(invcf$POS)
      endPos=rep(0,times=nrow(invcf))
      for(i in 1:nrow(invcf)){
        if(is.na(SVlen[i])){
          endPos[i]=SVend[i]
        }else{
          endPos[i]=startPos[i]+abs(SVlen[i])-1
        }
      }
      endPos=as.integer(endPos)
      
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
      
      startPos=as.integer(invcf$POS)
      
      # # if SVlen is NA, use END as endPos; otherwise use SVlen to calculate
      # this has the potential that POS+abs(SVLEN) exceeding chromosome length
      # endPos=rep(0,times=nrow(invcf))
      # for(i in 1:nrow(invcf)){
      #   if(is.na(SVlen[i])){
      #     endPos[i]=SVend[i]
      #   }else{
      #     endPos[i]=startPos[i]+abs(SVlen[i])-1
      #   }
      # }
      # endPos=as.integer(endPos)
      
      # use SV END directly, not using SVLEN
      endPos=as.integer(SVend)
      
      df=data.frame(chr=invcf$CHROM,startPos=startPos,endPos=endPos,name=SVtype,score=invcf$QUAL,stringsAsFactors=FALSE)
      
      ## check endPos>startPos
      switchInd=which(df$endPos<df$startPos)
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
      
      # filter the data when end exceeding chromsome size
      tmpChrSize=as.numeric(para$refLen[as.character(df$chr)])
      correctInd=which(df$startPos>tmpChrSize)
      if(length(correctInd)>0){
        warning(paste(length(correctInd), " of the CNVnator SVs in sample ",sample," have start position larger than chromosome size, remove these SVs",sep=""))
        df=df[-correctInd,]
        tmpChrSize=tmpChrSize[-correctInd]
      }
      correctInd=which(df$endPos>tmpChrSize)
      if(length(correctInd)>0){
        warning(paste(length(correctInd)," of the CNVnator SVs in sample ", sample," have end position larger than chromosome size, modify them",sep=""))
        df[correctInd,"endPos"]=tmpChrSize[correctInd]
      }
      
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



#### prepare_CNVnator
# to prepare CNVnator data: filtering, split by SV length bin and SV type
prepare_CNVnator <- function(para,sampleDir,Sample){
  
  ## for test, parameter
  #sampleDir=para$trainDataDir
  #Sample=para$trainSample
  
  tool="CNVnator"
  toolDir=paste(para$outDir,"/",para$CNVnatorDir,sep="")
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
      
      # score: natorP2
      scoreTmp=sapply(tmp,function(x) return(x[startsWith(x,"natorP2=")]))
      score=as.numeric(gsub("natorP2=","",as.character(scoreTmp)))
      
      startPos=as.integer(invcf$POS)
      
      # # if SVlen is NA, use END as endPos; otherwise use SVlen to calculate
      # this has the potential that POS+abs(SVLEN) exceeding chromosome length
      # endPos=rep(0,times=nrow(invcf))
      # for(i in 1:nrow(invcf)){
      #   if(is.na(SVlen[i])){
      #     endPos[i]=SVend[i]
      #   }else{
      #     endPos[i]=startPos[i]+abs(SVlen[i])-1
      #   }
      # }
      # endPos=as.integer(endPos)
      
      # use SV END directly, not using SVLEN
      endPos=as.integer(SVend)
      
      df=data.frame(chr=invcf$CHROM,startPos=startPos,endPos=endPos,name=SVtype,score=score,stringsAsFactors=FALSE)
      
      ## check endPos>startPos
      switchInd=which(df$endPos<df$startPos)
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
      
      # filter the data when end exceeding chromsome size
      tmpChrSize=as.numeric(para$refLen[as.character(df$chr)])
      correctInd=which(df$startPos>tmpChrSize)
      if(length(correctInd)>0){
        warning(paste(length(correctInd), " of the CNVnator SVs in sample ",sample," have start position larger than chromosome size, remove these SVs",sep=""))
        df=df[-correctInd,]
        tmpChrSize=tmpChrSize[-correctInd]
      }
      correctInd=which(df$endPos>tmpChrSize)
      if(length(correctInd)>0){
        warning(paste(length(correctInd)," of the CNVnator SVs in sample ",sample," have end position larger than chromosome size, modify them",sep=""))
        df[correctInd,"endPos"]=tmpChrSize[correctInd]
      }
      
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
      print(paste("Sample has been prepared, skip",sep=""))
    }
  }  # end for(sample in Sample)
} # end prepare_CNVnator function



#### prepare_delly
# to prepare delly data: filtering, split by SV length bin and SV type
prepare_delly <- function(para,sampleDir,Sample){
  
  ## for test, parameter
  #sampleDir=para$trainDataDir
  #Sample=para$trainSample
  
  tool="delly"
  toolDir=paste(para$outDir,"/",para$dellyDir,sep="")
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
      #SVLENtmp=sapply(tmp,function(x) return(x[startsWith(x,"SVLEN=")]))
      #SVlen=as.numeric(gsub("SVLEN=","",as.character(SVLENtmp)))
      
      # score: GQ
      tmpValue=strsplit(invcf[,10],split=":")
      GQ=as.numeric(sapply(tmpValue,function(x) return(x[3])))
     
      startPos=as.integer(invcf$POS)
      
      # # if SVlen is NA, use END as endPos; otherwise use SVlen to calculate
      # this has the potential that POS+abs(SVLEN) exceeding chromosome length
      # endPos=rep(0,times=nrow(invcf))
      # for(i in 1:nrow(invcf)){
      #   if(is.na(SVlen[i])){
      #     endPos[i]=SVend[i]
      #   }else{
      #     endPos[i]=startPos[i]+abs(SVlen[i])-1
      #   }
      # }
      # endPos=as.integer(endPos)
      
      # use SV END directly, not using SVLEN
      endPos=as.integer(SVend)
      
      df=data.frame(chr=invcf$CHROM,startPos=startPos,endPos=endPos,name=SVtype,score=GQ,stringsAsFactors=FALSE)
      
      ## filter GT, only keep 0/1 or 1/1, remove 0/0 or ./.
      if(para$filterDellyGenotype==TRUE){
        GT=sapply(tmpValue,function(x) return(x[1]))
        keepInd=which(GT%in%para$keepGenotype)
        df=df[keepInd,]
      }
      
      ## check endPos>startPos
      switchInd=which(df$endPos<df$startPos)
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
      
      # filter the data when end exceeding chromsome size
      tmpChrSize=as.numeric(para$refLen[as.character(df$chr)])
      correctInd=which(df$startPos>tmpChrSize)
      if(length(correctInd)>0){
        warning(paste(length(correctInd), " of the delly SVs in sample ",sample," have start position larger than chromosome size, remove these SVs",sep=""))
        df=df[-correctInd,]
        tmpChrSize=tmpChrSize[-correctInd]
      }
      correctInd=which(df$endPos>tmpChrSize)
      if(length(correctInd)>0){
        warning(paste(length(correctInd)," of the delly SVs in sample ",sample," have end position larger than chromosome size, modify them",sep=""))
        df[correctInd,"endPos"]=tmpChrSize[correctInd]
      }
      
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
      print(paste("Sample has been prepared, skip",sep=""))
    }
  }  # end for(sample in Sample)
} # end prepare_delly function


#### prepareData
# to prepare the data for 
prepareData <- function(para){
  
  system(paste("mkdir -p ",para$outDir,sep=""))
  system(paste("mkdir -p ",para$outDir,"/",para$trueDir,sep=""))
  system(paste("mkdir -p ",para$outDir,"/",para$breakdancerDir,sep=""))
  system(paste("mkdir -p ",para$outDir,"/",para$CNVnatorDir,sep=""))
  system(paste("mkdir -p ",para$outDir,"/",para$dellyDir,sep=""))
  
  refFile=para$reference
  outDir=para$outDir
  
  ##### for training
  print(paste("[",Sys.time(),"] Prepare training data",sep=""))
  sampleDir=para$trainData
  Sample=para$trainSample
  
  # truth
  print(paste("[",Sys.time(),"] prepare true training data",sep=""))
  prepare_true(para=para,sampleDir=sampleDir,Sample=Sample)
  
  # breakdancer
  print(paste("[",Sys.time(),"] prepare breakdancer training data",sep=""))
  prepare_breakdancer(para=para,sampleDir=sampleDir,Sample=Sample)
  
  # CNVnator
  print(paste("[",Sys.time(),"] prepare CNVnator training data",sep=""))
  prepare_CNVnator(para=para,sampleDir=sampleDir,Sample=Sample)
  
  # delly
  print(paste("[",Sys.time(),"] prepare delly training data",sep=""))
  prepare_delly(para=para,sampleDir=sampleDir,Sample=Sample)
  
  
  ### for testing
  print(paste("[",Sys.time(),"] Prepare testing data",sep=""))
  sampleDir=para$testData
  Sample=para$testSample
  
  # truth
  print(paste("[",Sys.time(),"] prepare true testing data",sep=""))
  prepare_true(para=para,sampleDir=sampleDir,Sample=Sample)
  
  # breakdancer
  print(paste("[",Sys.time(),"] prepare breakdancer testing data",sep=""))
  prepare_breakdancer(para=para,sampleDir=sampleDir,Sample=Sample)
  
  # CNVnator
  print(paste("[",Sys.time(),"] prepare CNVnator training data",sep=""))
  prepare_CNVnator(para=para,sampleDir=sampleDir,Sample=Sample)
  
  # delly
  print(paste("[",Sys.time(),"] prepare delly training data",sep=""))
  prepare_delly(para=para,sampleDir=sampleDir,Sample=Sample)
  
}
