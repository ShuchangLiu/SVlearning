## R source functions

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
    if(file.exists(paste(trueBEDallPre,"_",para$Type[1],".bed",sep=""))){
      print(paste("Sample has been processed, skip",sep=""))
      next
    }
    
    # no true file available, skip
    if(!file.exists(trueVCF)){
      print(paste("Sample doesn't have true file available, skip"))
      next
    }
    
    ### vcf to bed, with filtering
    
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
    
  }  # end for(sample in Sample)
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
  
  ##### for training
  if(para$knownModel==TRUE){
    print(paste("[",Sys.time(),"] Known model, skip training data preparation",sep=""))
  }else{
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
  }
  
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
  
  ## prepare the python input format
  
}


####### function for SV -> SVpiece
SV2piece <- function(SVpos, SVscore){
  # SVpos -- a SVnum*2 matrix for left and right position
  # SVscore -- a SVnum*callerNum matrix for each caller score 
  
  if(nrow(SVpos)!=nrow(SVscore)){
    stop("Error: The numbers of input SV are not equal for SVpos and SVscore")
  }
  
  if(length(SVpos)==0){
    warning("Zero SVs to be cut into pieces")
    piecePos=matrix(0,nrow=0,ncol=2)
    pieceScore=matrix(0,nrow=0,ncol=ncol(SVscore))
  }else if(nrow(SVpos)==1){  # only one line
    piecePos=SVpos
    pieceScore=SVscore
  }else{
    
    ##### sort the data by left position first, then right position
    # sort by right position
    sortInd1=order(SVpos[,2])
    SVpos=SVpos[sortInd1,]
    SVscore=SVscore[sortInd1,]
    
    # sort by left position
    sortInd2=order(SVpos[,1])
    SVpos=SVpos[sortInd2,]
    SVscore=SVscore[sortInd2,]
    
    ## skip remove duplication steps, 
    #since duplciated SVs will have the same mean value for pieces
    
    ##### SV to piece
    # initialize
    scoreNum=ncol(SVscore)
    pieceScore=c()
    piecePos=c()
    currentPiecePos=c(SVpos[1,1],0)   # current piece position, c(left, right)
    posHP=matrix(c(SVpos[1,2],1),nrow=1,ncol=2) # col: right position, index
    
    # iteration 
    for(i in 2:nrow(SVpos)){
      if(SVpos[i,1]!=SVpos[i-1,1]){
        while(SVpos[i,1]>min(posHP[,1])){ # is posHP is empty, then min(posHP) returns Inf, then the expression is false
          minInd=which.min(posHP[,1])
          if(length(piecePos)==0){
            preRight=-1  # a minus value
          }else{
            preRight=piecePos[nrow(piecePos),2]
          }
          if(posHP[minInd,1]!=preRight){   # check here later
            currentPieceScore=apply(matrix(SVscore[posHP[,2],],ncol=scoreNum),2,mean)  # use the mean score value, change later for other method
            pieceScore=rbind(pieceScore,currentPieceScore)
            currentPiecePos[2]=posHP[minInd,1]
            piecePos=rbind(piecePos,currentPiecePos)
            currentPiecePos=c(posHP[minInd,1]+1,0)
          }
          posHP=matrix(posHP[-minInd,],ncol=2)
        }
        currentPiecePos[2]=SVpos[i,1]-1
        piecePos=rbind(piecePos,currentPiecePos)
        currentPieceScore=apply(matrix(SVscore[posHP[,2],],ncol=scoreNum),2,mean) 
        pieceScore=rbind(pieceScore,currentPieceScore)
        currentPiecePos=c(SVpos[i,1],0)
      }
      posHP=rbind(posHP,matrix(c(SVpos[i,2],i),nrow=1,ncol=2))
    }
    
    # end: pop up all the heap items
    while(length(posHP)!=0){
      minInd=which.min(posHP[,1])
      if(length(piecePos)==0){
        preRight=-1  # a minus value
      }else{
        preRight=piecePos[nrow(piecePos),2]
      }
      if(posHP[minInd,1]!=preRight){   # check here later
        currentPieceScore=apply(matrix(SVscore[posHP[,2],],ncol=scoreNum),2,mean)  # use the mean score value, change later for other method
        pieceScore=rbind(pieceScore,currentPieceScore)
        currentPiecePos[2]=posHP[minInd,1]
        piecePos=rbind(piecePos,currentPiecePos)
        currentPiecePos=c(posHP[minInd,1]+1,0)
      }
      posHP=matrix(posHP[-minInd,],ncol=2)
    }
    
    # for those gap region
    naInd=which(is.na(pieceScore))
    if(length(naInd)>0){
      piecePos=matrix(piecePos[-naInd,],ncol=2)
      pieceScore=matrix(pieceScore[-naInd,],ncol=ncol(SVscore))
    }
    
    rownames(piecePos)=NULL
    rownames(pieceScore)=NULL
  }
  return(list(piecePos=piecePos,pieceScore=pieceScore))
}


#### formatPieceSVsample
# a function to be called by formatPieceSV to prepare the data
formatPieceSVsample <- function(para,sampleDir,Sample){
  
  pieceDir=paste(para$outDir,"/",para$pieceDir,sep="")
  
  #[[type]][[sample]][[as.integer(binlen)]][[as.character(chr)]]
  for(type in para$Type){
    
    if(type=="DEL"){
      binLen=para$DELbin
    }else if(type=="DUP"){
      binLen=para$DUPbin
    }else{
      stop(paste("Type ",type," not supporting yet!"))
    }
    binLenAll=c(binLen,Inf)
    
    for(sample in Sample){
      
      trueBEDtype=paste(para$outDir,"/",para$trueDir,"/",sample,"_",type,".bed",sep="")
      
      tmp=list.files(path=pieceDir,pattern=paste(type,"_",sample,sep=""))
      if(length(tmp)!=0){
        print(paste("type=",type," sample=",sample," has been pieced, skip",sep=""))
        next
      }else{
        print(paste("type=",type," sample=",sample,sep=""))
        ## read in the data
        # breakdancer
        toolBEDtype=paste(para$outDir,"/",para$breakdancerDir,"/",sample,"_",type,".bed",sep="")
        res1=read.table(toolBEDtype,sep="\t",header=F,as.is=T)
        
        # CNVnator
        toolBEDtype=paste(para$outDir,"/",para$CNVnatorDir,"/",sample,"_",type,".bed",sep="")
        res2=read.table(toolBEDtype,sep="\t",header=F,as.is=T)
        
        # delly
        toolBEDtype=paste(para$outDir,"/",para$dellyDir,"/",sample,"_",type,".bed",sep="")
        res3=read.table(toolBEDtype,sep="\t",header=F,as.is=T)
        
        # combine
        SVposAll=rbind(res1[,1:3],res2[,1:3],res3[,1:3])
        ## make this more flexiable to different number of tools in the future
        SVscoreAll=matrix(0,nrow=nrow(res1)+nrow(res2)+nrow(res3),ncol=3)
        SVscoreAll[1:nrow(res1),1]=as.numeric(res1[,5])
        SVscoreAll[((nrow(res1)+1):(nrow(res1)+nrow(res2))),2]=as.numeric(res2[,5])
        SVscoreAll[((nrow(res1)+nrow(res2)+1):(nrow(res1)+nrow(res2)+nrow(res3))),3]=as.numeric(res3[,5])
        
        SVlenAll=SVposAll[,3]-SVposAll[,2]
        
        for(binlenInd in 1:length(binLen)){
          binlen=binLen[binlenInd]
          
          for(chr in para$Chr){
            
            pieceFileXY=paste(pieceDir,"/",type,"_",sample,"_",as.integer(binlen),"_",chr,"_xy.txt",sep="")
            selectInd=which(SVlenAll>=binLenAll[binlenInd] & SVlenAll<binLenAll[binlenInd+1] & SVposAll==chr)
            if(length(selectInd)>0){
              pieceOut=SV2piece(SVpos=as.matrix(SVposAll[selectInd,2:3]),SVscore=matrix(SVscoreAll[selectInd,],ncol=ncol(SVscoreAll)))
              out=data.frame(chr,pieceOut$piecePos,type,pieceOut$pieceScore)
              
              # intersect with truth
              if(file.exists(trueBEDtype)){
                pieceFile=paste(pieceDir,"/",type,"_",sample,"_",as.integer(binlen),"_",chr,".bed",sep="")
                pieceFileMark=paste(pieceDir,"/",type,"_",sample,"_",as.integer(binlen),"_",chr,"_mark.txt",sep="")
                
                write.table(out,file=pieceFile,sep="\t",col.names=F,row.names=F,quote=F)
                s=paste(para$bedtools," intersect -a ", pieceFile, " -b ",trueBEDtype," -c > ",pieceFileMark,sep="")
                system(s)
                
                out=read.table(pieceFileMark,header=F,sep="\t",as.is=T) # bd, CNVnator, delly, true
                
                # remove tmp files
                system(paste("rm ",pieceFile,sep=""))
                system(paste("rm ",pieceFileMark,sep=""))
              }else{
                out=cbind(out,0)  # y=0 if no known truth
              }
              
              # correct the values for ML input
              # modify CNVnator
              ind1=out[,6]<0
              ind2=out[,6]==0
              out[ind1,6]=0
              out[ind2,6]=-1
              
              # modify y
              out[out[,8]>1,8]=1
              
              write.table(out,file=pieceFileXY,sep="\t",col.names=F,row.names=F,quote=F)
              
            }else{
              # write an empty file
              write.table(c(),file=pieceFileXY,sep="\t",col.names=F,row.names=F,quote=F)
            }
          } # end for(chr in para$Chr)
        } # end for(binlenInd in 1:length(binLen))
        
      } # end check file existence
    } # for(sample in Sample)
  }  # end for(type in para$Type)
}


#### formatPieceSV
# a function to format SV to pieces for training and testing data
formatPieceSV <- function(para){
  
  system(paste("mkdir -p ",para$outDir,sep=""))
  system(paste("mkdir -p ",para$outDir,"/",para$pieceDir,sep=""))
  
  ##### for training
  if(para$knownModel==TRUE){
    print(paste("[",Sys.time(),"] Known model, skip breaking training data into pieces",sep=""))
  }else{
    print(paste("[",Sys.time(),"] Breaking training data into pieces",sep=""))
    sampleDir=para$trainData
    Sample=para$trainSample
    formatPieceSVsample(para=para,sampleDir=sampleDir,Sample=Sample)
  }
  
  ##### for testing
  print(paste("[",Sys.time(),"] Breaking testing data into pieces",sep=""))
  sampleDir=para$testData
  Sample=para$testSample
  formatPieceSVsample(para=para,sampleDir=sampleDir,Sample=Sample)
}





