## R source functions

### checkFormatPara
# to check and format the paramters
checkFormatPara <- function(configFileName){
  
  configFile=read.table(configFileName,comment.char="#",sep="\t",as.is=T,header=F)
  
  para=list()
  for(i in 1:nrow(configFile)){
    para[[configFile[i,1]]]=unlist(strsplit(gsub(" ","",configFile[i,2]),split=","))
  }
  
  # change format
  para$Nrate=as.numeric(para$Nrate)
  para$vote=as.numeric(para$vote)
  para$DELbin=as.integer(para$DELbin)
  para$DUPbin=as.integer(para$DUPbin)
  para$Chr=as.character(para$Chr)
  para$pieceBalanceRate=as.numeric(para$pieceBalanceRate)
  para$PgapDis=as.numeric(para$PgapDis)
  para$SVgapDis=as.numeric(para$SVgapDis)
  para$markPieceRate=as.numeric(para$markPieceRate)
  
  
  # ref len file
  para$refLen=rep(0,times=length(para$Chr))
  names(para$refLen)=para$Chr
  tmp=read.table(para$refLenFile,sep="\t",header=F,as.is=T)
  for(chr in para$Chr){
    para$refLen[chr]=as.integer(tmp[which(tmp[,1]==chr),2])
  }
  
  ## training sample, change here in case don't want to use the full training samples
  para$trainSample=list.files(para$trainDataDir)
  ## testing sample, change here in case don't want to use the full testing samples
  para$testSample=list.files(para$testDataDir)
  
  ##### TBD: add code to check the parameter inputs, whether they satify the criteria
  
  # check para$MLmethod
  if(!(all(para$MLmethod%in%c("SVMradial","SVMpolynomial","RF","NN","LDA","adaboost")))){
    stop("Error: In config file, MLmethod can only choose from SVMpolynomial, SVMrdial, RF, NN, LDA, adaboost")
  }
  if(length(para$MLmethod)<para$vote){
    stop("Error: The number of MLmethod has to be equal to or greater than vote value.")
  }
  
  # check para$Type
  if(!all(para$Typ %in% c("DEL","DUP"))){
    stop("Error: In config file, Type can only choose from DEL and DUP")
  }
  
  ### to test all the MLmethod
  #para$MLmethod=c("NN","SVMpolynomial","SVMradial","LDA","RF","adaboost")
  
  # check scoreID  ## @@@@@@ TBD: change here: combine same toolName
  scoreIDtmp=list()
  allToolName=c()
  for(i in 1:length(para$scoreID)){
    tmp=unlist(strsplit(para$scoreID[i],split=":"))
    if(length(tmp)!=3){
      stop("Invalide paramter setting for scoreID")
    }
    if(!tmp[2]%in%c("6","8","10")){
      stop("scoreID only accept column 6, 8 and 10")
    }
    
    currentName=paste(tmp[1],tmp[3],sep="_")
    if(currentName%in%allToolName){
      # don't allow duplicate names
      warning(paste("Duplicate toolName and score for ",currentName," in para$scoreID"))  
      next
    }
    allToolName=c(allToolName,currentName)
    
    if(tmp[1]%in%names(scoreIDtmp)){
      # tool name already exit
      scoreIDtmp[[tmp[1]]]=rbind(scoreIDtmp[[tmp[1]]],data.frame(col=tmp[2],ID=tmp[3],stringsAsFactors=FALSE))
    }else{
      # new tool name
      scoreIDtmp[[tmp[1]]]=data.frame(col=tmp[2],ID=tmp[3],stringsAsFactors=FALSE)
    }
  }
  para$scoreID=scoreIDtmp
  allToolName=c()
  for(i in 1:length(para$scoreID)){
    allToolName=c(allToolName,paste(names(para$scoreID)[i],para$scoreID[[i]]$ID,sep="_"))
  }
  para$allToolName=allToolName
  
  
  print(para)
  
  return(para)
  
}


### installRpackage
# code to install R packages
installRpackage <- function(para){
  # library("nnet")  ## nnet
  # library("e1071")  ## SVM
  # library("MASS")  ## LDA
  # library("randomForest")  ## random forest
  # library("fastAdaboost")  ## adaboost
  
  print(paste("[",Sys.time(),"] Load R packge",sep=""))
  
  modelPackageMap=c("nnet","e1071","e1071","MASS","randomForest","fastAdaboost")
  names(modelPackageMap)=c("NN","SVMpolynomial","SVMradial","LDA","RF","adaboost")
  
  packageList=unique(modelPackageMap[para$MLmethod])
  
  for(packageName in packageList){
    print(packageName)
    if(!packageName%in%rownames(installed.packages())){
      install.packages(as.character(packageName))
    }
    library(as.character(packageName),character.only=T)
  }
  
}


########## prepare function ############

#### prepare_truth
# to prepare truth data: filtering, split by SV type
prepare_true <- function(para,sampleDir,Sample){
  for(sample in Sample){
    print(sample)
    
    trueVCFtmp=list.files(path=paste(sampleDir,"/",sample,sep=""),pattern="true.vcf")
    
    if(length(trueVCFtmp)==0){
      print(paste("Sample ",sample," doesn't have true file available, skip",sep=""))
      next
    }else if(length(trueVCFtmp)>1){
      stop(paste("Multiple true files for sample ",sample,sep=""))
    }
    
    trueVCF=paste(sampleDir,"/",sample,"/",trueVCFtmp,sep="")
    trueBEDallPre=paste(para$tmpDir,"/",para$trueDir,"/",sample,sep="")
    
    # if the file has been processed, no need to processed again
    if(file.exists(paste(trueBEDallPre,"_",para$Type[1],".bed",sep=""))){
      print(paste("Sample has been processed, skip",sep=""))
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
    
    df=data.frame(chr=invcf$CHROM,startPos=startPos,endPos=endPos,name=SVtype,stringsAsFactors=FALSE)
    df$startPos=as.integer(df$startPos)
    df$endPos=as.integer(df$endPos)
    
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
  toolDir=paste(para$tmpDir,"/",para$breakdancerDir,sep="")
  #trueDir=paste(para$tmpDir,"/",para$trueDir,sep="")
  
  for(sample in Sample){
    print(sample)
    
    toolVCFtmp=list.files(path=paste(sampleDir,"/",sample,sep=""),pattern=paste(tool,".vcf",sep=""))
    if(length(toolVCFtmp)==0){
      stop(paste("No VCF file available for tool=",tool,", sample=",sample,sep=""))
    }else if(length(toolVCFtmp)>1){
      stop(paste("Multiple VCF files available for tool=",tool,", sample=",sample,sep=""))
    }
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
      df$startPos=as.integer(df$startPos)
      df$endPos=as.integer(df$endPos)
      
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
      
      # filter the data when end exceeding chromsome size or start smaller than 0
      tmpChrSize=as.numeric(para$refLen[as.character(df$chr)])
      correctInd=which(df$startPos>tmpChrSize | df$startPos<=0)
      if(length(correctInd)>0){
        warning(paste(length(correctInd), " of the ",tool," SVs in sample ",sample," have start position larger than chromosome size OR smaller than 0, remove these SVs",sep=""))
        df=df[-correctInd,]
        tmpChrSize=tmpChrSize[-correctInd]
      }
      correctInd=which(df$endPos>tmpChrSize)
      if(length(correctInd)>0){
        warning(paste(length(correctInd)," of the ",tool," SVs in sample ", sample," have end position larger than chromosome size, modify them",sep=""))
        df[correctInd,"endPos"]=tmpChrSize[correctInd]
      }
      
      # write the original bed 
      for(type in para$Type){
        dfTmp=df[df$name==type,]
        fileTmp=paste(toolBEDpre,"_",type,"_ori.bed",sep="")
        write.table(dfTmp,file=fileTmp,sep="\t",col.names=F,row.names=F,quote=F)
      }
      
      # filter the N region
      if(para$filterNregion==TRUE){
        write.table(data.frame(df[,1],as.integer(df[,2]-1),as.integer(df[,3])),file=toolBED1,sep="\t",col.names=F, row.names=F, quote=F)
        Ncount=as.numeric(system(paste(para$bedtools," getfasta -tab -fi ",para$reference," -bed ",toolBED1," | awk -F \"\\t\" '{print $2}' | awk -F \"N|n\" '{print NF-1}'",sep=""),intern=T))
        if(length(Ncount)!=nrow(df)){
          stop(paste("File ",toolBED1, " is wrong input for getfasta tool for sample=",sample," tool=breakdancer",sep=""))
        }
        Npercent=Ncount/(df$endPos-df$startPos+1)
        keepInd=which(Npercent<para$Nrate)
        df=df[keepInd,]
      }
      
      # split the data by type
      for(type in para$Type){
        
        toolBEDtype=paste(toolBEDpre,"_",type,".bed",sep="")
        #trueBEDtype=paste(trueDir,"/",sample,"_",type,".bed",sep="")
        #toolTrueMarktype=paste(toolDir,"/",sample,"_",type,"_mark.txt",sep="")
        
        dfSub=df[df$name==type,]
        write.table(dfSub, file=toolBEDtype,sep="\t",col.names=F, row.names=F, quote=F)
        
        # intersect with truth
        #s=paste(para$bedtools," intersect -a ", toolBEDtype, " -b ",trueBEDtype," -f 0.5 -r -c > ",toolTrueMarktype,sep="")
        #system(s)
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
  toolDir=paste(para$tmpDir,"/",para$CNVnatorDir,sep="")
  #trueDir=paste(para$tmpDir,"/",para$trueDir,sep="")
  
  for(sample in Sample){
    print(sample)
    
    toolVCFtmp=list.files(path=paste(sampleDir,"/",sample,sep=""),pattern=paste(tool,".vcf",sep=""))
    if(length(toolVCFtmp)==0){
      stop(paste("No VCF file available for tool=",tool,", sample=",sample,sep=""))
    }else if(length(toolVCFtmp)>1){
      stop(paste("Multiple VCF files available for tool=",tool,", sample=",sample,sep=""))
    }
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
      df$startPos=as.integer(df$startPos)
      df$endPos=as.integer(df$endPos)
      
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
      
      # filter the data when end exceeding chromsome size or start smaller than 0
      tmpChrSize=as.numeric(para$refLen[as.character(df$chr)])
      correctInd=which(df$startPos>tmpChrSize | df$startPos<=0)
      if(length(correctInd)>0){
        warning(paste(length(correctInd), " of the ",tool," SVs in sample ",sample," have start position larger than chromosome size OR smaller than 0, remove these SVs",sep=""))
        df=df[-correctInd,]
        tmpChrSize=tmpChrSize[-correctInd]
      }
      correctInd=which(df$endPos>tmpChrSize)
      if(length(correctInd)>0){
        warning(paste(length(correctInd)," of the ",tool," SVs in sample ",sample," have end position larger than chromosome size, modify them",sep=""))
        df[correctInd,"endPos"]=tmpChrSize[correctInd]
      }
      
      # write the original bed 
      for(type in para$Type){
        dfTmp=df[df$name==type,]
        fileTmp=paste(toolBEDpre,"_",type,"_ori.bed",sep="")
        write.table(dfTmp,file=fileTmp,sep="\t",col.names=F,row.names=F,quote=F)
      }
      
      # filter the N region
      if(para$filterNregion==TRUE){
        write.table(data.frame(df[,1],as.integer(df[,2]-1),as.integer(df[,3])),file=toolBED1,sep="\t",col.names=F, row.names=F, quote=F)
        Ncount=as.numeric(system(paste(para$bedtools," getfasta -tab -fi ",para$reference," -bed ",toolBED1," | awk -F \"\\t\" '{print $2}' | awk -F \"N|n\" '{print NF-1}'",sep=""),intern=T))
        if(length(Ncount)!=nrow(df)){
          stop(paste("File ",toolBED1, " is wrong input for getfasta tool for sample=",sample," tool=CNVnator",sep=""))
        }
        Npercent=Ncount/(df$endPos-df$startPos+1)
        keepInd=which(Npercent<para$Nrate)
        df=df[keepInd,]
      }
      
      # split the data by type
      for(type in para$Type){
        
        toolBEDtype=paste(toolBEDpre,"_",type,".bed",sep="")
        #trueBEDtype=paste(trueDir,"/",sample,"_",type,".bed",sep="")
        #toolTrueMarktype=paste(toolDir,"/",sample,"_",type,"_mark.txt",sep="")
        
        dfSub=df[df$name==type,] # still working when it's empty
        write.table(dfSub, file=toolBEDtype,sep="\t",col.names=F, row.names=F, quote=F)
        
        # intersect with truth
        # s=paste(para$bedtools," intersect -a ", toolBEDtype, " -b ",trueBEDtype," -f 0.5 -r -c > ",toolTrueMarktype,sep="")
        # system(s)
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
  toolDir=paste(para$tmpDir,"/",para$dellyDir,sep="")
  # trueDir=paste(para$tmpDir,"/",para$trueDir,sep="")
  
  for(sample in Sample){
    print(sample)
    
    toolVCFtmp=list.files(path=paste(sampleDir,"/",sample,sep=""),pattern=paste(tool,".vcf",sep=""))
    if(length(toolVCFtmp)==0){
      stop(paste("No VCF file available for tool=",tool,", sample=",sample,sep=""))
    }else if(length(toolVCFtmp)>1){
      stop(paste("Multiple VCF files available for tool=",tool,", sample=",sample,sep=""))
    }
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
      
      # genotype: GT
      GT=sapply(tmpValue,function(x) return(x[1]))
      
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
      
      df=data.frame(chr=invcf$CHROM,startPos=startPos,endPos=endPos,name=SVtype,score=GQ, GT=GT,stringsAsFactors=FALSE)
      df$startPos=as.integer(df$startPos)
      df$endPos=as.integer(df$endPos)
      
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
      
      # filter the data when end exceeding chromsome size or start smaller than 0
      tmpChrSize=as.numeric(para$refLen[as.character(df$chr)])
      correctInd=which(df$startPos>tmpChrSize | df$startPos<=0)
      if(length(correctInd)>0){
        warning(paste(length(correctInd), " of the ",tool," SVs in sample ",sample," have start position larger than chromosome size OR smaller than 0, remove these SVs",sep=""))
        df=df[-correctInd,]
        tmpChrSize=tmpChrSize[-correctInd]
      }
      correctInd=which(df$endPos>tmpChrSize)
      if(length(correctInd)>0){
        warning(paste(length(correctInd)," of the ",tool," SVs in sample ",sample," have end position larger than chromosome size, modify them",sep=""))
        df[correctInd,"endPos"]=tmpChrSize[correctInd]
      }
      
      # write the original bed 
      for(type in para$Type){
        dfTmp=df[df$name==type,]
        fileTmp=paste(toolBEDpre,"_",type,"_ori.bed",sep="")
        write.table(dfTmp,file=fileTmp,sep="\t",col.names=F,row.names=F,quote=F)
      }
      
      ## filter GT, only keep 0/1 or 1/1, remove 0/0 or ./.
      if(para$filterDellyGenotype==TRUE){
        keepInd=which(df$GT%in%para$keepGenotype)
        df=df[keepInd,]
      }
      
      # filter the N region
      if(para$filterNregion==TRUE){
        write.table(data.frame(df[,1],as.integer(df[,2]-1),as.integer(df[,3])),file=toolBED1,sep="\t",col.names=F, row.names=F, quote=F)
        Ncount=as.numeric(system(paste(para$bedtools," getfasta -tab -fi ",para$reference," -bed ",toolBED1," | awk -F \"\\t\" '{print $2}' | awk -F \"N|n\" '{print NF-1}'",sep=""),intern=T))
        if(length(Ncount)!=nrow(df)){
          stop(paste("File ",toolBED1, " is wrong input for getfasta tool for sample=",sample," tool=delly",sep=""))
        }
        Npercent=Ncount/(df$endPos-df$startPos+1)
        keepInd=which(Npercent<para$Nrate)
        df=df[keepInd,]
      }
      
      # split the data by type
      for(type in para$Type){
        
        toolBEDtype=paste(toolBEDpre,"_",type,".bed",sep="")
        # trueBEDtype=paste(trueDir,"/",sample,"_",type,".bed",sep="")
        # toolTrueMarktype=paste(toolDir,"/",sample,"_",type,"_mark.txt",sep="")
        
        dfSub=df[df$name==type,]
        write.table(dfSub, file=toolBEDtype,sep="\t",col.names=F, row.names=F, quote=F)
        
        # intersect with truth
        # s=paste(para$bedtools," intersect -a ", toolBEDtype, " -b ",trueBEDtype," -f 0.5 -r -c > ",toolTrueMarktype,sep="")
        # system(s)
      }
    }else{
      print(paste("Sample has been prepared, skip",sep=""))
    }
  }  # end for(sample in Sample)
} # end prepare_delly function


#### prepare_general
# to prepare general VCF input: filtering, split by SV length bin and SV type
prepare_general <- function(para,sampleDir,Sample, toolIndex){
  
  ## for test, parameter
  # sampleDir=para$trainDataDir
  # Sample=para$trainSample
  # toolIndex=1
  
  tool=names(para$scoreID)[[toolIndex]]
  toolDir=paste(para$tmpDir,"/",para$toolDirPre,"_",tool,sep="")
  
  for(sample in Sample){
    print(sample)
    
    toolVCFtmp=list.files(path=paste(sampleDir,"/",sample,sep=""),pattern=paste(tool,".vcf",sep=""))
    if(length(toolVCFtmp)==0){
      stop(paste("No VCF file available for tool=",tool,", sample=",sample,sep=""))
    }else if(length(toolVCFtmp)>1){
      stop(paste("Multiple VCF files available for tool=",tool,", sample=",sample,sep=""))
    }
    toolVCF=paste(sampleDir,"/",sample,"/",toolVCFtmp,sep="")
    toolBED1=paste(toolDir,"/",sample,".1.bed",sep="")
    toolBEDpre=paste(toolDir,"/",sample,sep="")
    
    # if the file has been processed, no need to processed again
    if(!file.exists(paste(toolBEDpre,"_",para$Type[1],".bed",sep=""))){
      
      # read in vcf file
      invcf=read.table(toolVCF,sep="\t",as.is=T,comment.char="#",header=F)
      colnames(invcf)[1:8]=c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO")
      
      tmpINFO=strsplit(invcf$INFO,split=";")
      
      # SVtype
      SVTYPEtmp=sapply(tmpINFO,function(x) return(x[startsWith(x,"SVTYPE=")]))
      SVtype=gsub("SVTYPE=","",as.character(SVTYPEtmp))
      
      # END
      SVENDtmp=sapply(tmpINFO,function(x) return(x[startsWith(x,"END=")]))
      SVend=as.numeric(gsub("END=","",as.character(SVENDtmp)))
      
      # SVlen
      SVLENtmp=sapply(tmpINFO,function(x) return(x[startsWith(x,"SVLEN=")]))
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
      
      # get score
      score=matrix(0,nrow=nrow(invcf),ncol=nrow(para$scoreID[[toolIndex]]))
      for(scoreIDind in 1:ncol(score)){
        if(para$scoreID[[toolIndex]][scoreIDind,"col"]=="6"){
          ## QUAL column
          score[,scoreIDind]=as.numeric(invcf$QUAL)
        }else if(para$scoreID[[toolIndex]][scoreIDind,"col"]=="8"){
          ## INFO column
          scoreTmp=sapply(tmpINFO,function(x) return(x[startsWith(x,paste(para$scoreID[[toolIndex]][scoreIDind,"ID"],"=",sep=""))]))
          score[,scoreIDind]=as.numeric(gsub(paste(para$scoreID[[toolIndex]][scoreIDind,"ID"],"=",sep=""),"",as.character(scoreTmp)))
        }else if(para$scoreID[[toolIndex]][scoreIDind,"col"]=="10"){
          ## FORMAT sample 
          if(ncol(invcf)<10){
            stop(paste("VCF file doesn't have col 10 available for sample=",sample," tool=",tool,sep=""))
          }
          formatTmp=strsplit(invcf[,9],split=":")
          itemIndex=sapply(formatTmp,function(x) return(which(x==para$scoreID[[toolIndex]][scoreIDind,"ID"])))
          formatSampleTmp=strsplit(invcf[,10],split=":")
          score[,scoreIDind]=as.numeric(unlist(sapply(1:length(formatSampleTmp),function(i) return(formatSampleTmp[[i]][itemIndex[i]]))))
        }else{
          stop("scoreID only accept column 6, 8 and 10")
        }
      }
      # modify all the NA values into 1
      score[is.na(score)]=1
      
      # check GT & form data frame
      if(toupper(tool)=="DELLY" & para$filterDellyGenotype==TRUE){
        ## FORMAT sample 
        if(ncol(invcf)<10){
          stop(paste("VCF file doesn't have col 10 genotype for sample=",sample," tool=",tool,sep=""))
        }
        formatTmp=strsplit(invcf[,9],split=":")
        itemIndex=sapply(formatTmp,function(x) return(which(x=="GT")))
        formatSampleTmp=strsplit(invcf[,10],split=":")
        GT=as.character(unlist(sapply(1:length(formatSampleTmp),function(i) return(formatSampleTmp[[i]][itemIndex[i]]))))
        df=data.frame(chr=invcf$CHROM,startPos=startPos,endPos=endPos,name=SVtype,score=score,GT=GT,stringsAsFactors=FALSE)
      }else{
        df=data.frame(chr=invcf$CHROM,startPos=startPos,endPos=endPos,name=SVtype,score=score,stringsAsFactors=FALSE)
      }
      
      df$startPos=as.integer(df$startPos)
      df$endPos=as.integer(df$endPos)
      
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
      
      # filter the data when end exceeding chromsome size or start smaller than 0
      tmpChrSize=as.numeric(para$refLen[as.character(df$chr)])
      correctInd=which(df$startPos>tmpChrSize | df$startPos<=0)
      if(length(correctInd)>0){
        warning(paste(length(correctInd), " of the ",tool," SVs in sample ",sample," have start position larger than chromosome size OR smaller than 0, remove these SVs",sep=""))
        df=df[-correctInd,]
        tmpChrSize=tmpChrSize[-correctInd]
      }
      correctInd=which(df$endPos>tmpChrSize)
      if(length(correctInd)>0){
        warning(paste(length(correctInd)," of the ",tool," SVs in sample ", sample," have end position larger than chromosome size, modify them",sep=""))
        df[correctInd,"endPos"]=tmpChrSize[correctInd]
      }
      
      # write the original bed 
      for(type in para$Type){
        dfTmp=df[df$name==type,]
        fileTmp=paste(toolBEDpre,"_",type,"_ori.bed",sep="")
        write.table(dfTmp,file=fileTmp,sep="\t",col.names=F,row.names=F,quote=F)
      }
      
      ## filter GT, only keep 0/1 or 1/1, remove 0/0 or ./.
      if(toupper(tool)=="DELLY" & para$filterDellyGenotype==TRUE){
        keepInd=which(df$GT%in%para$keepGenotype)
        df=df[keepInd,]
      }
      
      # filter the N region
      if(para$filterNregion==TRUE){
        write.table(data.frame(df[,1],as.integer(df[,2]-1),as.integer(df[,3])),file=toolBED1,sep="\t",col.names=F, row.names=F, quote=F)
        Ncount=as.numeric(system(paste(para$bedtools," getfasta -tab -fi ",para$reference," -bed ",toolBED1," | awk -F \"\\t\" '{print $2}' | awk -F \"N|n\" '{print NF-1}'",sep=""),intern=T))
        if(length(Ncount)!=nrow(df)){
          stop(paste("File ",toolBED1, " is wrong input for getfasta tool for sample=",sample," tool=",tool,sep=""))
        }
        Npercent=Ncount/(df$endPos-df$startPos+1)
        keepInd=which(Npercent<para$Nrate)
        df=df[keepInd,]
      }
      
      # split the data by type
      for(type in para$Type){
        
        toolBEDtype=paste(toolBEDpre,"_",type,".bed",sep="")
        #trueBEDtype=paste(trueDir,"/",sample,"_",type,".bed",sep="")
        #toolTrueMarktype=paste(toolDir,"/",sample,"_",type,"_mark.txt",sep="")
        
        dfSub=df[df$name==type,]
        write.table(dfSub, file=toolBEDtype,sep="\t",col.names=F, row.names=F, quote=F)
        
        # intersect with truth
        #s=paste(para$bedtools," intersect -a ", toolBEDtype, " -b ",trueBEDtype," -f 0.5 -r -c > ",toolTrueMarktype,sep="")
        #system(s)
      }
    }else{
      print(paste("Sample has been processed, skip",sep=""))
    }
  }  # end for(sample in Sample)
} # end prepare_general function


#### prepareData
# to prepare the data for  
prepareData <- function(para){
  
  system(paste("mkdir -p ",para$tmpDir,sep=""))
  system(paste("mkdir -p ",para$tmpDir,"/",para$trueDir,sep=""))
  
  for(toolIndex in 1:length(para$scoreID)){
    system(paste("mkdir -p ",para$tmpDir,"/",para$toolDirPre,"_",names(para$scoreID)[[toolIndex]],sep=""))
  }
  
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
    
    # tool
    for(toolIndex in 1:length(para$scoreID)){
      print(paste("[",Sys.time(),"] prepare ",names(para$scoreID)[[toolIndex]]," training data",sep=""))
      prepare_general(para=para,sampleDir=sampleDir,Sample=Sample,toolIndex=toolIndex)
    }
    
  }
  
  ### for testing
  print(paste("[",Sys.time(),"] Prepare testing data",sep=""))
  sampleDir=para$testData
  Sample=para$testSample
  
  # truth
  # in case testing data have truth available, for future testing results evaluation
  print(paste("[",Sys.time(),"] prepare true testing data",sep=""))
  prepare_true(para=para,sampleDir=sampleDir,Sample=Sample)
  
  # tool
  for(toolIndex in 1:length(para$scoreID)){
    print(paste("[",Sys.time(),"] prepare ",para$scoreID[[toolIndex]][["toolName"]]," testing data",sep=""))
    prepare_general(para=para,sampleDir=sampleDir,Sample=Sample,toolIndex=toolIndex)
  }
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
    
    scoreNum=ncol(SVscore)
    
    ##### sort the data by left position first, then right position
    # sort by right position
    sortInd1=order(SVpos[,2])
    SVpos=SVpos[sortInd1,]
    SVscore=matrix(SVscore[sortInd1,],ncol=scoreNum)
    
    # sort by left position
    sortInd2=order(SVpos[,1])
    SVpos=SVpos[sortInd2,]
    SVscore=matrix(SVscore[sortInd2,],ncol=scoreNum)
    
    ## skip remove duplication steps, 
    #since duplciated SVs will have the same mean value for pieces
    
    ##### SV to piece
    # initialize
    pieceScore=c()
    piecePos=c()
    currentPiecePos=c(SVpos[1,1],0)   # current piece position, c(left, right)
    posHP=matrix(c(SVpos[1,2],1),nrow=1,ncol=2) # col: right position, index
    
    # iteration 
    for(i in 2:nrow(SVpos)){
      if(SVpos[i,1]!=SVpos[i-1,1]){
        while(SVpos[i,1]>min(posHP[,1])){ # if posHP is empty, then min(posHP) returns Inf, then the expression is false
          minInd=which.min(posHP[,1])
          if(length(piecePos)==0){
            preRight=-1  # a minus value
          }else{
            preRight=piecePos[nrow(piecePos),2]
          }
          if(posHP[minInd,1]!=preRight){   # check here later
            currentPieceScore=apply(matrix(SVscore[posHP[,2],],ncol=scoreNum),2,mean,na.rm=T)  # use the mean score value, change later for other method
            pieceScore=rbind(pieceScore,currentPieceScore)
            currentPiecePos[2]=posHP[minInd,1]
            piecePos=rbind(piecePos,currentPiecePos)
            currentPiecePos=c(posHP[minInd,1]+1,0)
          }
          posHP=matrix(posHP[-minInd,],ncol=2)
        }
        if(SVpos[i,1]>currentPiecePos[1]){
          currentPiecePos[2]=SVpos[i,1]-1
          piecePos=rbind(piecePos,currentPiecePos)
          currentPieceScore=apply(matrix(SVscore[posHP[,2],],ncol=scoreNum),2,mean,na.rm=T) 
          pieceScore=rbind(pieceScore,currentPieceScore)
          currentPiecePos=c(SVpos[i,1],0)
        }
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
        currentPieceScore=apply(matrix(SVscore[posHP[,2],],ncol=scoreNum),2,mean,na.rm=T)  # use the mean score value, change later for other method
        pieceScore=rbind(pieceScore,currentPieceScore)
        currentPiecePos[2]=posHP[minInd,1]
        piecePos=rbind(piecePos,currentPiecePos)
        currentPiecePos=c(posHP[minInd,1]+1,0)
      }
      posHP=matrix(posHP[-minInd,],ncol=2)
    }
    
    # for those gap region
    naInd=apply(is.na(pieceScore),1,all)
    piecePos=piecePos[!naInd,]
    pieceScore=pieceScore[!naInd,]
    
    rownames(piecePos)=NULL
    rownames(pieceScore)=NULL
  }
  return(list(piecePos=piecePos,pieceScore=pieceScore))
}

### mergeSV
# to merge the SVs from different bins together
mergeSV <- function(inPos, inScore, SVgapDis=1){
  # inPos -input matrix, these data should be within the same chromosome
  #   inPos[,1] for left position
  #   inPos[,2] for right position
  # inScore - score vector corresponding with inPos, output score from piece2SV score
  # SVgapDis - gap distance to merge the SVs
  # output-outPos - left and right position for merged SV
  # output-outScore - start:end:pieceScore, comma separated if multiple pieces or multiple SVs
  
  if(nrow(inPos)!=length(inScore)){
    stop("Error: The numbers of input piece are not equal for inPos and inScore")
  }
  
  if(length(inPos)==0){
    warning("Zero pieces to be merge into SV")
    outPos=matrix(0,nrow=0,ncol=2)
    outScore=c()
  }else if(nrow(inPos)==1){  # only one line
    outPos=matrix(as.integer(inPos),nrow=nrow(inPos),ncol=ncol(inPos))
    outPos=inPos
    outScore=inScore
  }else{
    inPos=matrix(as.integer(inPos),nrow=nrow(inPos),ncol=ncol(inPos))
    N=length(inScore)
    
    ## sort the input data
    # by right position
    ind1=order(inPos[,2])
    inPos=inPos[ind1,]
    inScore=inScore[ind1]
    
    # by left position
    ind2=order(inPos[,1])
    inPos=inPos[ind2,]
    inScore=inScore[ind2]
    
    #print(cbind(inPos,inScore))
    
    ## initialization
    outPos=c()
    outScore=c()
    
    currentSVpos=c(inPos[1,1],0)
    currentSVscore=inScore[1]
    maxRight=inPos[1,2]
    
    ## iteration
    for(i in 2:N){
      if(inPos[i,1]-maxRight > SVgapDis){  ## not merge, creat a new SV
        # add the current SV
        currentSVpos[2]=maxRight
        outPos=rbind(outPos,currentSVpos)
        outScore=c(outScore,currentSVscore)
        
        # creat a new SV
        currentSVpos=c(inPos[i,1],0)
        currentSVscore=inScore[i]
        maxRight=inPos[i,2]
        
      }else{  # merge
        maxRight=max(maxRight,inPos[i,2])
        currentSVscore=paste(currentSVscore,inScore[i],sep=",")
      }
    }
    
    ## end, complete the last SV
    currentSVpos[2]=maxRight
    outPos=rbind(outPos,currentSVpos)
    outScore=c(outScore,currentSVscore)
  }
  
  rownames(outPos)=NULL
  return(list(outPos=outPos,outScore=outScore))
}


#### formatPieceSVsample
# a function to be called by formatPieceSV to prepare the data
formatPieceSVsample <- function(para,sampleDir,Sample){
  
  pieceDir=paste(para$tmpDir,"/",para$pieceDir,sep="")
  
  #[[type]][[sample]][[as.character(as.integer(binlen))]][[as.character(chr)]]
  for(type in para$Type){
    
    if(type=="DEL"){
      binLen=para$DELbin
    }else if(type=="DUP"){
      binLen=para$DUPbin
    }else{
      stop(paste("Type ",type," is not supported yet!"))
    }
    binLenAll=c(binLen,Inf)
    
    for(sample in Sample){
      
      trueBEDtype=paste(para$tmpDir,"/",para$trueDir,"/",sample,"_",type,".bed",sep="")
      
      tmp=list.files(path=pieceDir,pattern=paste(type,"_",sample,sep=""))
      if(length(tmp)!=0){
        print(paste("type=",type," sample=",sample," has been pieced, skip",sep=""))
        next
      }else{
        print(paste("type=",type," sample=",sample,sep=""))
        
        ## read in and combine the data 
        SVposAll=matrix(0,nrow=0,ncol=3)
        SVscoreAll=matrix(0,nrow=0,ncol=length(para$allToolName))
        
        scoreIndex=0
        for(toolIndex in 1:length(para$scoreID)){
          toolBEDtype=paste(para$tmpDir,"/",para$toolDirPre,"_",names(para$scoreID)[toolIndex],"/",sample,"_",type,".bed",sep="")
          scoreNum=nrow(para$scoreID[[toolIndex]])
          if(file.info(toolBEDtype)$size!=0){
            res=read.table(toolBEDtype,sep="\t",header=F,as.is=T)
            SVposAll=rbind(SVposAll,res[,1:3])
            tmp=matrix(NA,nrow=nrow(res),ncol=length(para$allToolName))  # default is NA
            for(sIndex in 1:scoreNum){
              tmp[,scoreIndex+sIndex]=as.numeric(res[,4+sIndex])
            }
            # TBD add score modification: -log2 or some other operation
            SVscoreAll=rbind(SVscoreAll,tmp)
          }
          scoreIndex=scoreIndex+scoreNum
        }
        
        # modify CNVnator to avaoid some extreme value
        # ind1=res2[,5]<0
        # ind2=res2[,5]==0
        # res2[ind1,5]=0
        # res2[ind2,5]=-1
        
        SVlenAll=as.integer(SVposAll[,3]-SVposAll[,2])
        
        for(binlenInd in 1:length(binLen)){
          binlen=binLen[binlenInd]
          
          for(chr in para$Chr){
            
            pieceFileXY=paste(pieceDir,"/",type,"_",sample,"_",as.character(as.integer(binlen)),"_",chr,"_xy.txt",sep="")
            selectInd=which(SVlenAll>=binLenAll[binlenInd] & SVlenAll<binLenAll[binlenInd+1] & SVposAll[,1]==chr)
            if(length(selectInd)>0){
              pieceOut=SV2piece(SVpos=as.matrix(SVposAll[selectInd,2:3]),SVscore=matrix(SVscoreAll[selectInd,],ncol=ncol(SVscoreAll)))
              pieceOut$pieceScore[is.na(pieceOut$pieceScore)]=0  # set NA to be 0
              out=data.frame(chr,pieceOut$piecePos,type,pieceOut$pieceScore,stringsAsFactors=FALSE)
              out[,2]=as.integer(out[,2])
              out[,3]=as.integer(out[,3])
              
              # intersect with truth
              if(para$markPieceRate==0 & file.exists(trueBEDtype)){ 
                # at least 1bp overlap is marked as y=1
                pieceFile=paste(pieceDir,"/",type,"_",sample,"_",as.character(as.integer(binlen)),"_",chr,".bed",sep="")
                pieceFileMark=paste(pieceDir,"/",type,"_",sample,"_",as.character(as.integer(binlen)),"_",chr,"_mark.txt",sep="")
                
                write.table(out,file=pieceFile,sep="\t",col.names=F,row.names=F,quote=F)
                s=paste(para$bedtools," intersect -a ", pieceFile, " -b ",trueBEDtype," -c > ",pieceFileMark,sep="")
                system(s)
                
                out=read.table(pieceFileMark,header=F,sep="\t",as.is=T) # chr, start, end, type, bd, CNVnator, delly, true
                out[,2]=as.integer(out[,2])
                out[,3]=as.integer(out[,3])
                
                # correct the values for ML input
                # modify y
                yIndex=4+length(para$allToolName)+1
                out[out[,yIndex]>1,yIndex]=1  # change here is multiple scores for one caller
                
                # remove tmp files
                system(paste("rm ",pieceFile,sep=""))
                system(paste("rm ",pieceFileMark,sep=""))
                
              }else if(para$markPieceRate>0 & file.exists(trueBEDtype)){
                # at least markPiecerRate of the piece region has to overlap with truth to mark as y=1
                
                # prepare nonOverlap true file
                trueNonOverlap=paste(para$tmpDir,"/",para$trueDir,"/",sample,"_",type,"_nonOverlap.bed",sep="")
                if(file.exists(trueNonOverlap)==F){
                  
                  if(file.info(trueBEDtype)$size==0){
                    # empty file
                    write.table(c(),file=trueNonOverlap,sep="\t",col.names=F,row.names=F,quote=F)
                  }else{
                    inTrue=read.table(trueBEDtype,header=F,sep="\t",as.is=T)
                    
                    outTrue=c()
                    for(chrTrue in para$Chr){
                      subInd=which(inTrue[,1]==chrTrue)
                      if(length(subInd)==0){
                        next
                      }
                      mergeTrue=mergeSV(inPos=matrix(as.integer(unlist(inTrue[subInd,2:3])),ncol=2),inScore=rep(1,times=length(subInd)),SVgapDis=1)
                      outTrue=rbind(outTrue,data.frame(chr=chrTrue,startPos=mergeTrue$outPos[,1],endPos=mergeTrue$outPos[,2],type=type,stringsAsFactors=FALSE))
                    }
                    
                    write.table(outTrue,file=trueNonOverlap,col.names=F,row.names=F,quote=F,sep="\t")
                  }
                }
                
                # prepare piece file
                pieceFile=paste(pieceDir,"/",type,"_",sample,"_",as.character(as.integer(binlen)),"_",chr,".bed",sep="")
                write.table(out,file=pieceFile,sep="\t",col.names=F,row.names=F,quote=F)
                
                # intersect
                pieceFileMark=paste(pieceDir,"/",type,"_",sample,"_",as.character(as.integer(binlen)),"_",chr,"_mark.txt",sep="")
                s=paste(para$bedtools," intersect -a ", pieceFile, " -b ",trueNonOverlap," -wao > ",pieceFileMark,sep="")
                system(s)
                
                interFile=read.table(pieceFileMark,header=F,sep="\t",as.is=T) 
                # chr, start, end, type, score matrix, trueChr, trueStart, trueEnd, trueType, overlapBP
                interFile[,2]=as.integer(interFile[,2])
                interFile[,3]=as.integer(interFile[,3])
                
                pieceMark=apply(interFile[,2:3],1,paste,collapse="_")
                interFileSplit=split(interFile,f=pieceMark)
                overlapLen=sapply(interFileSplit,function(x) sum(x[,ncol(interFile)]))
                pieceFullLen=sapply(interFileSplit,function(x) as.integer(x[1,3]-x[1,2]+1) )
                selectedInd=which(overlapLen/pieceFullLen >= para$markPieceRate)
                
                outTmp=data.frame(t(sapply(interFileSplit,function(x) return(x[1,1:(4+length(para$allToolName))]))))
                out=data.frame(chr=chr,startPos=as.integer(outTmp[,2]),endPos=as.integer(outTmp[,3]),type=type)
                for(sIndex in 1:length(para$allToolName)){
                  out=cbind(out,as.numeric(outTmp[,4+sIndex]))
                }
                colnames(out)=NULL
                
                out=cbind(out,0)
                yInd=4+length(para$allToolName)+1
                
                if(length(selectedInd)>0){
                  out[selectedInd,yInd]=1
                }
                
                # sort by position
                out=out[order(out[,2]),]
                
                # remove tmp files
                system(paste("rm ",pieceFile,sep=""))
                system(paste("rm ",pieceFileMark,sep=""))
                
              }else{
                # y=0 if no known truth
                out=cbind(out,0)  
              }
              
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
  
  system(paste("mkdir -p ",para$tmpDir,sep=""))
  system(paste("mkdir -p ",para$tmpDir,"/",para$pieceDir,sep=""))
  
  ##### for training
  if(para$knownModel==TRUE){
    print(paste("[",Sys.time(),"] Known model, skip breaking training data into pieces",sep=""))
  }else{
    print(paste("[",Sys.time(),"] Breaking training data into pieces",sep=""))
    Sample=para$trainSample
    formatPieceSVsample(para=para,Sample=Sample)
  }
  
  ##### for testing
  print(paste("[",Sys.time(),"] Breaking testing data into pieces",sep=""))
  Sample=para$testSample
  formatPieceSVsample(para=para,Sample=Sample)
}


#### applyTrain
# a function to apply each ML model training individually
applyTrain <- function(para,x,y,type,binlen,chr){
  # x - training x matrix, piece x features
  # y - training y, vector with length to be number of pieces
  # output: outModel - return the model list
  
  ## check input format
  if(nrow(x)!=length(y)){
    stop(paste("nrow(x) != length(y)",sep=""))
  }
  
  ## model initialization
  outModel=list()
  for(modelName in para$MLmethod){
    outModel[[modelName]]=NA
  }
  
  ## check y label
  if(length(y)==0){
    # no supporting training data
    print(paste("No supporting training data for type=",type,", binlen=",as.character(as.integer(binlen)),", chr=",chr,sep=""))
    return(outModel)
  }else if(length(unique(y))==1){
    # only one y label
    print(paste("Only y=",unique(y)," for type=",type,", binlen=",as.character(as.integer(binlen)),", chr=",chr,sep=""))
    return(outModel)
  }else if (sum(y==0)==1){ # otherwise 'sample' function will get error
    print(paste("Not enough number of y=0 cases for type=", type,", binlen=",as.character(as.integer(binlen)),", chr=",chr,sep=""))
    return(outModel)
  }else if(sum(y==1)==1){ # # otherwise 'sample' function will get error
    print(paste("Not enough number of y=1 cases for type=", type,", binlen=",as.character(as.integer(binlen)),", chr=",chr,sep=""))
    return(outModel)
  }
  
  # balance y=0 and y=1 to make num(y=1)/num(y=0) = para$pieceBalanceRate
  y0Ind=which(y==0)
  y1Ind=which(y==1)
  if(length(y0Ind)*para$pieceBalanceRate>length(y1Ind)){
    # increase y1 number 
    y1Ind=c(y1Ind,sample(y1Ind,size=max(0,round(length(y0Ind)*para$pieceBalanceRate-length(y1Ind))),replace=T)) # at least one copy
  }else{
    # increase y0 number
    y0Ind=c(y0Ind,sample(y0Ind,size=max(0,round(length(y1Ind)/para$pieceBalanceRate-length(y0Ind))),replace=T)) # at least one copy
  }
  
  x=x[c(y0Ind,y1Ind),]
  y=y[c(y0Ind,y1Ind)]
  colnames(x)=NULL
  x=as.matrix(x)
  
  ## to avoid constant features
  ## TBD: check this part later, either assign some random values OR remove the constant features
  xx=x
  for(ylabel in c(0,1)){
    yInd=which(y==ylabel)
    sd0Ind=which(apply(matrix(x[yInd,],nrow=length(yInd),ncol=ncol(x)),2,sd)<1E-10) # very small sd
    if(length(sd0Ind)>0){
      #print(c(ylabel,sd0Ind))
      for(j in sd0Ind){
        xx[yInd,j]=xx[yInd,j]+ rnorm(length(yInd),mean=0,sd=(mean(x[yInd,j])/10+0.01))   # add some noise
      }
    }
  }
  
  # model 1: SVM radial
  if("SVMradial"%in%para$MLmethod){
    outModel[["SVMradial"]] = svm(x=xx, y=as.factor(y), kernel="radial", cost = 100, probability=TRUE)  # default cost
  }
  
  # model 2: SVM polynomial
  if("SVMpolynomial"%in%para$MLmethod){
    outModel[["SVMpolynomial"]] = svm(x=xx, y=as.factor(y), kernel="polynomial", cost = 100, probability=TRUE)  # default cost
  }
  
  # model 3: RF
  if("RF"%in%para$MLmethod){
    outModel[["RF"]]=randomForest(x=x,y=as.factor(y))
  }
  
  # model 4: NN
  if("NN"%in%para$MLmethod){
    yy=cbind(y,1-y)
    outModel[["NN"]] = nnet(x, yy,size = 2, rang = 0.1, decay = 5e-4, maxit = 200)  # adjust parameters later
  }
  
  # model 5: LDA
  if("LDA"%in%para$MLmethod){
    outModel[["LDA"]] = lda(x=xx,grouping=as.factor(y))
  }
  
  # model 6: adaboost
  if("adaboost"%in%para$MLmethod){
    prepData=data.frame(x,Y=y)
    prepData$Y=factor(prepData$Y)
    outModel[["adaboost"]] = adaboost(Y ~., prepData,nIter=100)
  }
  
  return(outModel)
}


#### trainModel
# a function to train the SVlearning model
trainModel <- function(para, Sample){
  
  print(paste("[",Sys.time(),"] Train the model",sep=""))
  
  #Sample=para$trainSample
  
  model=list()
  #[[type]][[as.character(as.integer(binlen))]][[as.character(chr)]][[modelName]]
  
  for(type in para$Type){
    
    model[[type]]=list()
    
    if(type=="DEL"){
      binLen=para$DELbin
    }else if(type=="DUP"){
      binLen=para$DUPbin
    }else{
      stop(paste("Type ",type," not supporting yet!"))
    }
    binLenAll=c(binLen,Inf)
    
    for(binlen in binLen){
      print(paste("type=",type,", binlen=",as.character(as.integer(binlen)),sep=""))
      model[[type]][[as.character(as.integer(binlen))]]=list()
      
      if(para$binByChr==T){
        # train the model by chromosome
        
        for(chr in para$Chr){
          ### prepare training data
          x=c()
          y=c()
          
          for(sample in Sample){
            pieceFile=paste(para$tmpDir,"/",para$pieceDir,"/",type,"_",sample,"_",as.character(as.integer(binlen)),"_",chr,"_xy.txt",sep="")
            if(file.info(pieceFile)$size!=0){
              res=read.table(pieceFile,header=F,sep="\t",as.is=T)
              x=rbind(x,res[,5:(4+length(para$allToolName))])
              y=c(y,res[,ncol(res)])
            }
          }
          
          ### train the model
          model[[type]][[as.character(as.integer(binlen))]][[as.character(chr)]]=applyTrain(para,x,y,type,binlen,chr)
          
        } # end for(chr in para$Chr)
      }else{
        # pool all the chromosomes together
        
        ### prepare training data
        x=c()
        y=c()
        
        for(chr in para$Chr){
          for(sample in Sample){
            pieceFile=paste(para$tmpDir,"/",para$pieceDir,"/",type,"_",sample,"_",as.character(as.integer(binlen)),"_",chr,"_xy.txt",sep="")
            if(file.info(pieceFile)$size!=0){
              res=read.table(pieceFile,header=F,sep="\t",as.is=T)
              x=rbind(x,res[,5:(4+length(para$allToolName))])
              y=c(y,res[,ncol(res)])
            }
          }
        } 
        
        ### train the model
        model[[type]][[as.character(as.integer(binlen))]][["allChr"]]=applyTrain(para,x,y,type,binlen,chr="allChr")
      }
      
    } # end for(binlen in binLen)
  }  # end for(type in para$Type)
  
  # save the model
  system(paste("mkdir -p ",para$modelDir,sep=""))
  save(model,file=paste(para$modelDir,"/",para$modelFile,sep=""))
}


### piece2SV
# to merge the piece to SV by gapDis
piece2SV <- function(piecePos, pieceScore,gapDis=1){
  # piecePos -input matrix, these data should be sorted and nonoverlapped, within the same chromosome
  #   piecePos[,1] for left position
  #   piecePos[,2] for right position
  # pieceScore - score matrix corresponding with piecePos, cbind(callerScore, MLpredict)
  # gapDis - gap distance to merge the pieces
  # output-SVpos - left and right position for SV
  # output-SVscore - start:end:pieceScore, comma separated if multiple pieces
  
  if(nrow(piecePos)!=nrow(pieceScore)){
    stop("Error: The numbers of input piece are not equal for piecePos and pieceScore")
  }
  
  if(length(piecePos)==0){
    warning("Zero pieces to be merge into SV")
    SVpos=matrix(0,nrow=0,ncol=2)
    SVscore=c()
  }else if(nrow(piecePos)==1){  # only one line
    piecePos=matrix(as.integer(piecePos),nrow=nrow(piecePos),ncol=ncol(piecePos))
    SVpos=piecePos
    SVscore=paste(c(piecePos[1,],pieceScore[1,]),collapse=":")
  }else{
    piecePos=matrix(as.integer(piecePos),nrow=nrow(piecePos),ncol=ncol(piecePos))
    N=nrow(pieceScore)
    
    ## initialization
    SVpos=c()
    SVscore=c()
    
    currentSVpos=c(piecePos[1,1],0)
    currentSVscore=paste(c(piecePos[1,],pieceScore[1,]),collapse=":")
    
    ## iteration
    for(i in 2:N){
      if(piecePos[i,1]-piecePos[i-1,2]>gapDis){
        # complete the current SV
        currentSVpos[2]=piecePos[i-1,2]
        SVpos=rbind(SVpos,currentSVpos)
        SVscore=c(SVscore,currentSVscore)
        # start a new SV
        currentSVpos=c(piecePos[i,1],0)
        currentSVscore=paste(c(piecePos[i,],pieceScore[i,]),collapse=":")
      }else{
        # keep adding to the current SV
        currentSVscore=paste(currentSVscore,",",paste(c(piecePos[i,],pieceScore[i,]),collapse=":"),sep="")
      }
    }
    
    ## end, add the last piece and complete the SV
    currentSVpos[2]=piecePos[N,2]
    SVpos=rbind(SVpos,currentSVpos)
    SVscore=c(SVscore,currentSVscore)
    
  }
  rownames(SVpos)=NULL
  return(list(SVpos=SVpos,SVscore=SVscore))
}


#### modelPredict
# to predict the label based on the model
modelPredict <- function(para, Sample){
  
  print(paste("[",Sys.time(),"] Predict testing data",sep=""))
  
  ### load the reference length file
  refLen=read.table(para$refLenFile,header=F,sep="\t",as.is=T)
  
  #### load the model 
  load(paste(para$modelDir,"/",para$modelFile,sep=""))
  #[[type]][[as.character(as.integer(binlen))]][[as.character(chr)]][[modelName]]
  
  modelNum=length(para$MLmethod)
  
  #### prediction within each bin
  for(sample in Sample){
    
    print(paste("[",Sys.time(),"] sample=",sample,sep=""))
    
    SVall=data.frame(stringsAsFactors=FALSE)
    
    for(type in para$Type){
      
      print(paste("sample=",sample,", type=",type,sep=""))
      
      if(type=="DEL"){
        binLen=para$DELbin
      }else if(type=="DUP"){
        binLen=para$DUPbin
      }else{
        stop(paste("Type ",type," not supporting yet!"))
      }
      binLenAll=c(binLen,Inf)
      
      for(chr in para$Chr){
        
        if(para$binByChr==T){
          chrName=as.character(chr)
        }else{
          chrName="allChr"
        }
        
        # initialization, to record SVpos and SVscore within each chr
        SVpos1=c()
        SVscore1=c()
        
        for(binlen in binLen){
          pieceFileXY=paste(para$tmpDir,"/",para$pieceDir,"/",type,"_",sample,"_",as.character(as.integer(binlen)),"_",as.character(chr),"_xy.txt",sep="")
          
          if(file.exists(pieceFileXY)==FALSE){
            stop(paste("No piece file for sample=",sample,", type=",", binlen=",as.character(as.integer(binlen)),", chr=",as.character(chr),". Please run function formatPieceSV() first.",sep=""))
          }
          
          if(file.info(pieceFileXY)$size==0){ # no testing data
            next
          }
          
          res=read.table(pieceFileXY,header=F,sep="\t",as.is=T)
          
          x=as.matrix(res[,5:(4+length(para$allToolName))])
          colnames(x)=NULL
          y=as.numeric(res[,ncol(res)]) 
          
          ## a matrix to record predict results for each ML model
          y.predict=matrix(0,nrow=nrow(res),ncol=modelNum+1)
          colnames(y.predict)=c(para$MLmethod,"vote")
          
          # model 1: SVM radial
          modelName="SVMradial"
          if(modelName%in%para$MLmethod){
            currentModel=model[[type]][[as.character(as.integer(binlen))]][[chrName]][[modelName]]
            if(is.null(currentModel) | all(is.na(currentModel))){
              print(paste("No model for sample=",sample,", type=",type,", binlen=",binlen,", chr=",chr,", model=",modelName,", predict all as 0",sep=""))
            }else{
              y.predict[,modelName] = as.numeric(predict(currentModel, x, decision.values=T))-1
            }
          }
          
          # model 2: SVM polynomial
          modelName="SVMpolynomial"
          if(modelName%in%para$MLmethod){
            currentModel=model[[type]][[as.character(as.integer(binlen))]][[chrName]][[modelName]]
            if(is.null(currentModel) | all(is.na(currentModel))){
              print(paste("No model for sample=",sample,", type=",type,", binlen=",binlen,", chr=",chr,", model=",modelName,", predict all as 0",sep=""))
            }else{
              y.predict[,modelName] = as.numeric(predict(currentModel, x, decision.values=T))-1
            }
          }
          
          # model 3: RF
          modelName="RF"
          if(modelName%in%para$MLmethod){
            currentModel=model[[type]][[as.character(as.integer(binlen))]][[chrName]][[modelName]]
            if(is.null(currentModel) | all(is.na(currentModel))){
              print(paste("No model for sample=",sample,", type=",type,", binlen=",binlen,", chr=",chr,", model=",modelName,", predict all as 0",sep=""))
            }else{
              y.predict[,modelName] = as.numeric(predict(currentModel,x))-1
            }
          }
          
          # model 4: NN
          modelName="NN"
          if(modelName%in%para$MLmethod){
            currentModel=model[[type]][[as.character(as.integer(binlen))]][[chrName]][[modelName]]
            if(is.null(currentModel) | all(is.na(currentModel))){
              print(paste("No model for sample=",sample,", type=",type,", binlen=",binlen,", chr=",chr,", model=",modelName,", predict all as 0",sep=""))
            }else{
              preds=predict(currentModel, x)
              y.predict[,modelName] = as.numeric(preds[,1]>preds[,2])
            }
          }
          
          # model 5: LDA
          modelName="LDA"
          if(modelName%in%para$MLmethod){
            currentModel=model[[type]][[as.character(as.integer(binlen))]][[chrName]][[modelName]]
            if(is.null(currentModel) | all(is.na(currentModel))){
              print(paste("No model for sample=",sample,", type=",type,", binlen=",binlen,", chr=",chr,", model=",modelName,", predict all as 0",sep=""))
            }else{
              y.predict[,modelName] = as.numeric(predict(currentModel,x)$class)-1
            }
          }
          
          # model 6: adaboost
          modelName="adaboost"
          if(modelName%in%para$MLmethod){
            currentModel=model[[type]][[as.character(as.integer(binlen))]][[chrName]][[modelName]]
            if(is.null(currentModel) | all(is.na(currentModel))){
              print(paste("No model for sample=",sample,", type=",type,", binlen=",binlen,", chr=",chr,", model=",modelName,", predict all as 0",sep=""))
            }else{
              y.predict[,modelName] = as.numeric(predict(currentModel,data.frame(x))$class)-1
            }
          }
          
          ### vote to get final prediction
          y.predict[,"vote"]=as.numeric(apply(matrix(y.predict[,1:modelNum],ncol=modelNum),1,sum)>=para$vote)
          
          ### merge piece to SV
          selectInd=which(y.predict[,"vote"]==1)
          if(length(selectInd)!=0){
            # piece2SV, output score with format, start:end:callerScore:predict
            tmp=piece2SV(piecePos=as.matrix(res[selectInd,c(2,3)]), 
                         pieceScore=cbind(matrix(x[selectInd,],nrow=length(selectInd),ncol=ncol(x)),
                                          matrix(y.predict[selectInd,1:modelNum],nrow=length(selectInd),ncol=modelNum)),
                         gapDis=para$PgapDis)
            SVpos1=rbind(SVpos1,tmp$SVpos)
            SVscore1=c(SVscore1,tmp$SVscore)
          }
          
        } # end for(binlen in binLen)
        
        if(length(SVscore1)>0){
          if(para$SVgapDis>0){
            # merge the SVs within given sample, type, chr
            tmp=mergeSV(inPos=SVpos1,inScore=SVscore1,SVgapDis=para$SVgapDis)
            SVall=rbind(SVall, data.frame(chr=chr,startPos=tmp$outPos[,1],endPos=tmp$outPos[,2],type=type,info=tmp$outScore,stringsAsFactors=FALSE))
          }else{  # SVgapDis==0
            # don't merge SV from different bins
            SVall=rbind(SVall, data.frame(chr=chr,startPos=SVpos1[,1],endPos=SVpos1[,2],type=type,info=SVscore1,stringsAsFactors=FALSE))
          }
        }
        
      } # end for(chr in para$Chr)
      
      ## any code for sample & type bin, can be added here
      
    } # end for (type in para$Type)
    
    ### prepare for output
    if(nrow(SVall)>0){
      SVall$chr=factor(SVall$chr,levels=para$Chr)
      
      ### combine all into one file, sort
      # by startPos
      ind=order(SVall$startPos)
      SVall=SVall[ind,]
      # by chr
      ind=order(SVall$chr)
      SVall=SVall[ind,]
      
      # get reference
      out=data.frame(SVall$chr,startPos=SVall$startPos,endPos=SVall$startPos+1,stringsAsFactors=FALSE)
      out$startPos=as.integer(out$startPos)
      out$endPos=as.integer(out$endPos)
      system(paste("mkdir -p ",para$tmpDir,"/",para$shortBedDir,sep=""))
      outFile=paste(para$tmpDir,"/",para$shortBedDir,"/",sample,".short.bed",sep="")
      write.table(out,file=outFile,row.names=F,col.names=F,sep="\t",quote=F)
      REF=as.character(system(paste(para$bedtools," getfasta -tab -fi ",para$reference," -bed ",outFile," | awk -F \"\\t\" '{print $2}' ",sep=""),intern=T))
      if(length(REF)!=nrow(SVall)){
        stop("Errors for reference: non-equal number of reference and SV, SV position may exceed chromosome length.")
      }
      ## TBD: some called regions with multiple starting N or end N strings, do some final filtering
      
      # ID name
      ID=paste("SVlearing",1:nrow(SVall),sep="_")
      
      # ALT
      ALT=paste("<",SVall$type,">",sep="")
      
      # INFO
      INFO=paste("SVTYPE=",SVall$type,";SVLEN=",as.character(as.integer(SVall$endPos-SVall$startPos+1)),";END=",as.character(as.integer(SVall$endPos)),";IMPRECISE",sep="")
    }
    
    toolNameString=paste(para$allToolName,collapse=":")
    
    # BED format
    if("BED"%in%para$outFormat){
      system(paste("mkdir -p ",para$outDir,"/",para$bedDir,sep=""))
      outFile=paste(para$outDir,"/",para$bedDir,"/",sample,".bed",sep="")
      
      if(nrow(SVall)==0){
        cat(paste("#chrom","chromStart","chromEnd","name","REF","ALT",
                  paste("PSTART:PEND:",toolNameString,":",paste(para$MLmethod,collapse=":"),sep=""),sep="\t"),file=outFile,sep="\n")
      }else{
        outSVall=data.frame(chrom=SVall$chr, chromStart=SVall$startPos, chromEnd=SVall$endPos, name=ID, 
                            REF=REF, ALT=ALT, INFO=SVall$info, stringsAsFactors=FALSE) 
        outSVall$chromStart=as.integer(outSVall$chromStart)
        outSVall$chromEnd=as.integer(outSVall$chromEnd)
        colnames(outSVall)[7]=paste("PSTART:PEND:",toolNameString,":",paste(para$MLmethod,collapse=":"),sep="")
        
        cat(paste("#",paste(colnames(outSVall),collapse="\t"),sep=""),file=outFile,sep="\n")
        write.table(outSVall,file=outFile,row.names=F,col.names=F,sep="\t",quote=F,append=T)
      }
      
    }  # end  if("BED"%in%para$outFormat)
    
    # VCF format
    if("VCF"%in%para$outFormat){
      
      system(paste("mkdir -p ",para$outDir,"/",para$vcfDir,sep=""))
      outFile=paste(para$outDir,"/",para$vcfDir,"/",sample,".vcf",sep="")
      
      cat("##fileformat=VCF",file=outFile,sep="\n")
      cat(paste("##fileDate=",Sys.Date(),sep=""),file=outFile,sep="\n",append=TRUE)
      cat(paste("##reference=",para$reference,sep=""),file=outFile,sep="\n",append=TRUE)
      
      cat("##FILTER=<ID=LowQual,Description=\"Low Quality\">",file=outFile,sep="\n",append=TRUE)
      cat("##FILTER=<ID=PASS,Description=\"All filters passed\">",file=outFile,sep="\n",append=TRUE)
      
      cat("##ALT=<ID=DEL,Description=\"Deletion\">",file=outFile,sep="\n",append=TRUE)
      cat("##ALT=<ID=DUP,Description=\"Duplication\">",file=outFile,sep="\n",append=TRUE)
      
      cat("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">",file=outFile,sep="\n",append=TRUE)
      cat("##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">",file=outFile,sep="\n",append=TRUE)
      cat("##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">",file=outFile,sep="\n",append=TRUE)
      cat("##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise structural variation\">",file=outFile,sep="\n",append=TRUE)
      
      cat("##FORMAT=<ID=PSTART,Number=1,Type=Integer,Description=\"Start position of the combining pieces\">",file=outFile,sep="\n",append=TRUE)
      cat("##FORMAT=<ID=PEND,Number=1,Type=Integer,Description=\"End position of the combining pieces\">",file=outFile,sep="\n",append=TRUE)
      
      for(toolIndex in 1:length(para$scoreID)){
        for(sInd in 1:nrow(para$scoreID[[toolIndex]])){
          cat(paste("##FORMAT=<ID=",names(para$scoreID)[[toolIndex]],",Number=1,Type=Double,Description=\"",names(para$scoreID)[[toolIndex]]," ",para$scoreID[[toolIndex]][sInd,"ID"],"\">",sep=""),file=outFile,sep="\n",append=TRUE)
        }
      }
      
      for(methodName in para$MLmethod){
        cat(paste("##FORMAT=<ID=",methodName,",Number=1,Type=Integer,Description=\"",methodName," prediction\">",sep=""),file=outFile,sep="\n",append=TRUE)
      }
      
      for(refInd in 1:nrow(refLen)){
        cat(paste("##contig=<ID=",refLen[refInd,1],",len=",refLen[refInd,2],">",sep=""),file=outFile,sep="\n",append=TRUE)
      }
      
      if(nrow(SVall)==0){
        cat(paste("#CHROM","POS","ID","REF","ALT", "QUAL","FILTER","INFO","FORMAT",sample,sep="\t"),file=outFile,sep="\n",append=TRUE)
      }else{
        
        outSVall=data.frame(CHROM=SVall$chr, POS=SVall$startPos, ID=ID, REF=REF, ALT=ALT, QUAL=".",FILTER="PASS", INFO=INFO,
                            FORMAT=paste("PSTART:PEND:",toolNameString,":",paste(para$MLmethod,collapse=":"),sep=""),
                            formatInfo=SVall$info, stringsAsFactors=FALSE)
        colnames(outSVall)[10]=sample
        outSVall$POS=as.integer(outSVall$POS)
        
        cat(paste("#",paste(colnames(outSVall),collapse="\t"),sep=""),file=outFile,sep="\n",append=TRUE)
        write.table(outSVall,file=outFile,append=T,row.names=F,col.names=F,sep="\t",quote=F)
      }
      
    } # if("VCF"%in%para$outFormat)
    
  }  #end for(sample in SampleTest)
  
}



### eva
# to evaluate the prediction performance compared with truth
# for piece level evaluation
eva<-function(yTrue, yPredict){
  # yTrue and yPredict need to be binary value 0,1
  tab=table(factor(yTrue,levels=c(0,1)),factor(yPredict,levels=c(0,1)))
  TP=tab[2,2]
  FP=tab[1,2]
  TN=tab[1,1]
  FN=tab[2,1]
  sen=TP/(TP+FN)
  spe=TN/(TN+FP)
  acc=(TN+TP)/length(yTrue)
  youden=sen+spe-1
  mcc=(as.double(TP*TN) -as.double(FP*FN)) /sqrt(as.double(TP+FP)*as.double(TP+FN)*as.double(TN+FP)*as.double(TN+FN))
  
  out=c(TP,FP,TN,FN,sen,spe,acc,youden,mcc)
  names(out)=c("TP","FP","TN","FN","Sen","Spe","Acc","Youden","MCC")
  return(out)
}


### callerSVevaluation
# to evaluate the caller results compared with true set
# SV level evaluation, precision and recall
callerSVevaluation <- function(para,Sample){
  # Sample=para$testSample
  
  print(paste("[",Sys.time(),"] SV level evaluation",sep=""))
  
  system(paste("mkdir -p ",para$outDir,"/",para$evaDir,sep=""))
  
  ## evaluation initialization
  colName=c("Predict","True","TrueHitPredict","Precision","PredictHitTrue","Recall")
  evaMat0=matrix(0,nrow=length(Sample),ncol=length(colName))
  rownames(evaMat0)=Sample
  colnames(evaMat0)=colName
  
  ## tool caller evaluation
  for(tool in para$evaCaller){
    print(paste("Evaluating ",tool,sep=""))
    
    for(type in para$Type){
      
      evaMat=evaMat0
      
      for(sample in Sample){
        
        ## true file
        trueFile=paste(para$tmpDir,"/",para$trueDir,"/",sample,"_",type,".bed",sep="")
        if(file.exists(trueFile)==F){
          print(paste("No true file for sample=",sample,", type=",type,", skip evaluation. Suggestion: please try prepareData function first.",sep=""))
          next
        }
        
        if(file.info(trueFile)$size==0){
          print(paste("Empty True file for sample=",sample,", type=",type,", skip evaluation.",sep=""))
          next
        }
        
        ## tool file
        if(tool=="SVlearning"){
          prepBed=paste(para$tmpDir,"/",para$SVlearningDir,"/",sample,"_",type,"_ori.bed",sep="")
          if(!file.exists(prepBed)){
            system(paste("mkdir -p ",para$tmpDir,"/",para$SVlearningDir,sep=""))
            # prepare SV learning files
            bedFile=paste(para$outDir,"/",para$bedDir,"/",sample,".bed",sep="")
            vcfFile=paste(para$outDir,"/",para$vcfDir,"/",sample,".vcf",sep="")
            if(file.exists(bedFile)){
              # read in bed to prepare
              res=read.table(bedFile,sep="\t",as.is=T,header=F)
              resSub=res[which(res[,6]==paste("<",type,">",sep="")),]
              if(nrow(resSub)>0){
                outBed=data.frame(resSub[,1:3],type=type,stringsAsFactors=F)
                outBed[,2]=as.integer(outBed[,2])
                outBed[,3]=as.integer(outBed[,3])
                write.table(outBed,file=prepBed,sep="\t",col.names=F,row.names=F,quote=F)
              }
            }else if(file.exists(vcfFile)){
              # read in vcf to prepare
              res=read.table(vcfFile,sep="\t",as.is=T,header=F,comment.char="#")
              resSub=res[which(res[,5]==paste("<",type,">",sep="")),]
              if(nrow(resSub)>0){
                tmp=strsplit(resSub[,8],split=";")
                SVENDtmp=sapply(tmp,function(x) return(x[startsWith(x,"END=")]))
                SVend=as.integer(gsub("END=","",as.character(SVENDtmp)))
                outVCF=data.frame(resSub[,1:2],SVend,type=type,stringsAsFactors=F)
                outVCF[,2]=as.integer(outVCF[,2])
                write.table(outVCF,file=prepBed,sep="\t",col.names=F,row.names=F,quote=F)
              }
            }else{
              stop(paste("No BED or VCF available for SVlearning sample=",sample,", type=",type,sep=""))
            }
          }
          toolDir=paste(para$tmpDir,"/",para$SVlearningDir,sep="")
          
        }else if(tool %in% names(para$scoreID)){
          toolDir=paste(para$tmpDir,"/",para$toolDirPre,"_",tool,sep="")
        }else{
          stop(paste("Wrong tool name ",tool,sep=""))
        }
        
        toolFile=paste(toolDir,"/",sample,"_",type,"_ori.bed",sep="")
        
        if(file.exists(toolFile)==F){
          if(tool!="SVlearning"){
            print(paste("No tool file for sample=",sample,", type=",type,", skip evaluation. Suggestion: please try prepareData function first.",sep=""))
          }else{
            print(paste("No SVlearning BED or VCF file for sample=",sample,", type=",type,", skip evaluation. Suggesting: please do the prediction first.",sep=""))
          }
          next
        }
        
        if(file.info(trueFile)$size==0){
          print(paste("Empty tool file for sample=",sample,", type=",type,", skip evaluation.",sep=""))
          next
        }
        
        ## mark file 1: true hit predict
        markFile1=paste(toolDir,"/",sample,"_",type,"_mark1.txt",sep="")
        
        # intersect with truth
        s=paste(para$bedtools," intersect -a ", toolFile, " -b ",trueFile," -f ",para$evaOverlapRate," -r -wa -c > ",markFile1,sep="")
        system(s)
        
        # read in mark file 1 to evaluate
        res=read.table(markFile1,header=F,sep="\t",as.is=T)
        evaMat[sample,"Predict"]=nrow(res)
        evaMat[sample,"TrueHitPredict"]=sum(res[,ncol(res)]>0)
        evaMat[sample,"Precision"]=evaMat[sample,"TrueHitPredict"]/evaMat[sample,"Predict"]
        
        ## mark file 2: predict hit true
        markFile2=paste(toolDir,"/",sample,"_",type,"_mark2.txt",sep="")
        
        # intersect with truth
        s=paste(para$bedtools," intersect -a ", trueFile, " -b ",toolFile," -f ",para$evaOverlapRate," -r -wa -c > ",markFile2,sep="")
        system(s)
        
        # read in the mark file 2 to evaluate
        res=read.table(markFile2,header=F,sep="\t",as.is=T)
        evaMat[sample,"True"]=nrow(res)
        evaMat[sample,"PredictHitTrue"]=sum(res[,ncol(res)]>0)
        evaMat[sample,"Recall"]=evaMat[sample,"PredictHitTrue"]/evaMat[sample,"True"]
        evaMat[is.na(evaMat)]=0 # in case 0 denominator
      } # end for(sample in Sample)
      
      # write evaMat to file
      evaFile=paste(para$outDir,"/",para$evaDir,"/",tool,"_",type,".evaSV.txt",sep="")
      write.table(data.frame(Sample=Sample,evaMat),file=evaFile,col.names=T,row.names=F,sep="\t",quote=F)
      
    } # end for(type in para$Type)
  }  # end for(tool in para$evaCaller)
  
}



### callerBPevaluation
# to evaluate the caller results compared with true set
# Base pair level evaluation, precision and recall
callerBPevaluation <- function(para,Sample){
  # Sample=para$testSample
  
  print(paste("[",Sys.time(),"] Base pair level evaluation",sep=""))
  
  system(paste("mkdir -p ",para$outDir,"/",para$evaDir,sep=""))
  
  ## evaluation initialization
  colName=c("Predict","True","Overlap","Precision","Recall")
  evaMat0=matrix(0,nrow=length(Sample),ncol=length(colName))
  rownames(evaMat0)=Sample
  colnames(evaMat0)=colName
  
  ## tool caller evaluation
  for(tool in para$evaCaller){
    print(paste("Evaluating ",tool,sep=""))
    
    for(type in para$Type){
      
      evaMat=evaMat0
      
      for(sample in Sample){
        
        ## true file
        trueFile=paste(para$tmpDir,"/",para$trueDir,"/",sample,"_",type,".bed",sep="")
        if(file.exists(trueFile)==F){
          print(paste("No true file for sample=",sample,", type=",type,", skip evaluation. Suggestion: please try prepareData function first.",sep=""))
          next
        }
        
        if(file.info(trueFile)$size==0){
          print(paste("Empty True file for sample=",sample,", type=",type,", skip evaluation.",sep=""))
          next
        }
        
        ## tool file
        if(tool=="SVlearning"){
          prepBed=paste(para$tmpDir,"/",para$SVlearningDir,"/",sample,"_",type,"_ori.bed",sep="")
          if(!file.exists(prepBed)){
            system(paste("mkdir -p ",para$tmpDir,"/",para$SVlearningDir,sep=""))
            # prepare SV learning files
            bedFile=paste(para$outDir,"/",para$bedDir,"/",sample,".bed",sep="")
            vcfFile=paste(para$outDir,"/",para$vcfDir,"/",sample,".vcf",sep="")
            if(file.exists(bedFile)){
              # read in bed to prepare
              res=read.table(bedFile,sep="\t",as.is=T,header=F)
              resSub=res[which(res[,6]==paste("<",type,">",sep="")),]
              if(nrow(resSub)>0){
                outBed=data.frame(resSub[,1:3],type=type,stringsAsFactors=F)
                outBed[,2]=as.integer(outBed[,2])
                outBed[,2]=as.integer(outBed[,3])
                write.table(outBed,file=prepBed,sep="\t",col.names=F,row.names=F,quote=F)
              }
            }else if(file.exists(vcfFile)){
              # read in vcf to prepare
              res=read.table(vcfFile,sep="\t",as.is=T,header=F,comment.char="#")
              resSub=res[which(res[,5]==paste("<",type,">",sep="")),]
              if(nrow(resSub)>0){
                tmp=strsplit(resSub[,8],split=";")
                SVENDtmp=sapply(tmp,function(x) return(x[startsWith(x,"END=")]))
                SVend=as.integer(gsub("END=","",as.character(SVENDtmp)))
                outVCF=data.frame(resSub[,1:2],SVend,type=type,stringsAsFactors=F)
                outVCF[,2]=as.integer(outVCF[,2])
                write.table(outVCF,file=prepBed,sep="\t",col.names=F,row.names=F,quote=F)
              }
            }else{
              stop(paste("No BED or VCF available for SVlearning sample=",sample,", type=",type,sep=""))
            }
          }
          toolDir=paste(para$tmpDir,"/",para$SVlearningDir,sep="")
          
        }else if(tool %in% names(para$scoreID)){
          toolDir=paste(para$tmpDir,"/",para$toolDirPre,"_",tool,sep="")
        }else{
          stop(paste("Wrong tool name ",tool,sep=""))
        }
        
        toolFile=paste(toolDir,"/",sample,"_",type,"_ori.bed",sep="")
        
        if(file.exists(toolFile)==F){
          if(tool!="SVlearning"){
            print(paste("No tool file for sample=",sample,", type=",type,", skip evaluation. Suggestion: please try prepareData function first.",sep=""))
          }else{
            print(paste("No SVlearning BED or VCF file for sample=",sample,", type=",type,", skip evaluation. Suggesting: please do the prediction first.",sep=""))
          }
          next
        }
        
        if(file.info(trueFile)$size==0){
          print(paste("Empty tool file for sample=",sample,", type=",type,", skip evaluation.",sep=""))
          next
        }
        
        ## prepare data
        toolSV=read.table(toolFile,sep="\t",as.is=T,header=F)
        trueSV=read.table(trueFile,sep="\t",as.is=T,header=F)
        
        toolLen=0  # total tool SV length in bp
        trueLen=0  # total true SV length in bp
        overlapLen=0  # total overlapped length in bp
        
        for(chr in para$Chr){
          # subfile for each chr
          
          # set tool score for first col & true score for second col 
          # in case the scoring algorithm is changing. Everything larger than 0 can indicate 
          toolPos=as.matrix(toolSV[toolSV[,1]==chr,2:3])
          toolScore=matrix(0,nrow=nrow(toolPos),ncol=2)
          toolScore[,1]=1  # tool score for first col
          
          truePos=as.matrix(trueSV[trueSV[,1]==chr,2:3])
          trueScore=matrix(0,nrow=nrow(truePos),ncol=2)
          trueScore[,2]=1  # true score for second col
          
          # break into pieces
          SVpos=rbind(toolPos,truePos)
          SVscore=rbind(toolScore,trueScore)
          if(nrow(SVpos)==0){
            next
          }
          pieceRes=SV2piece(SVpos=SVpos, SVscore=SVscore)
          
          # count base pair
          pieceLen=pieceRes$piecePos[,2]-pieceRes$piecePos[,1]+1
          toolInd=which(pieceRes$pieceScore[,1]>0)
          trueInd=which(pieceRes$pieceScore[,2]>0)
          overlapInd=which(pieceRes$pieceScore[,1]>0 & pieceRes$pieceScore[,2]>0)
          
          toolLen=toolLen + sum(pieceLen[toolInd])
          trueLen=trueLen + sum(pieceLen[trueInd])
          overlapLen=overlapLen + sum(pieceLen[overlapInd])
        }
        
        ## evaluation
        evaMat[sample,"Predict"]=toolLen
        evaMat[sample,"True"]=trueLen
        evaMat[sample,"Overlap"]=overlapLen
        evaMat[sample,"Precision"]=overlapLen/toolLen
        evaMat[sample,"Recall"]=overlapLen/trueLen
        evaMat[is.na(evaMat)]=0 # in case 0 denominator
      } # end for(sample in Sample)
      
      # write evaMat to file
      evaFile=paste(para$outDir,"/",para$evaDir,"/",tool,"_",type,".evaBP.txt",sep="")
      write.table(data.frame(Sample=Sample,evaMat),file=evaFile,col.names=T,row.names=F,sep="\t",quote=F)
      
    } # end for(type in para$Type)
  }  # end for(tool in para$evaCaller)
  
}



### callerSVevaluationVCF
# to evaluate the caller results compared with true set
# SV level evaluation, precision and recall
# take in VCF file directly, not formated files
callerSVevaluationVCF <- function(para,toolVCF,inputToolDir,toolName,trueVCF,inputTrueDir,trueName,Sample){
  
  # toolVCF - bool whether tool is in VCF format 
  # inputToolDir - directory for tool VCF files
  # toolName - name of the tool
  # trueVCF - bool whether true is in VCF format
  # inputTrueDir - directory for true VCF files
  # trueName - name of the truth
  # Sample - total sample to be compared
  
  print(paste("[",Sys.time(),"] SV level evaluation on VCF files",sep=""))
  
  system(paste("mkdir -p ",para$outDir,"/",para$evaDir,sep=""))
  
  if(toolVCF==TRUE){
    system(paste("mkdir -p ",para$tmpDir,"/",toolName,sep=""))
  }
  
  if(trueVCF==TRUE){
    system(paste("mkdir -p ",para$tmpDir,"/",trueName,sep=""))
  }
  
  ## evaluation initialization
  colName=c("Predict","True","TrueHitPredict","Precision","PredictHitTrue","Recall")
  evaMat0=matrix(0,nrow=length(Sample),ncol=length(colName))
  rownames(evaMat0)=Sample
  colnames(evaMat0)=colName
  
  ## tool caller evaluation
  tool=toolName
  print(paste("Evaluating ",tool,sep=""))
  
  for(type in para$Type){
    
    evaMat=evaMat0
    
    for(sample in Sample){
      
      ## true file
      prepBed=paste(para$tmpDir,"/",trueName,"/",sample,"_",type,".bed",sep="")
      if(trueVCF==TRUE){
        # true file in VCF format
        if(!file.exists(prepBed)){
          # prepare SV learning files
          tmp=system(paste("ls ",inputTrueDir,"/",sample,"*vcf",sep=""),intern=T)
          if(length(tmp)==1){
            vcfFile=tmp
          }else{
            stop(paste("Multiple true ",trueName," OR no VCF files for sample=",sample,sep=""))
          }
          
          # read in vcf to prepare
          res=read.table(vcfFile,sep="\t",as.is=T,header=F,comment.char="#")
          resSub=res[which(res[,5]==paste("<",type,">",sep="")),]
          if(nrow(resSub)>0){
            tmp=strsplit(resSub[,8],split=";")
            SVENDtmp=sapply(tmp,function(x) return(x[startsWith(x,"END=")]))
            SVend=as.integer(gsub("END=","",as.character(SVENDtmp)))
            df=data.frame(chr=resSub[,1],startPos=as.integer(resSub[,2]),endPos=SVend,type=type,stringsAsFactors=F)
            
            ## check endPos>startPos
            switchInd=which(df$endPos<df$startPos)
            tmp=df$startPos[switchInd]
            df$startPos[switchInd]=df$endPos[switchInd]
            df$endPos[switchInd]=tmp
            
            # filter the data by chr
            keepInd=which(as.character(df$chr)%in%para$Chr)
            df=df[keepInd,]
            
            # filter the data with start=end
            keepInd=which(df$startPos < df$endPos)
            df=df[keepInd,]
            
            # filter the data when end exceeding chromosome size
            tmpChrSize=as.numeric(para$refLen[as.character(df$chr)])
            correctInd=which(df$startPos>tmpChrSize)
            if(length(correctInd)>0){
              warning(paste(length(correctInd), " of the ",tool," SVs in sample ",sample," have start position larger than chromosome size, remove these SVs",sep=""))
              df=df[-correctInd,]
              tmpChrSize=tmpChrSize[-correctInd]
            }
            correctInd=which(df$endPos>tmpChrSize)
            if(length(correctInd)>0){
              warning(paste(length(correctInd)," of the ",tool," SVs in sample ", sample," have end position larger than chromosome size, modify them",sep=""))
              df[correctInd,"endPos"]=tmpChrSize[correctInd]
            }
            
            write.table(df,file=prepBed,sep="\t",col.names=F,row.names=F,quote=F)
          }
        }  # end if(!file.exists(prepBed))
      } # end if(trueVCF==TRUE)
      
      trueFile=prepBed
      trueDir=paste(para$tmpDir,"/",trueName,sep="")
      
      if(file.exists(trueFile)==F){
        print(paste("No true ",trueName," file for sample=",sample,", type=",type,", skip evaluation.",sep=""))
        next
      }
      
      if(file.info(trueFile)$size==0){
        print(paste("Empty true ", trueName," file for sample=",sample,", type=",type,", skip evaluation.",sep=""))
        next
      }
      
      
      ## tool file
      prepBed=paste(para$tmpDir,"/",toolName,"/",sample,"_",type,"_ori.bed",sep="")
      
      if(toolVCF==TRUE){
        # tool file in VCF format
        
        if(!file.exists(prepBed)){
          # prepare SV learning files
          tmp=system(paste("ls ",inputToolDir,"/",sample,"*vcf",sep=""),intern=T)
          if(length(tmp)==1){
            vcfFile=tmp
          }else{
            stop(paste("Multiple tool ", toolName, " OR no VCF files for sample=",sample,sep=""))
          }
          
          # read in vcf to prepare
          res=read.table(vcfFile,sep="\t",as.is=T,header=F,comment.char="#")
          
          if(toolName=="manta"){
            svtype=sapply(strsplit(res[,3],split=":"), function(x) return(gsub("Manta","",x[1])))
            resSub=res[which(svtype==type),]
          }else{
            resSub=res[which(res[,5]==paste("<",type,">",sep="")),]
          }
          
          
          
          if(nrow(resSub)>0){
            tmp=strsplit(resSub[,8],split=";")
            SVENDtmp=sapply(tmp,function(x) return(x[startsWith(x,"END=")]))
            SVend=as.integer(gsub("END=","",as.character(SVENDtmp)))
            df=data.frame(chr=resSub[,1],startPos=as.integer(resSub[,2]),endPos=SVend,type=type,stringsAsFactors=F)
            
            ## check endPos>startPos
            switchInd=which(df$endPos<df$startPos)
            tmp=df$startPos[switchInd]
            df$startPos[switchInd]=df$endPos[switchInd]
            df$endPos[switchInd]=tmp
            
            # filter the data by chr
            keepInd=which(as.character(df$chr)%in%para$Chr)
            df=df[keepInd,]
            
            # filter the data with start=end
            keepInd=which(df$startPos < df$endPos)
            df=df[keepInd,]
            
            # filter the data when end exceeding chromosome size
            tmpChrSize=as.numeric(para$refLen[as.character(df$chr)])
            correctInd=which(df$startPos>tmpChrSize)
            if(length(correctInd)>0){
              warning(paste(length(correctInd), " of the ",tool," SVs in sample ",sample," have start position larger than chromosome size, remove these SVs",sep=""))
              df=df[-correctInd,]
              tmpChrSize=tmpChrSize[-correctInd]
            }
            correctInd=which(df$endPos>tmpChrSize)
            if(length(correctInd)>0){
              warning(paste(length(correctInd)," of the ",tool," SVs in sample ", sample," have end position larger than chromosome size, modify them",sep=""))
              df[correctInd,"endPos"]=tmpChrSize[correctInd]
            }
            
            write.table(df,file=prepBed,sep="\t",col.names=F,row.names=F,quote=F)
          }
        }  # end if(!file.exists(prepBed))
        
      } # end if(toolVCF==TRUE)
      
      toolFile=prepBed
      toolDir=paste(para$tmpDir,"/",toolName,sep="")
      
      if(file.exists(toolFile)==F){
        print(paste("No tool ",toolName," file for sample=",sample,", type=",type,", skip evaluation.",sep=""))
        next
      }
      
      if(file.info(trueFile)$size==0){
        print(paste("Empty tool ",toolName," file for sample=",sample,", type=",type,", skip evaluation.",sep=""))
        next
      }
      
      ##### mark to compare
      ## mark file 1: true hit predict
      markFile1=paste(toolDir,"/",sample,"_",type,"_mark1.txt",sep="")
      
      # intersect with truth
      s=paste(para$bedtools," intersect -a ", toolFile, " -b ",trueFile," -f ",para$evaOverlapRate," -r -wa -c > ",markFile1,sep="")
      system(s)
      
      # read in mark file 1 to evaluate
      res=read.table(markFile1,header=F,sep="\t",as.is=T)
      evaMat[sample,"Predict"]=nrow(res)
      evaMat[sample,"TrueHitPredict"]=sum(res[,ncol(res)]>0)
      evaMat[sample,"Precision"]=evaMat[sample,"TrueHitPredict"]/evaMat[sample,"Predict"]
      
      ## mark file 2: predict hit true
      markFile2=paste(toolDir,"/",sample,"_",type,"_mark2.txt",sep="")
      
      # intersect with truth
      s=paste(para$bedtools," intersect -a ", trueFile, " -b ",toolFile," -f ",para$evaOverlapRate," -r -wa -c > ",markFile2,sep="")
      system(s)
      
      # read in the mark file 2 to evaluate
      res=read.table(markFile2,header=F,sep="\t",as.is=T)
      evaMat[sample,"True"]=nrow(res)
      evaMat[sample,"PredictHitTrue"]=sum(res[,ncol(res)]>0)
      evaMat[sample,"Recall"]=evaMat[sample,"PredictHitTrue"]/evaMat[sample,"True"]
      evaMat[is.na(evaMat)]=0 # in case 0 denominator
    } # end for(sample in Sample)
    
    # write evaMat to file
    evaFile=paste(para$outDir,"/",para$evaDir,"/",toolName,"_",trueName,"_",type,".evaSV.txt",sep="")
    write.table(data.frame(Sample=Sample,evaMat),file=evaFile,col.names=T,row.names=F,sep="\t",quote=F)
    
  } # end for(type in para$Type)
  
}



