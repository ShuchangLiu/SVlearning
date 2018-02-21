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
  para$MLmethod=c("NN","SVMpolynomial","SVMradial","LDA","RF","adaboost")
  
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
    trueVCF=paste(sampleDir,"/",sample,"/",trueVCFtmp,sep="")
    trueBEDallPre=paste(para$tmpDir,"/",para$trueDir,"/",sample,sep="")
    
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
    
    df=data.frame(chr=invcf$CHROM,startPos=startPos,endPos=endPos,name=SVtype,stringsAsFactors=FALSE)
    
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
        warning(paste(length(correctInd), " of the ",tool," SVs in sample ",sample," have start position larger than chromosome size, remove these SVs",sep=""))
        df=df[-correctInd,]
        tmpChrSize=tmpChrSize[-correctInd]
      }
      correctInd=which(df$endPos>tmpChrSize)
      if(length(correctInd)>0){
        warning(paste(length(correctInd)," of the ",tool," SVs in sample ", sample," have end position larger than chromosome size, modify them",sep=""))
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
        warning(paste(length(correctInd), " of the ",tool," SVs in sample ",sample," have start position larger than chromosome size, remove these SVs",sep=""))
        df=df[-correctInd,]
        tmpChrSize=tmpChrSize[-correctInd]
      }
      correctInd=which(df$endPos>tmpChrSize)
      if(length(correctInd)>0){
        warning(paste(length(correctInd)," of the ",tool," SVs in sample ",sample," have end position larger than chromosome size, modify them",sep=""))
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
        #trueBEDtype=paste(trueDir,"/",sample,"_",type,".bed",sep="")
        #toolTrueMarktype=paste(toolDir,"/",sample,"_",type,"_mark.txt",sep="")
        
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
        warning(paste(length(correctInd), " of the ",tool," SVs in sample ",sample," have start position larger than chromosome size, remove these SVs",sep=""))
        df=df[-correctInd,]
        tmpChrSize=tmpChrSize[-correctInd]
      }
      correctInd=which(df$endPos>tmpChrSize)
      if(length(correctInd)>0){
        warning(paste(length(correctInd)," of the ",tool," SVs in sample ",sample," have end position larger than chromosome size, modify them",sep=""))
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


#### prepareData
# to prepare the data for 
prepareData <- function(para){
  
  system(paste("mkdir -p ",para$tmpDir,sep=""))
  system(paste("mkdir -p ",para$tmpDir,"/",para$trueDir,sep=""))
  system(paste("mkdir -p ",para$tmpDir,"/",para$breakdancerDir,sep=""))
  system(paste("mkdir -p ",para$tmpDir,"/",para$CNVnatorDir,sep=""))
  system(paste("mkdir -p ",para$tmpDir,"/",para$dellyDir,sep=""))
  
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
  # in case testing data have truth available, for future testing results evaluation
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
  
  pieceDir=paste(para$tmpDir,"/",para$pieceDir,sep="")
  
  #[[type]][[sample]][[as.character(as.integer(binlen))]][[as.character(chr)]]
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
      
      trueBEDtype=paste(para$tmpDir,"/",para$trueDir,"/",sample,"_",type,".bed",sep="")
      
      tmp=list.files(path=pieceDir,pattern=paste(type,"_",sample,sep=""))
      if(length(tmp)!=0){
        print(paste("type=",type," sample=",sample," has been pieced, skip",sep=""))
        next
      }else{
        print(paste("type=",type," sample=",sample,sep=""))
        ## read in the data
        # breakdancer
        toolBEDtype=paste(para$tmpDir,"/",para$breakdancerDir,"/",sample,"_",type,".bed",sep="")
        res1=read.table(toolBEDtype,sep="\t",header=F,as.is=T)
        
        # CNVnator
        toolBEDtype=paste(para$tmpDir,"/",para$CNVnatorDir,"/",sample,"_",type,".bed",sep="")
        res2=read.table(toolBEDtype,sep="\t",header=F,as.is=T)
        
        # modify CNVnator to avaoid some extreme value
        ind1=res2[,5]<0
        ind2=res2[,5]==0
        res2[ind1,5]=0
        res2[ind2,5]=-1
        
        # delly
        toolBEDtype=paste(para$tmpDir,"/",para$dellyDir,"/",sample,"_",type,".bed",sep="")
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
            
            pieceFileXY=paste(pieceDir,"/",type,"_",sample,"_",as.character(as.integer(binlen)),"_",chr,"_xy.txt",sep="")
            selectInd=which(SVlenAll>=binLenAll[binlenInd] & SVlenAll<binLenAll[binlenInd+1] & SVposAll[,1]==chr)
            if(length(selectInd)>0){
              pieceOut=SV2piece(SVpos=as.matrix(SVposAll[selectInd,2:3]),SVscore=matrix(SVscoreAll[selectInd,],ncol=ncol(SVscoreAll)))
              out=data.frame(chr,pieceOut$piecePos,type,pieceOut$pieceScore,stringsAsFactors=FALSE)
              
              # intersect with truth
              if(file.exists(trueBEDtype)){
                pieceFile=paste(pieceDir,"/",type,"_",sample,"_",as.character(as.integer(binlen)),"_",chr,".bed",sep="")
                pieceFileMark=paste(pieceDir,"/",type,"_",sample,"_",as.character(as.integer(binlen)),"_",chr,"_mark.txt",sep="")
                
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
  
  system(paste("mkdir -p ",para$tmpDir,sep=""))
  system(paste("mkdir -p ",para$tmpDir,"/",para$pieceDir,sep=""))
  
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



#### trainModel
# a function to train the SVlearning model
trainModel <- function(para, Sample){
  
  print(paste("[",Sys.time(),"] Train the model",sep=""))
  
  #Sample=para$trainSample
  
  model=list()
  #[[type]][[as.character(as.integer(binlen))]][[as.character(chr)]][[modelName]]
  
  modelNA=list()
  for(modelName in para$MLmethod){
    modelNA[[modelName]]=NA
  }
  
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
      
      for(chr in para$Chr){
        
        ### prepare training data
        x=c()
        y=c()
        
        for(sample in Sample){
          pieceFile=paste(para$tmpDir,"/",para$pieceDir,"/",type,"_",sample,"_",as.character(as.integer(binlen)),"_",chr,"_xy.txt",sep="")
          if(file.info(pieceFile)$size!=0){
            res=read.table(pieceFile,header=F,sep="\t",as.is=T)
            x=rbind(x,res[,5:7])
            y=c(y,res[,8])
          }
        }
        
        model[[type]][[as.character(as.integer(binlen))]][[as.character(chr)]]=modelNA
        
        if(length(y)==0){
          # no supporting training data
          print(paste("No supporting training data for type=",type,", binlen=",as.character(as.integer(binlen)),", chr=",chr,sep=""))
          next
        }else if(length(unique(y))==1){
          # only one y label
          print(paste("Only y=",unique(y)," for type=",type,", binlen=",as.character(as.integer(binlen)),", chr=",chr,sep=""))
          next
        }else if (sum(y==0)==1){
          print(paste("Not enough number of y=0 cases for type=", type,", binlen=",as.character(as.integer(binlen)),", chr=",chr,sep=""))
          next
        }else if(sum(y==1)==1){
          print(paste("Not enough number of y=1 cases for type=", type,", binlen=",as.character(as.integer(binlen)),", chr=",chr,sep=""))
          next
        }
        
        # make y=0 and y=1 equally
        y0Ind=which(y==0)
        y1Ind=which(y==1)
        if(length(y0Ind)>length(y1Ind)){
          # then make the y1Ind the same number as y0Ind
          y1Ind=sample(y1Ind,size=length(y0Ind),replace=T)
        }else if(length(y0Ind)<length(y1Ind)){
          # make the y0Ind the same number as y1Ind
          y0Ind=sample(y0Ind,size=length(y1Ind),replace=T)
        }
        
        x=x[c(y0Ind,y1Ind),]
        y=y[c(y0Ind,y1Ind)]
        colnames(x)=NULL
        x=as.matrix(x)
        
        ### train the model
        # to avoid constant features
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
          model[[type]][[as.character(as.integer(binlen))]][[as.character(chr)]][["SVMradial"]] = svm(x=xx, y=as.factor(y), kernel="radial", cost = 100, probability=TRUE)  # default cost
        }
        
        # model 2: SVM polynomial
        if("SVMpolynomial"%in%para$MLmethod){
          model[[type]][[as.character(as.integer(binlen))]][[as.character(chr)]][["SVMpolynomial"]] = svm(x=xx, y=as.factor(y), kernel="polynomial", cost = 100, probability=TRUE)  # default cost
        }
        
        # model 3: RF
        if("RF"%in%para$MLmethod){
          model[[type]][[as.character(as.integer(binlen))]][[as.character(chr)]][["RF"]]=randomForest(x=x,y=as.factor(y))
        }
        
        # model 4: NN
        if("NN"%in%para$MLmethod){
          yy=cbind(y,1-y)
          model[[type]][[as.character(as.integer(binlen))]][[as.character(chr)]][["NN"]] = nnet(x, yy,size = 2, rang = 0.1, decay = 5e-4, maxit = 200)  # can adjust parameters later
        }
        
        # model 5: LDA
        if("LDA"%in%para$MLmethod){
          model[[type]][[as.character(as.integer(binlen))]][[as.character(chr)]][["LDA"]] = lda(x=xx,grouping=as.factor(y))
        }
        
        # model 6: adaboost
        if("adaboost"%in%para$MLmethod){
          prepData=data.frame(x,Y=y)
          prepData$Y=factor(prepData$Y)
          model[[type]][[as.character(as.integer(binlen))]][[as.character(chr)]][["adaboost"]] = adaboost(Y ~., prepData,nIter=100)
        }
        
      } # end for(chr in para$Chr)
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
  #   pieceScore[,2] for right position
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
    SVpos=piecePos
    SVscore=paste(c(piecePos[1,],pieceScore[1,]),collapse=":")
  }else{
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
  callerNum=3  # currently three caller, change this later
  
  #### prediction within each bin
  for(sample in Sample){
    
    print(paste("[",Sys.time(),"] sample=",sample,sep=""))
    
    SVall=data.frame(stringsAsFactors=FALSE)
    
    for(type in para$Type){
      
      print(paste("sample=",sample,", type=",type,sep=""))
      
      SVtype=data.frame(stringsAsFactors=FALSE)
      
      if(type=="DEL"){
        binLen=para$DELbin
      }else if(type=="DUP"){
        binLen=para$DUPbin
      }else{
        stop(paste("Type ",type," not supporting yet!"))
      }
      binLenAll=c(binLen,Inf)
      
      for(binlen in binLen){
        
        for(chr in para$Chr){
          
          pieceFileXY=paste(para$tmpDir,"/",para$pieceDir,"/",type,"_",sample,"_",as.character(as.integer(binlen)),"_",as.character(chr),"_xy.txt",sep="")
          
          if(file.info(pieceFileXY)$size==0){ # no testing data
            next
          }
          
          res=read.table(pieceFileXY,header=F,sep="\t",as.is=T)
          
          x=as.matrix(res[,5:7])
          colnames(x)=NULL
          y=as.numeric(res[,8])
          
          ## a matrix to record predict results for each ML model
          y.predict=matrix(0,nrow=nrow(res),ncol=modelNum+1)
          colnames(y.predict)=c(para$MLmethod,"vote")
          
          # model 1: SVM radial
          modelName="SVMradial"
          if(modelName%in%para$MLmethod){
            currentModel=model[[type]][[as.character(as.integer(binlen))]][[as.character(chr)]][[modelName]]
            if(is.null(currentModel) | all(is.na(currentModel))){
              print(paste("No model for sample=",sample,", type=",type,", binlen=",binlen,", chr=",chr,", model=",modelName,", predict all as 0",sep=""))
            }else{
              y.predict[,modelName] = as.numeric(predict(currentModel, x, decision.values=T))-1
            }
          }
          
          # model 2: SVM polynomial
          modelName="SVMpolynomial"
          if(modelName%in%para$MLmethod){
            currentModel=model[[type]][[as.character(as.integer(binlen))]][[as.character(chr)]][[modelName]]
            if(is.null(currentModel) | all(is.na(currentModel))){
              print(paste("No model for sample=",sample,", type=",type,", binlen=",binlen,", chr=",chr,", model=",modelName,", predict all as 0",sep=""))
            }else{
              y.predict[,modelName] = as.numeric(predict(currentModel, x, decision.values=T))-1
            }
          }
          
          # model 3: RF
          modelName="RF"
          if(modelName%in%para$MLmethod){
            currentModel=model[[type]][[as.character(as.integer(binlen))]][[as.character(chr)]][[modelName]]
            if(is.null(currentModel) | all(is.na(currentModel))){
              print(paste("No model for sample=",sample,", type=",type,", binlen=",binlen,", chr=",chr,", model=",modelName,", predict all as 0",sep=""))
            }else{
              y.predict[,modelName] = as.numeric(predict(currentModel,x))-1
            }
          }
          
          # model 4: NN
          modelName="NN"
          if(modelName%in%para$MLmethod){
            currentModel=model[[type]][[as.character(as.integer(binlen))]][[as.character(chr)]][[modelName]]
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
            currentModel=model[[type]][[as.character(as.integer(binlen))]][[as.character(chr)]][[modelName]]
            if(is.null(currentModel) | all(is.na(currentModel))){
              print(paste("No model for sample=",sample,", type=",type,", binlen=",binlen,", chr=",chr,", model=",modelName,", predict all as 0",sep=""))
            }else{
              y.predict[,modelName] = as.numeric(predict(currentModel,x)$class)-1
            }
          }
          
          # model 6: adaboost
          modelName="adaboost"
          if(modelName%in%para$MLmethod){
            currentModel=model[[type]][[as.character(as.integer(binlen))]][[as.character(chr)]][[modelName]]
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
                         gapDis=para$gapDis)
            SVtype=rbind(SVtype,data.frame(chr=chr,startPos=tmp$SVpos[,1],endPos=tmp$SVpos[,2],type=type,info=tmp$SVscore,stringsAsFactors=FALSE))
          }
          
        } # end for(chr in para$Chr)
      } # end for(binlen in binLen)
      
      ## any code for sample & type bin, can be added here
      
      SVall=rbind(SVall, SVtype)
      
    } # end for (type in para$Type)
    
    ### prepare for output
    
    SVall$chr=factor(SVall$chr,levels=para$Chr)
    
    if(nrow(SVall)>0){
      ### combine all into one file, sort
      # by startPos
      ind=order(SVall$startPos)
      SVall=SVall[ind,]
      # by chr
      ind=order(SVall$chr)
      SVall=SVall[ind,]
      
      # get reference
      out=data.frame(SVall$chr,SVall$startPos,SVall$startPos+1,stringsAsFactors=FALSE)
      system(paste("mkdir -p ",para$tmpDir,"/",para$shortBedDir,sep=""))
      outFile=paste(para$tmpDir,"/",para$shortBedDir,"/",sample,".short.bed",sep="")
      write.table(out,file=outFile,row.names=F,col.names=F,sep="\t",quote=F)
      REF=as.character(system(paste(para$bedtools," getfasta -tab -fi ",para$reference," -bed ",outFile," | awk -F \"\\t\" '{print $2}' ",sep=""),intern=T))
      if(length(REF)!=nrow(SVall)){
        stop("Errors for reference: non-euqal number of reference and SV, SV position may exceed chromosome length.")
      }
      ## TBD: some called regions with multiple starting N or end N strings, do some final fitering
      
      # ID name
      ID=paste("SVlearing",1:nrow(SVall),sep="_")
      
      # ALT
      ALT=paste("<",SVall$type,">",sep="")
      
      # INFO
      INFO=paste("SVTYPE=",SVall$type,";SVLEN=",SVall$endPos-SVall$startPos+1,";END=",SVall$endPos,";IMPRECISE",sep="")
    }
    
    # BED format
    if("BED"%in%para$outFormat){
      system(paste("mkdir -p ",para$outDir,"/",para$bedDir,sep=""))
      outFile=paste(para$outDir,"/",para$bedDir,"/",sample,".bed",sep="")
      
      if(nrow(SVall)==0){
        cat(paste("#chrom","chromStart","chromEnd","name","REF","ALT",
                  paste("PSTART:PEND:breakdancer:CNVnator:delly:",paste(para$MLmethod,collapse=":"),sep=""),sep="\t"),file=outFile,sep="\n")
      }else{
        outSVall=data.frame(chrom=SVall$chr, chromStart=SVall$startPos, chromEnd=SVall$endPos, name=ID, 
                            REF=REF, ALT=ALT, INFO=SVall$info, stringsAsFactors=FALSE) 
        colnames(outSVall)[7]=paste("PSTART:PEND:breakdancer:CNVnator:delly:",paste(para$MLmethod,collapse=":"),sep="")
        
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
      cat("##FORMAT=<ID=breakdancer,Number=1,Type=Double,Description=\"Breakdancer quality score\">",file=outFile,sep="\n",append=TRUE)
      cat("##FORMAT=<ID=CNVnator,Number=1,Type=Double,Description=\"CNVnator natorP2 value\">",file=outFile,sep="\n",append=TRUE)
      cat("##FORMAT=<ID=delly,Number=1,Type=Double,Description=\"Delly genotype quality score\">",file=outFile,sep="\n",append=TRUE)
      
      for(methodName in para$MLmethod){
        cat(paste("##FORMAT=<ID=",methodName,",Number=1,Type=Integer,Description=\"",methodName," prediction\">",sep=""),file=outFile,sep="\n",append=TRUE)
      }
      
      for(refInd in 1:nrow(refLen)){
        cat(paste("##contig=<ID=",refLen[refInd,1],",len=",refLen[refInd,2],">",sep=""),file=outFile,sep="\n",append=TRUE)
      }
      
      if(nrow(SVall)==0){
        cat(paste("#CHROM","POS","ID","REF","ALT", "QUAL","FILTER","INFO","FORMAT",sample,sep="\t"),file=outFile,sep="\n")
      }else{
        
        outSVall=data.frame(CHROM=SVall$chr, POS=SVall$startPos, ID=ID, REF=REF, ALT=ALT, QUAL=".",FILTER="PASS", INFO=INFO,
                            FORMAT=paste("PSTART:PEND:breakdancer:CNVnator:delly:",paste(para$MLmethod,collapse=":"),sep=""),
                            formatInfo=SVall$info, stringsAsFactors=FALSE)
        colnames(outSVall)[10]=sample
        
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
        trueFile=paste(para$tmp,"/",para$trueDir,"/",sample,"_",type,".bed",sep="")
        if(file.exists(trueFile)==F){
          print(paste("No true file for sample=",sample,", type=",type,", skip evaluation. Suggestion: please try prepareData function first.",sep=""))
          next
        }
        
        if(file.info(trueFile)$size==0){
          print(paste("Empty True file for sample=",sample,", type=",type,", skip evaluation.",sep=""))
          next
        }
        
        ## tool file
        if(tool=="breakdancer"){
          toolDir=paste(para$tmp,"/",para$breakdancerDir,sep="")
        }else if(tool=="CNVnator"){
          toolDir=paste(para$tmp,"/",para$CNVnatorDir,sep="")
        }else if(tool=="delly"){
          toolDir=paste(para$tmp,"/",para$dellyDir,sep="")
        }else if(tool=="SVlearning"){
          prepBed=paste(para$tmpDir,"/",para$SVlearningDir,"/",sample,"_",type,".bed",sep="")
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
                write.table(outBed,file=prepBed,sep="\t",col.names=F,row.names=F,quote=F)
              }
            }else if(file.exists(vcfFile)){
              # read in vcf to prepare
              res=read.table(vcfFile,sep="\t",as.is=T,header=F,comment.char="#")
              resSub=res[which(res[,5]==paste("<",type,">",sep="")),]
              if(nrow(resSub)>0){
                tmp=strsplit(resSub[,8],split=";")
                SVENDtmp=sapply(tmp,function(x) return(x[startsWith(x,"END=")]))
                SVend=as.numeric(gsub("END=","",as.character(SVENDtmp)))
                outVCF=data.frame(resSub[,1:2],SVend,type=type,stringsAsFactors=F)
                write.table(outVCF,file=prepBed,sep="\t",col.names=F,row.names=F,quote=F)
              }
            }
          }
          toolDir=paste(para$tmpDir,"/",para$SVlearningDir,sep="")
          
        }else{
          stop(paste("Wrong tool name ",tool,sep=""))
        }
        
        toolFile=paste(toolDir,"/",sample,"_",type,".bed",sep="")
        
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
        s=paste(para$bedtools," intersect -a ", toolFile, " -b ",trueFile," -f 0.5 -r -c > ",markFile1,sep="")
        system(s)
        
        # read in mark file 1 to evaluate
        res=read.table(markFile1,header=F,sep="\t",as.is=T)
        evaMat[sample,"Predict"]=nrow(res)
        evaMat[sample,"TrueHitPredict"]=sum(res[,ncol(res)]>0)
        evaMat[sample,"Precision"]=evaMat[sample,"TrueHitPredict"]/evaMat[sample,"Predict"]
        
        ## mark file 2: predict hit true
        markFile2=paste(toolDir,"/",sample,"_",type,"_mark2.txt",sep="")
        
        # intersect with truth
        s=paste(para$bedtools," intersect -a ", trueFile, " -b ",toolFile," -f 0.5 -r -c > ",markFile2,sep="")
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
        trueFile=paste(para$tmp,"/",para$trueDir,"/",sample,"_",type,".bed",sep="")
        if(file.exists(trueFile)==F){
          print(paste("No true file for sample=",sample,", type=",type,", skip evaluation. Suggestion: please try prepareData function first.",sep=""))
          next
        }
        
        if(file.info(trueFile)$size==0){
          print(paste("Empty True file for sample=",sample,", type=",type,", skip evaluation.",sep=""))
          next
        }
        
        ## tool file
        if(tool=="breakdancer"){
          toolDir=paste(para$tmp,"/",para$breakdancerDir,sep="")
        }else if(tool=="CNVnator"){
          toolDir=paste(para$tmp,"/",para$CNVnatorDir,sep="")
        }else if(tool=="delly"){
          toolDir=paste(para$tmp,"/",para$dellyDir,sep="")
        }else if(tool=="SVlearning"){
          prepBed=paste(para$tmpDir,"/",para$SVlearningDir,"/",sample,"_",type,".bed",sep="")
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
                write.table(outBed,file=prepBed,sep="\t",col.names=F,row.names=F,quote=F)
              }
            }else if(file.exists(vcfFile)){
              # read in vcf to prepare
              res=read.table(vcfFile,sep="\t",as.is=T,header=F,comment.char="#")
              resSub=res[which(res[,5]==paste("<",type,">",sep="")),]
              if(nrow(resSub)>0){
                tmp=strsplit(resSub[,8],split=";")
                SVENDtmp=sapply(tmp,function(x) return(x[startsWith(x,"END=")]))
                SVend=as.numeric(gsub("END=","",as.character(SVENDtmp)))
                outVCF=data.frame(resSub[,1:2],SVend,type=type,stringsAsFactors=F)
                write.table(outVCF,file=prepBed,sep="\t",col.names=F,row.names=F,quote=F)
              }
            }
          }
          toolDir=paste(para$tmpDir,"/",para$SVlearningDir,sep="")
          
        }else{
          stop(paste("Wrong tool name ",tool,sep=""))
        }
        
        toolFile=paste(toolDir,"/",sample,"_",type,".bed",sep="")
        
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




