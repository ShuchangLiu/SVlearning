###  main script for SVlearning

module load gcc/4.9.2
module load bedtools/2.27.0
module load R/3.4.1
cd /projects/liusi/DP/Gitpackage/SVlearning/scripts/

  
rm(list=ls())

## read in config file to set parameters
#args=commandArgs(TRUE)
#configFileName=args[1]
configFileName="../example/code/config"  ## change this part later-on
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


print(para)

## load source function
source(paste(para$scriptDir,"/Rscource.r",sep=""))
## TBD load in python function

## prepare data from vcf to python input file
prepareData(para)

## SV -> piece-SV by python code
# call python code

## piece-SV data format
formatPieceSV(para)

## train the model or not
if(para$knownModel==FALSE){
  # no known model, need to train the model
  trainModel(para)
}


## apply the model to testing data
modelPreidct(para)





