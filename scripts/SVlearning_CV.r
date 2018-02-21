### cross validation for SVlearning
# 26 sample as training + 1 testing for each CV

## system setting
# module load gcc/4.9.2
# module load bedtools/2.27.0
# module load R/3.4.1
# cd /projects/liusi/DP/Gitpackage/SVlearning/scripts/


rm(list=ls())

args=commandArgs(TRUE)

## load source function
sourceFile=args[2]
source(sourceFile)
#source(paste(para$scriptDir,"/Rsource.r",sep=""))


## read in config file to set parameters
configFileName=args[1]
#configFileName="../example/code/config"  ## change here, for testing 

## check and format the input paramter
para=checkFormatPara(configFileName)

## install R package 
installRpackage(para)

## prepare data from vcf to python input file
prepareData(para)

## piece-SV data format
formatPieceSV(para)

### cross validation
set.seed(06032)
totalSample=para$trainSample
for(i in 1:length(totalSample)){
  print(paste("[",Sys.time(),"] CV=1",i,", testing sample=",totalSample[i],sep=""))
  
  ## train the model
  para$modelFile=paste("model",i,".rdata",sep="")
  trainModel(para,Sample=totalSample[-i])
  
  ## apply the model to testing data
  modelPredict(para,Sample=totalSample[i])
  
}

## caller evaluation
if(para$reportTestEvaluation==TRUE){
  callerSVevaluation(para,Sample=para$testSample)
  callerBPevaluation(para,Sample=para$testSample)
}

## remove tmp file
if(para$rmTmpFile==TRUE){
  system(paste("rm -rf ",para$tmpDir,sep=""))
}



