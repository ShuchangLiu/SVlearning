###  main script for SVlearning

## system setting
# module load gcc/4.9.2
# module load bedtools/2.27.0
# module load R/3.4.1
# cd /projects/liusi/DP/Gitpackage/SVlearning/scripts/

  
rm(list=ls())

## load source function
source(paste(para$scriptDir,"/Rsource.r",sep=""))

## read in config file to set parameters
#args=commandArgs(TRUE)
#configFileName=args[1]
configFileName="../example/code/config"  ## change here, for testing 

## check and format the input paramter
para=checkFormatPara(configFileName)

## install R package 
installRpackage(para)

## prepare data from vcf to python input file
prepareData(para)

## piece-SV data format
formatPieceSV(para)

## train the model or not
if(para$knownModel==FALSE){
  # no known model, need to train the model
  set.seed(06032)
  trainModel(para,Sample=para$trainSample)
}

## apply the model to testing data
modelPredict(para,Sample=para$testSample)

## caller evaluation
if(para$reportTestEvaluation==TRUE){
  callerEvaluation(para,Sample=para$testSample)
}




