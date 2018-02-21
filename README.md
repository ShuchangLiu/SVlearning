# SVlearning
To merge the SV callers by Machine Learning Method

SVlearning is a Linux/Unix system based R script pipeline to merge multiple SV callers into one final call set by machine learning methods. The program takes in BreakDancer, CNVnator and Delly VCF files, and can apply and combine machine learning models choosing from SVM with radial kernel, SVM with polynomial kernel, Random Forest, Neural Network, Linear discriminant analysis and Adaboost.

Please read the whole README file before applying the pipeline.

## Requirements

* R 3.4.1 or greater
* R packages: nnet, e1071, MASS, randomForest, fastAdaboost (the pipeline can automatically install these packages at the first time)
* gcc 4.8 or greater
* [bedtools](http://bedtools.readthedocs.io/en/latest/index.html) 2.27.0 or greater

## Installation from GitHub

```bash
git clone https://github.com/ShuchangLiu/SVlearning.git
```

## Usage
### Input data preparation
#### Caller VCF files

Please refer to each SV caller website for details
* [BreakDancer](https://github.com/genome/breakdancer) (v1.4.5)
* [CNVnator](https://github.com/abyzovlab/CNVnator) (v0.3.3)
* [Delly](https://github.com/dellytools/delly) (v2)

Alternatively, please refer to our tool [SVE (Structural Variation Engine)](https://github.com/TheJacksonLaboratory/SVE) for easier running on all these three callers.

#### VCF file organization
For example, the VCF file is organized as follows,

* vcfFiles/sample1/sample1_true.vcf (required for training data, optional for testing data if reportTestEvaluation==FALSE in config file)
* vcfFiles/sample1/sample1_breakdancer.vcf
* vcfFiles/sample1/sample1_CNVnator.vcf
* vcfFiles/sample1/sample1_delly.vcf
* vcfFiles/sample2/sample2_true.vcf (required for training data, optional for testing data if reportTestEvaluation==FALSE in config file)
* vcfFiles/sample2/sample2_breakdancer.vcf
* vcfFiles/sample2/sample2_CNVnator.vcf
* vcfFiles/sample2/sample2_delly.vcf
* ... all the other samples

#### Config file 
All the parameter settings are included in config file. This file is tab separated with pound sign as comment. Please refer to the config file in the example folder and modify it accordingly. 

### Run the pipeline

Users can adjust the pipeline SVlearning.r and SVlearning_CV.r accordingly. 

For the general model training and sample prediction pipeline, please refer to the following command

```bash
Rscript path/to/SVlearning.r path/to/config path/to/Rsource.r
```

For the cross validation (CV), please refer to the following command

```bash
Rscript path/to/SVlearning_CV.r path/to/config path/to/Rsource.r
```

### Output file
The output can be in both BED and VCF format (set outFormat to be BED,VCF).

#### BED format output
The files will be located at outDir/bedDir as illustrated in config file. The file is tab separabted and the format is as follows

Header | Illustration
--- | ---
chrom | chromosome
chromStart | start position 
chromEnd | end position
name | SVlearning SV ID
REF | reference for chromStart position
ALT | SV type, either <DEL> or <DUP> for deletion and duplication
INFO-PSTART | start position of the composing piece
INFO-PEND | end position of the composing piece
INFO-breakdancer | breakdancer QUAL score for the composing piece
INFO-CNVnator | CNVnator NatorP2 score for the composing piece
INFO-delly | delly genotype quality (GT) score for the composing piece
INFO-NN | Neural network prediction (0: predict as false; 1: predict as true)
INFO-SVMpolynomial | SVM with polynomial kernel prediction (0: predict as false; 1: predict as true)
INFO-SVMradial | SVM with radial kernel prediction (0: predict as false; 1: predict as true)
INFO-LDA | Linear discriminant analysis prediction (0: predict as false; 1: predict as true)
INFO-RF | Random Forest prediction (0: predict as false; 1: predict as true)
INFO-adaboost | adaboost prediction (0: predict as false; 1: predict as true)

The INFO part is colon separated within each comoposing piece and comma separated between pieces if a SV is composed by multiple pieces.

#### VCF format output
The files will be located at outDir/vcfDir as illustrated in config file. The file is illustrated by the header of the VCF file. The FORMAT column is the same as BED INFO column. 

## Algorithm details

* Model file. If the model is unknown (knownModel==FALSE), the pipeline will train the model by training data; alternatively, users can set knownModel==TRUE and provide model file (modelDir/modelFile) to skip the training step. In this pipeline, we prepared a ready-made model based on 27 1000G samples based on their phase 4 true set. If the users don't have training data available, this model file is suggested as default one.
* Working pipeline: SVlearning will first prepare the data and cut the the SVs into pieces based on all the margins of SVs called by all the three callers. And then each piece is regarded as one unit to apply machine learning method to predict its true of false. Then neighboring predicted-true pieces will merge back into final SVs.
* Machine learning (ML) models. The ML models choosing from SVM with radial kernel, SVM with polynomial kernel, Random Forest, Neural Network, Linear discriminant analysis and Adaboost.
* Voting models. The pipeline takes the majority vote (cutoff determined by 'vote' in config file) among the above ML methods.

## Contact

Please contact Silvia Liu (silvia dot shuchang dot liu at gmail dot com) for bug reporting or question discussion.


