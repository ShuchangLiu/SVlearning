# SVlearning
To integrate the SV call sets by Machine Learning Method

SVlearning is a Linux/Unix system based R script pipeline to integrate multiple SV call sets into one final call set by machine learning methods. Samples are divided into two groups: well-studied training data with known true SV callings, and new unknown samples to be explored. Both training and new samples will first be pre-processed separately by filtering, score extraction, categorization and segmentation. Then the SVlearning models are trained from pilot training samples, and finally applied into new unknown samples for structural variation prediction. If a known model is provided, the training steps will be skipped and the model can be used for new sample prediction directly. Six machine learning algorithms are available for the current SVlearning pipeline: (1) support vector machines with radial kernel, (2) support vector machines with polynomial kernel, (3) random forest, (4) neural network, (5) linear discriminant analysis and (6) Adaptive Boosting. Users can choose any one or more algorithm combinations from them. 

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

Please refer to each SV caller website for details. The software can take in any caller combination. Take the combination of Breakdancer, CNVnator and delly as an example.

* [BreakDancer](https://github.com/genome/breakdancer) (v1.4.5)
* [CNVnator](https://github.com/abyzovlab/CNVnator) (v0.3.3)
* [Delly](https://github.com/dellytools/delly) (v2)

Alternatively, please refer to our tool [SVE (Structural Variation Engine)](https://github.com/TheJacksonLaboratory/SVE) for easier running on all these three callers (and other four callers).

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
All the parameter settings are included in config file. This file is "=" separated with pound sign as comment. Please refer to the config file in the example folder and modify it accordingly. Space will not be counted.

#### Model file
If the model file is not provided, users need to set config file parameter knownModel=FALSE and provide the training data. The pipeline will train and generate a model file. Otherwise, if the model file is provided, users can set knownModel=TRUE and provide the model file path in the configuration file.

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
The output can be in both BED and VCF format (set outFormat to be BED,VCF) in the output directory.


#### BED format output
The files will be located at outDir/bedDir as illustrated in config file. The file is tab separabted and the format is as follows

Header | Illustration
--- | ---
chrom | chromosome
chromStart | start position 
chromEnd | end position
name | SVlearning SV ID
REF | reference for chromStart position
ALT | SV type, either 'DEL' or 'DUP' for deletion and duplication
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

The INFO part is colon separated within each comoposing piece, and comma separated between pieces if a SV is composed by multiple pieces.

#### VCF format output
The files will be located at outDir/vcfDir as illustrated in config file. The file is illustrated by the header of the VCF file. The FORMAT column is the same as BED INFO column. 

#### Evaluation 
In the config file, if users set reportTestEvaluation = TRUE and provide the new data (testing data) true VCF call set, then SVlearning will evaluate the performance in the 'evaDir' directory. If reportTestEvaluation = FALSE, then no testing data true set required.

## Simple example
Please check SVlearning/example folder to test the pipeline. Users can go to SVlearning/example/code for the command guide and modify the config file accordingly. The testing data are VCF files called from BAM files in 1000GP project. Users can try their own data following the framework of the simple example.

## Algorithm details

* **Model file:** If the model is unknown (knownModel=FALSE), the pipeline will train the model by training data; alternatively, users can set knownModel=TRUE and provide model file (modelDir/modelFile) to skip the training step. In this pipeline, we prepared a ready-made model based on 27 1000G WGS samples and their phase 3 true sets. If the users do not have training data available, this model file is suggested as default one.
* **Working pipeline:** SVlearning will first prepare the data and cut the the SVs into pieces based on all the margins of SVs called by all the three callers. And then each piece is regarded as one unit to apply machine learning method to predict its true of false. Then neighboring predicted-true pieces will merge back into final SVs.
* **Machine learning (ML) models:** The ML models choosing from SVM with radial kernel, SVM with polynomial kernel, Random Forest, Neural Network, Linear discriminant analysis and Adaboost.
* **Voting models:** The pipeline takes the majority vote (cutoff determined by 'vote' in config file) among the above ML methods.

## Contact

Please contact Silvia Liu (silvia dot shuchang dot liu at gmail dot com) for bug reporting or question discussion.


