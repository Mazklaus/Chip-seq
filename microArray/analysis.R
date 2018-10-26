
#########################################################
## This pipeline analyse GeneChip data from affymetrix ##
#########################################################

##### Packages installation #####

install.packages("hexbin")
install.packages("statmod")
source("https://bioconductor.org/biocLite.R")
biocLite("affy")
biocLite("limma")

##### Packages loading #####

library(affy)
library(limma)
library(hexbin)
library(statmod)

##### Reading CEL files #####

#maybe try to make a script usable in command line ans thus make pathToData as a parameter
pathToData <- "~/Desktop/microArrayExperiment/GSE99802/"
setwd(pathToData) ## We are obligated to go to the directory since ReadAffy seems not able to work outside of it

celFileNames <- list.celfiles()
data <- ReadAffy(filenames = celFileNames,compress = TRUE)
sampleSheet <- read.csv("samplesheet.csv.txt",stringsAsFactors = FALSE)
sampleSheet$treatment <- as.factor(sampleSheet$treatment)
sampleSheet$status <- as.factor(sampleSheet$status)
sampleSheet$patientID <- as.factor(sampleSheet$patientID)

name <- strsplit(sampleNames(data),"_")
name <- as.data.frame(t(data.frame(name))[,1])
sampleNames(data) <- name[,1]


##### Quality control #####

data.expr <- exprs(data)
plot(hexplom(log2(data.expr))) #skeep this step during test since it last a very long time (using almost all the memory of the pc)

##### Normalisation #####

data.rma <- rma(data)

##### LIMMA analysis #####

##Computing
# Create a design matrix for treatments effects
design.trt <- model.matrix(~0+sampleSheet$treatment)

# Compute blocks for each genes (block = patients)
corfit <- duplicateCorrelation(data.rma,design.trt,block = sampleSheet$patientID)

# correlation histogram
hist(tanh(corfit$atanh.correlations))

#Compute the pooled sample variance and the sample mean by genes
fitTrtMean <- lmFit(data.rma, design.trt, block = sampleSheet$patientID, correlation = corfit$consensus.correlation)

##Create the coefficient matrix for the contrast
colnames(design.trt) <- c("trtDrug","trtPlacebo")
contrast.matrix = makeContrasts(trtDrug - trtPlacebo,levels = design.trt)

## Estimates the constrast
fit.contrast <- contrasts.fit(fitTrtMean, contrast.matrix)
