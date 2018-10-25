
#########################################################
## This pipeline analyse GeneChip data from affymetrix ##
#########################################################

##### Packages installation #####

install.packages("hexbin")
source("https://bioconductor.org/biocLite.R")
biocLite("affy")
biocLite("limma")

##### Packages loading #####

library(affy)
library(limma)
library(hexbin)

##### Reading CEL files #####

#maybe try to make a script usable in command line ans thus make pathToData as a parameter
pathToData <- "~/Desktop/microArrayExperiment/GSE99802/"
setwd(pathToData) ## We are obligated to go to the directory since ReadAffy seems not able to work outside of it

celFileNames <- list.celfiles()
data <- ReadAffy(filenames = celFileNames,compress = TRUE)

##### Quality control #####

data.expr <- exprs(data)
plot(hexplom(log2(data.expr))) #skeep this step during test since it last a very long time (using almost all the memory of the pc)

##### Normalisation #####

data.rma <- rma(data)

##### LIMMA analysis #####

