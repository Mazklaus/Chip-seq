
#########################################################
## This pipeline analyse GeneChip data from affymetrix ##
#########################################################

##### Packages installation #####

source("https://bioconductor.org/biocLite.R")
biocLite("affy")
biocLite("limma")

##### Packages loading #####

library(affy)
library(limma)

##### Reading CEL files #####

#maybe try to make a script usable in command line ans thus make pathToData as a parameter
pathToData <- "~/../../mnt/c/Users/team5.DESKTOP-GQ2GHLO/Desktop/travail_abel/projet/dataMining_IL22BP/RNA-seq_data/GSE99802"

celFileNames <- list.celfiles(pathToData)
data <- ReadAffy(filenames = celFileNames, celfile.path = pathToData,compress = TRUE)
