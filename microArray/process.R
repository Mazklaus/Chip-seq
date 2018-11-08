# This pipeline is made to process CELL file and analyse them #
# This pipeline is design to handle old affymetrix array #
# to see how to handle the new one refer to https://wiki.bits.vib.be/index.php/Analyze_your_own_microarray_data_in_R/Bioconductor #

##### Package installation ####

source("http://www.bioconductor.org/biocLite.R")

biocLite("affy")
biocLite("affyPLM")
biocLite("XML") # if don't work use before : sudo apt-get install libxml2-dev
biocLite("RCurl") # if don't work use before : sudo apt-get install libcurl4-openssl-dev
biocLite("simpleaffy")
biocLite("affydata")
biocLite("ArrayExpress")
biocLite("limma")
biocLite("Biobase")
biocLite("Biostrings")
biocLite("genefilter")

install.packages("gcrma")

#### Package load ####

library(affy)
library(affyPLM)
library(simpleaffy)
library(affydata)
library(ArrayExpress)
library(limma)
library(Biobase)
library(Biostrings)
library(genefilter)
library(gcrma)

#### Getting the CELL file data ####

# specify your working directory #
wd = "~/Desktop/microArrayExperiment/"
setwd(wd)

## Specify where are stored the CELL fies ##
celpath = "~/Desktop/microArrayExperiment/GSE99802/"

# import CELL files #
data <- ReadAffy(celfile.path = celpath)

## Specify the dirrectory to the sample sheet ##
sspath <- "~/Desktop/microArrayExperiment/GSE99802/samplesheet.csv.txt"
sampleSheet <- read.csv(sspath)

# Changing sample name in the data to better one #
ph <- data@phenoData

for (i in seq(1,nrow(sampleSheet))){
  sampleSheet$sampleName[i] <- paste(sampleSheet$patientID[i],sampleSheet$time[i],sampleSheet$status[i],sampleSheet$treatment[i],sep = "_")
}

ph@data[,1] <- sampleSheet$sampleName

##### Quality control ####

# microarray image #
dir.create("./QC")
dir.create("./QC/miccroarray_picture")

for (i in 1:nrow(sampleSheet)) {
  name = paste(wd,"QC/miccroarray_picture/","image",i,".jpg",sep="")
  jpeg(name)
  image(data[,i],main=ph@data$sample[i])
  dev.off()
}

# Chip pseudo image #
dir.create("./QC/pseudo_image")

Pset <- fitPLM(data,output.param=list(varcov="none"))

for (i in 1:nrow(sampleSheet)) {
  name = paste(wd,"QC/pseudo_image/","pseudoimage",i,".jpg",sep="")
  jpeg(name)
  image(Pset,which=i,main=ph@data$sample[i])
  dev.off()
}

rm(Pset)

# intensity histogram #
dir.create("./QC/intensity_histo")

for (i in 1:nrow(sampleSheet)) {
  name = paste(wd,"QC/intensity_histo/","histogram",i,".jpg",sep="")
  jpeg(name)
  hist(data[,i],lwd=2,which='pm',ylab='Density',xlab='Log2 intensities',main=ph@data$sample[i])
  dev.off()
}

color <- c()

for (i in 1:nrow(sampleSheet)){
  if(sampleSheet$status[i] == "NL"){
    color[i] <- "green"
  } else {
    color[i] <- "red"
  }
}

jpeg(paste(wd,"QC/intensity_histo/","histogram_all.jpg",sep=""))
hist(data[,1:nrow(sampleSheet)],lwd=2,which='pm',col=color,ylab='Density',xlab='Log2 intensities',main='Histogram of raw data')
dev.off()

# box plot #
dir.create("./QC/boxplot")

name = paste(wd,"QC/boxplot/boxplot.jpg",sep="")
jpeg(name)
boxplot(data,which='pm',col='red',names=ph@data$sample) 
dev.off()

name = paste(wd,"QC/boxplot/boxplotnorm.jpg",sep="")
jpeg(name)
boxplot(data,which='pm',col='red',names=ph@data$sample) 
dev.off()

# MA plot #
dir.create("./QC/MAplot")
dir.create("./QC/MAplotnorm")

for (i in 1:nrow(sampleSheet)) {
  name = paste(wd,"QC/MAplot/","MAplot",i,".jpg",sep="")
  jpeg(name)
  MAplot(data,which=i)
  dev.off()
}

for (i in 1:nrow(sampleSheet)) {
  name = paste(wd,"QC/MAplotnorm/","MAplotnorm",i,".jpg",sep="")
  jpeg(name)
  MAplot(data.rma,which=i)
  dev.off()
}

## Measures quality ## The easiest way will be to create a report

# recuperation of the quality measures #
data.qc <- qc(data)

# average background intensity #
avbg(data.qc)

# scale factors #
sfs(data.qc)

# Call rate # i will extend this step to automatically remove sample with a call rate < 20%
percent.present(data.qc)

# 3'/5' ratio #
ratios(data.qc)

##### Data normalization #####

my.affinity.info <- compute.affinities.local(data)
data.gcrma <- gcrma(data,fast=TRUE,affinity.info = my.affinity.info)
