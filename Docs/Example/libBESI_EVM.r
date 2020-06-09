#!/usr/bin/env Rscript


#==============================================================================
#     libBESI_EVM.r
#     Copyright (C) 2020  Scott P Morton (spm3c at mtmail.mtsu.edu)
# 
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
# 
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
# 
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <http://www.gnu.org/licenses/>.
# 
#==============================================================================


# BESI/EVM analysis tool library 
#   author: Scott P. Morton
#   Center for Computational Science
#   Middle Tennessee State University
#   12/30/2018   

#### required external software ####
#
# MAFFT script is coded against version 7.273 from
# http://mafft.cbrc.jp/alignment/software/
#
# RAxml script is coded against version 8.2.11 from
# https://github.com/stamatak/standard-RAxML
#
# WebLogo 3.6.0
# https://github.com/WebLogo/weblogo
#
# ImageMagick binaries
#
# VMD
#
#### ####

# Load required libraries
library(tools)
library(rgl)
library(colorspace)
#library(colorRamps)
library(fields)
library(ape)
require("seqinr")
library(ips)
library(phytools)
library(SDMTools)
library(plotrix)

# Used for EVM to align a common residue structure across all varients
#### Define HXB2 initial data ####
HXB2CG <-paste(">HXB2\nMRVKEKYQHLWRWGWRWGTMLLGMLMICSATEKLWVTVYYGVPVWKEATTTLFCASDAKAYDTEVHNVWATHA",
               "CVPTDPNPQEVVLVNVTENFNMWKNDMVEQMHEDIISLWDQSLKPCVKLTPLCVSLKCTDLKNDTNTNSSSGRMIMEKG",
               "EIKNCSFNISTSIRGKVQKEYAFFYKLDIIPIDNDTTSYKLTSCNTSVITQACPKVSFEPIPIHYCAPAGFAILKCNNK",
               "TFNGTGPCTNVSTVQCTHGIRPVVSTQLLLNGSLAEEEVVIRSVNFTDNAKTIIVQLNTSVEINCTRPNNNTRKRIRIQ",
               "RGPGRAFVTIGKIGNMRQAHCNISRAKWNNTLKQIASKLREQFGNNKTIIFKQSSGGDPEIVTHSFNCGGEFFYCNSTQ",
               "LFNSTWFNSTWSTEGSNNTEGSDTITLPCRIKQIINMWQKVGKAMYAPPISGQIRCSSNITGLLLTRDGGNSNNESEIF",
               "RPGGGDMRDNWRSELYKYKVVKIEPLGVAPTKAKRRVVQREKR", sep="")
# HXB2 vloop indices (inclusive)
VLOOPS <- list(c(131,157),c(157,196),c(296,331),c(385,418),c(460,470))

#### Define Residue conversion tables ####
# because seqinr::a and seqinr::aaa are case sensitive and mixed case and I use all caps

# Returns the single char representation of an amino acid using the three char
Res1 = c('ASX', 'CYS', 'ASP', 'SER', 'GLN', 'LYS',
         'ILE', 'PRO', 'THR', 'PHE', 'ASN', 
         'GLY', 'HIS', 'LEU', 'ARG', 'TRP', 
         'ALA', 'VAL', 'GLU', 'TYR', 'MET')
names(Res1) <- c('B', 'C', 'D', 'S', 'Q', 'K',
                 'I', 'P', 'T', 'F', 'N', 
                 'G', 'H', 'L', 'R', 'W', 
                 'A', 'V', 'E', 'Y', 'M')
# Returns the three char representation of an amino acid from the single char
Res3 = c('B','C', 'D', 'S', 'Q', 'K',
         'I', 'P', 'T', 'F', 'N', 
         'G', 'H', 'L', 'R', 'W', 
         'A', 'V', 'E', 'Y', 'M')
names(Res3) <- c('ASX', 'CYS', 'ASP', 'SER', 'GLN', 'LYS',
                 'ILE', 'PRO', 'THR', 'PHE', 'ASN', 
                 'GLY', 'HIS', 'LEU', 'ARG', 'TRP', 
                 'ALA', 'VAL', 'GLU', 'TYR', 'MET')

#### Define required Data Structures ####
# Generate a global target directory (subfolder where all data is stored and written)
targetDir = sprintf("%s/Analysis",targetBase)

# Get types table and check for required columns
typesFile <- sprintf("%s/Types/types.txt",targetDir)
myTypes <- read.table(typesFile,stringsAsFactors = FALSE)
if(length(myTypes) < 5){
  for(i in length(myTypes):5){
    switch (i,
      myTypes[,2] <- "A",
      myTypes[,3] <- "TF",
      myTypes[,4] <- "R",
      myTypes[,5] <- "1"
    )
  }
}

if (doSens){
  sensFile <- sprintf("%s/Types/sensitivity.txt",targetDir)
  mySensitivity <- read.table(sensFile,stringsAsFactors = FALSE)
}else{
  mySensitivity <- list()
}

# Set the global sequence list
seqList <- myTypes$V1
# Set the global sequence list length
num_Sequences <- length(seqList)
# Determine unique clades in sequences
myClades <- unique(myTypes[,2])
# Determine unique subclasses in sequences
mySubClasses <- unique(myTypes[,3])
# Detemrine unique roles D/R
myDnR <- unique(myTypes[,4])
# Determine unique Groupings
myGroupings <- unique(myTypes[,5])

# Build pH scale and subsets
pH = seq(3,9,0.1)
pH.subsetLow <- which(as.double(pH)<4.5 & as.double(pH)>=3.5)
pH.subsetHigh <- which(as.double(pH)<8.0 & as.double(pH)>=7.0)

densityFiles <- list("bound-unbound-total.mean","bound-total.mean","unbound-total.mean")
listNames <- list("bound-unbound","bound","unbound")
ratioFiles <- list("bound-unbound-total.ratio","bound-total.ratio","unbound-total.ratio")
densityTitles <- list(expression(paste("Mean ",Delta," Electrostatic Density(bound-unbound) vs. pH Level")),
                      expression(paste("Mean Electrostatic Density(bound) vs. pH Level")),
                      expression(paste("Mean Electrostatic Density(unbound) vs. pH Level")))

ylabels <- list(expression(paste(Delta," ESP (kT/e)")),
                expression(paste("ESP (kT/e)")),
                expression(paste("ESP (kT/e)")))

ratioTitles <- list(expression(paste("Ratio Electrostatic Density(bound-unbound) vs. pH Level")),
                    expression(paste("Ratio Electrostatic Density(bound) vs. pH Level")),
                    expression(paste("Ratio Electrostatic Density(unbound) vs. pH Level")))

sensitivtyTitles <- list(expression(paste("Mean Sensitivity as the Difference Between pH 4.5-6 and pH 7-8.5 (bound-unbound)")),
                         expression(paste("Mean Sensitivity as the Difference Between pH 4.5-6 and pH 7-8.5 (bound)")),
                         expression(paste("Mean Sensitivity as the Difference Between pH 4.5-6 and pH 7-8.5 (unbound)")))

# variables associated with BESI  
all_gp120_PCA <- list(list(),list(),list())
names(all_gp120_PCA) <- listNames
  # Returns a contrast acceptable color palette for donor representations
DonorColors <- colorRampPalette(c("#0090ff","lightgreen"))
BESI_Scores <- NULL
# variables associated with Fingerprints
all_gp120_Density <- list(list(),list(),list())
names(all_gp120_Density) <- listNames
# Variable associated with Binding Energies
BEscale <- c(NULL,NULL) # min/max
BEallBnd_raw <- list()
BEallUBnd_raw <- list()
BEallBnd_means <- list()
BEallUBnd_means <- list()
ComplexTypes <- NULL
# Variables associated with Residue Analaysis
SeqAlignmentResidue <- list()
Residue_Data_Bnd_raw <- list()
Residue_Data_UnBnd_raw <- list()
Residue_Data_Bnd_median <- list()
Residue_Data_UnBnd_median <- list()
Residue_Variance_Data_All <- array()
Residue_Variance_Data_Clades <- list()
Residue_Variance_Data_Classes <- list()
VarianceStandardDeviation <- NULL
VarianceCutoff <- NULL

# Print Titles and Subtitles??
TnS <- FALSE

# Do inline sequence names in BESI Score Bar Graph? TRUE /FALSE
# Does not work well with large lists of sequences or sequences with long names
ISN <- FALSE

#### Define Functions ####
  # methods of global structure,var,object manipulation
  # assign("a", "new", envir = .GlobalEnv)
  # myobject <- get(the_global_object,envir = .GlobalEnv)

AlignHXB2 <- function() { # Creates two files HXB2-'sequence' .fasta and .ali
  initLoopDataFiles()
  HXBConversionFile <- sprintf("%s/Datasets/HXB2_EVM_Selections.txt",targetDir)
  cat("",file = HXBConversionFile,append = FALSE)
  tmpSeq <- ""
  EVMList <- readLines(sprintf("%s/Datasets/Sequence_EVM_Selections.txt",targetDir))
  for (D in 1:length(Dirs)) {
    files <- dir(path = sprintf("%s/seqfiles",Dirs[[D]]),pattern="^.*\\.seq$",full.names = TRUE)
    for(z in 1:length(files)){
      tmp <- strsplit(files[z],'/')
      tmpSeq <- readLines(files[z], n=1,warn = FALSE)
      Sqnc <- tmp[[1]][length(tmp[[1]])]
      seqName <- strsplit(Sqnc[1],".seq")
      outFile <- sprintf("%s/HXB2/HXB2-%s.fasta",targetDir,seqName)
      outFileAligned <- sprintf("%s/HXB2/HXB2-%s.ali",targetDir,seqName)
      output <- sprintf("%s\n>%s\n%s",HXB2CG,seqName,tmpSeq)
      cat(output,file = outFile)
      cmdLine <- sprintf("mafft-einsi --quiet --op 2.0 %s > %s",outFile,outFileAligned)
      system(cmdLine)
      ALI <- read.fasta(outFileAligned, as.string = TRUE,seqtype = "AA",set.attributes = FALSE)
      HXB2 <- ALI[[1]]
      thisSeq <- ALI[[2]]
      resIDHXB <- vector()
      resNumH <- 1
      seqIdx <- 1
      i <- 1
      while(i <= stringi::stri_length(HXB2[[1]])){
        resEx <- 97
        while (stringi::stri_sub(HXB2,i,i) == "-") {
          resIDHXB <- append(resIDHXB,sprintf("%s%s",resNumH,intToUtf8(resEx)),after = length(resIDHXB))
          resEx <- resEx +1
          i <- i + 1
        }
        if (i > stringi::stri_length(HXB2[[1]])) {
          break()
        }
        while(stringi::stri_sub(thisSeq,i,i) == "-"){
          resNumH <- resNumH + 1
          i <- i + 1
          if (i > stringi::stri_length(HXB2[[1]])) {
            break()
          }
        }
        resIDHXB <- append(resIDHXB,sprintf("%s",resNumH),after = length(resIDHXB))
        resNumH <- resNumH + 1
        i <- i + 1
      }
      cat(resIDHXB,file= sprintf("%s/HXB2/%s.hxb2cg",targetDir,seqName))
      sels <- EVMList[which(EVMList == seqName)+1]
      mysels <- strsplit(sels," ")
      #HXBConversion[[seqName[[1]]]] <- resIDHXB[as.integer(mysels[[1]])]
      cat(sprintf("%s\n",seqName[[1]]),file = HXBConversionFile,append = TRUE)
      cat(resIDHXB[as.integer(mysels[[1]])],file = HXBConversionFile,append = TRUE)
      cat("\n",file = HXBConversionFile,append = TRUE)
      #WriteVLoopSeq(seqName[[1]],HXB2,thisSeq,tmpSeq)
    }
  }
  #plot_VloopData()
}
AlignMAFFT <- function(){
  cmdLine <- sprintf("mafft-linsi %s/seqfiles/all_seqs.fasta > %s/seqfiles/all_seqs.ali",targetDir,targetDir)
  system(cmdLine)
  for(i in SeriesFix(1,length(myGroupings))){
    srcFile <- sprintf("%s/seqfiles/group_%s_seqs.fasta",targetDir,myGroupings[i])
    dstFile <- sprintf("%s/seqfiles/group_%s_seqs.ali",targetDir,myGroupings[i])
    cmdLine <- sprintf("mafft-linsi %s > %s",srcFile,dstFile)
    system(cmdLine)
  }
}
buildDirectory <- function(){
  targetDirs <- list("Additional_Plots", # contains any additional plots such as EVM screeplot
                     "BE", # Contains binding energies plots and data
                     "BESI", # Contains BESI related plots
                     "PCA", # Contains all PCA dimensionality reduction plots (First 2 PC's)
                     "Datasets", # All correlated data saved in various formats
                     "EVM_images", # EVM images
                     "FingerPrints", # Electrophoretic Fingerprints (Stieh et al.)
                     "HXB2", # HXB2 alignment data
                     "Logos", # Logos plot data
                     "Residues", # Residue plots
                     "seqfiles", # sequence files
                     "VLoops", # Vloop plots and data
                     "Types"  # Place your types file here
  )
  for (i in 1:length(targetDirs)) {
    TDir <- sprintf("%s/Analysis/%s",targetBase,targetDirs[[i]])
    dir.create(TDir,'recursive'=T)
  }
}
buildFastas <- function() {
  targetFile <- sprintf("%s/seqfiles/all_seqs.fasta",targetDir)
  cat("",file = targetFile,append = FALSE)
  for (D in 1:length(Dirs)) {
    setwd(Dirs[[D]])
    # Read in the source sequences
    seqList.current <- list.files(path = "seqfiles/.", pattern = "^.*\\.seq$")
    seqListAB <- list.files(path = "ABseqfiles/.", pattern = "^.*\\.seq$")
    for (i in 1:length(seqList.current)) {
      # Trim the sequence filename
      thisSeq <- file_path_sans_ext(seqList.current[i])
      thisData <- stringi::stri_read_lines(fname = sprintf("seqfiles/%s.seq",thisSeq))
      cat(sprintf(">%s\n",thisSeq),file = targetFile,append = TRUE)
      cat(sprintf("%s\n",thisData),file = targetFile,append = TRUE)
    }
  }
  allFasta <- read.fasta(targetFile, as.string = TRUE,seqtype = "AA",set.attributes = FALSE)
  for(i in SeriesFix(1,length(myGroupings))){
    targetFile <- sprintf("%s/seqfiles/group_%s_seqs.fasta",targetDir,myGroupings[i])
    write.fasta(allFasta[myTypes[,1][which(myTypes[,5] == myGroupings[i])]],myTypes[,1][which(myTypes[,5] == myGroupings[i])],targetFile)
  }
}
checkAllNaN <- function(thisData){
  if(min(thisData,na.rm = TRUE) == Inf){
    thisData[which(thisData == "NaN")] <- 0
  }
  return(thisData)
}
doCSA <- function(eval,control){
  return(crossprod(eval,control)/sqrt(crossprod(eval) * crossprod(control)))
}
generateVarianceStats <- function(){
  sdeviation <- get("VarianceStandardDeviation", envir = .GlobalEnv)
  varCutoff <- get("VarianceCutoff", envir = .GlobalEnv)
  sdeviation <- sd(Residue_Variance_Data_All);
  varCutoff <- ceiling(sdeviation/2);
  while (TRUE){
    selections = which(Residue_Variance_Data_All > varCutoff)
    for(i in 1:(length(SeqAlignmentResidue))){
      tmp <- strsplit(SeqAlignmentResidue[[i]],"")
      seqName = labels(SeqAlignmentResidue[i])
      outdata <- tmp[[1]][selections]
      if ('-' %in% outdata){
        varCutoff = varCutoff + 1;
        break;
      }
    }
    if(i == length(SeqAlignmentResidue)) {
      break
    }
  }
  percVar <- (sum(Residue_Variance_Data_All[selections])/sum(Residue_Variance_Data_All))*100;
  percRes <- (length(outdata)/length(tmp[[1]]))*100;
  info <- sprintf('%-30s= %3.1f\n%-30s= %3.1f\n%-30s= %3.1f\n%-30s= %3.1f\n%-30s= %3.1f\n%-30s= %3.1f\n'
                  ,'Standard Deviation',sdeviation,
                  '1/2 Standard Deviation',sdeviation/2,
                  'Number of Selected Residues',length(outdata),
                  'Variance cutoff selected',varCutoff,
                  '% of variance selected',percVar,
                  '% of residues selected',percRes)
  fileName <- sprintf("%s/Datasets/EVM_Statistics.txt",targetDir)
  cat('',file = fileName,append = FALSE)
  cat(info,file = fileName,sep = '',append = TRUE)
  cat('\n',file = fileName,append = TRUE)
  
  SPL_Data = sprintf("--stacks-per-line %d",length(outdata))
  
  assign("SPL",SPL_Data, envir = .GlobalEnv)
  assign("VarianceStandardDeviation",sdeviation, envir = .GlobalEnv)
  assign("VarianceCutoff",varCutoff, envir = .GlobalEnv)
}
getBEGraphScale <- function(thisBEData){
  ul <- get("BEscale",envir = .GlobalEnv)
  if(is.null(ul)){
    ul <- c(100000,-100000)
  }
  for (i in 1:length(thisBEData)){
    if(max(thisBEData[[i]]) > ul[2]){
      ul[2] <- max(thisBEData[[i]])
    }
    if(min(thisBEData[[i]]) < ul[1]){
      ul[1] <- min(thisBEData[[i]])
    }
  }
  assign("BEscale",ul , envir = .GlobalEnv)
}
inferTree <- function(){
  setwd(sprintf("%s/seqfiles",targetDir))
  cmdLine <- sprintf("%s -m PROTCATHIVW -s all_seqs.ali -p 12345 -n BESI",raxmlBin)
  system(cmdLine)
  for(i in SeriesFix(1,length(myGroupings))){
    srcFile <- sprintf("%s/seqfiles/group_%s_seqs.ali",targetDir,myGroupings[i])
    cmdLine <- sprintf("%s -m PROTCATHIVW -s %s -p 12345 -n %s",raxmlBin,srcFile,myGroupings[i])
    system(cmdLine)
  }
  setwd(targetBase)
}
initLoopDataFiles <- function() {
  for(i in 1:5){
    fileLoopData <- sprintf("%s/Datasets/V%s-length_data.dat",targetDir,i)
    cat("",file = fileLoopData,append = FALSE)
  }
}
LoadChargeData_EFP_PCA <- function(){
  # Loads charge data, generates EFP graphs, PCA graphs and stores PCA data
  gp120sPCA <- get("all_gp120_PCA",envir = .GlobalEnv)
  gp120sDns <- get("all_gp120_Density",envir = .GlobalEnv)
  for (D in 1:length(Dirs)) {
    setwd(Dirs[[D]])
    # Read in the source sequences
    seqList.current <- list.files(path = "seqfiles/.", pattern = ".seq")
    
    for (i in 1:length(seqList.current)) {
      #idx <- idx + 1
      # Trim the sequence filename
      thisSeq <- file_path_sans_ext(seqList.current[i])
      
      for (j in 1:length(densityFiles)) {
        # Read in data and plot it!
        densityFile <- sprintf("Structures/%s/results/%s",thisSeq,densityFiles[j])
        densityData <- read.table(densityFile)
        densityData <- `colnames<-`(densityData,pH)
        plot_Fingerprint(targetDir,thisSeq,densityData,strsplit(densityFiles[[j]],'.mean')[[1]],j)
        # Populate the data structure
        gp120sDns[[listNames[[j]]]][[thisSeq]] <- densityData
        # Calculate PCA and populate the data structure
        PCAData <- prcomp(densityData,retx=TRUE, center=TRUE, scale=TRUE)
        gp120sPCA[[listNames[[j]]]][[thisSeq]] <- PCAData
      }
    }
  }
  setwd(Dirs[[1]])
  assign("all_gp120_PCA",gp120sPCA , envir = .GlobalEnv)
  assign("all_gp120_Density",gp120sDns , envir = .GlobalEnv)
  # Save datasets for later use
  saveRDS(gp120sPCA,file = sprintf("%s/Datasets/all_gp120_PCA.bin",targetDir))
  saveRDS(gp120sDns,file = sprintf("%s/Datasets/all_gp120_Density.bin",targetDir))
  PerformPCAReductionDns()
}
loadBEData_BE_CA <- function(){
  # Loads BE data, graphs individual BE plots, writes comparitive analysis data to disk
  BEallBndr <- get("BEallBnd_raw",envir = .GlobalEnv)
  BEallUBndr <- get("BEallUBnd_raw",envir = .GlobalEnv)
  BEallBndm <- get("BEallBnd_means",envir = .GlobalEnv)
  BEallUBndm <- get("BEallUBnd_means",envir = .GlobalEnv)
  # Setup the output files for comparitive analysis
  cmplxTypes <- NULL
  fileNameBnd <- sprintf("%s/BE/_Bound_pH%s_pH%s_Data.csv",targetDir, myCols[1], myCols[2])
  fileNameUBnd <- sprintf("%s/BE/_UBound_pH%s_pH%s_Data.csv",targetDir, myCols[1], myCols[2])
  cat("All units of measure in kJ/mol\n",file = fileNameUBnd,append = FALSE)
  cat("All units of measure in kJ/mol\n",file = fileNameBnd,append = FALSE)
  for (D in 1:length(Dirs)) {
    setwd(Dirs[[D]])
    # Read in the source sequences
    seqList <- list.files(path = "seqfiles/.", pattern = "^.*\\.seq$")
    seqListAB <- list.files(path = "ABseqfiles/.", pattern = "^.*\\.seq$")
    # Format the the tables for the comparitive analysis
    cat(",",paste(file_path_sans_ext(seqList),collapse=",,"),"\n",file = fileNameUBnd,append = TRUE)
    cat(",",paste(file_path_sans_ext(seqList),collapse=",,"),"\n",file = fileNameBnd,append = TRUE)
    for (x in 1:length(seqList)) {
      cat(sprintf(",pH %s,pH %s",myCols[1],myCols[2]),file = fileNameUBnd,append = TRUE)
      cat(sprintf(",pH %s,pH %s",myCols[1],myCols[2]),file = fileNameBnd,append = TRUE)
    }
    for (j in 1:length(seqListAB)) {
      ABSeq <- file_path_sans_ext(seqListAB[j])
      cat('\n',file = fileNameBnd,append = TRUE)
      cat('\n',file = fileNameUBnd,append = TRUE)
      
      # output the complex at hand for the comparitive analysis
      cat(sprintf("%s",ABSeq),file = fileNameBnd,append = TRUE)
      cat(sprintf("%s",ABSeq),file = fileNameUBnd,append = TRUE)
      
      for (i in 1:length(seqList)) {
        Seq <- file_path_sans_ext(seqList[i])
        thisInfo <- myTypes[which(myTypes[,1] == Seq),]
        thisComplex <- sprintf("Complex-%s___%s",Seq,ABSeq)
        cmplxTypes <- rbind(cmplxTypes,cbind(thisComplex[[1]],thisInfo$V2,thisInfo$V3))
        if(dir.exists(sprintf("%s/Structures/%s",Dirs[[D]],thisComplex))){
          
          # Read in data
          BEFileBnd <- sprintf("Structures/%s/results/bound_BEnergies.dat",thisComplex)
          BEDataBnd <- read.table(BEFileBnd)
          BEFileUBnd <- sprintf("Structures/%s/results/unbound_BEnergies.dat",thisComplex)
          BEDataUBnd <- read.table(BEFileUBnd)
          
          # Add column names and calculate the mean
          BEDataBnd <- `colnames<-`(BEDataBnd,pH)
          BEDataUBnd <- `colnames<-`(BEDataUBnd,pH)
          BEbndMeans <- colMeans(BEDataBnd)
          BEunbndMeans <- colMeans(BEDataUBnd)
          
          # Accumulate the raw and means data 
          BEallBndr[[thisComplex]] <- BEDataBnd
          BEallUBndr[[thisComplex]] <- BEDataUBnd
          BEallBndm[[thisComplex]] <- BEbndMeans
          BEallUBndm[[thisComplex]] <- BEunbndMeans
          
          # Colaborative analysis output (The 'Josh' Spreadsheet)
          outdataBnd = BEbndMeans[myCols]
          outdataUBnd = BEunbndMeans[myCols]
          cat(sprintf(",%.3f,%.3f",outdataBnd[1],outdataBnd[2]),file = fileNameBnd,append = TRUE)
          cat(sprintf(",%.3f,%.3f",outdataUBnd[1],outdataUBnd[2]),file = fileNameUBnd,append = TRUE)
        }
      }
    }
  }
  setwd(Dirs[[1]])
  assign("BEallBnd_raw",BEallBndr,envir = .GlobalEnv)
  assign("BEallUBnd_raw",BEallUBndr,envir = .GlobalEnv)
  assign("BEallBnd_means",BEallBndm,envir = .GlobalEnv)
  assign("BEallUBnd_means",BEallUBndm,envir = .GlobalEnv)
  assign("ComplexTypes",cmplxTypes,envir = .GlobalEnv)
  saveRDS(BEallBndr,file = sprintf("%s/Datasets/BE_all_Bnd_raw.bin",targetDir))
  saveRDS(BEallUBndr,file = sprintf("%s/Datasets/BE_all_UBnd_raw.bin",targetDir))
  saveRDS(BEallBndm,file = sprintf("%s/Datasets/BE_all_Bnd_means.bin",targetDir))
  saveRDS(BEallUBndm,file = sprintf("%s/Datasets/BE_all_UBnd_means.bin",targetDir))
  saveRDS(cmplxTypes,file = sprintf("%s/Datasets/ComplexTypes.bin",targetDir))
  
  if(uniformBEScale){
    getBEGraphScale(BEallBndm)
    getBEGraphScale(BEallUBndm)
  }
  for(i in 1:length(BEallUBndm)){
    thisComplex <- labels(BEallUBndm[i])
    # A common subtitle for this complex
    thisSubTitle <- sprintf("%s",labels(BEallUBndm[i]))
    
    thisPlot <- sprintf("%s/BE/%s_BE_Bnd.eps",targetDir,thisComplex)
    thisTitle <- "Binding Energies (Bound)"
    plotBindingEnergies(BEallBndm[[i]],thisPlot,thisTitle,thisSubTitle,BEscale)
    
    thisPlot <- sprintf("%s/BE/%s_BE_UBnd.eps",targetDir,thisComplex)
    thisTitle <- "Binding Energies (UnBound)"
    plotBindingEnergies(BEallUBndm[[i]],thisPlot,thisTitle,thisSubTitle,BEscale)
    
    thisPlot <- sprintf("%s/BE/%s_BE_Bnd-UBnd.eps",targetDir,thisComplex)
    thisTitle <- "Binding Energies (Bound - Unbound)"
    plotBindingEnergies(BEallBndm[[i]] - BEallUBndm[[i]],thisPlot,thisTitle,thisSubTitle,BEscale)
  }
}
Load_RDS_Data <- function(){
  gp120sPCA <- readRDS(file = sprintf("%s/Datasets/all_gp120_PCA.bin",targetDir))
  gp120sDns <- readRDS(file = sprintf("%s/Datasets/all_gp120_Density.bin",targetDir))
  ResBnd <- readRDS(file=sprintf('%s/Datasets/Residue_Data_Bnd_raw.bin',targetDir))
  ResUBnd <- readRDS(file=sprintf('%s/Datasets/Residue_Data_UnBnd_raw.bin',targetDir))
  mySeqAli <- readRDS(file=sprintf('%s/Datasets/SeqAlignmentResidue.bin',targetDir))
  varAll <- readRDS(file=sprintf('%s/Datasets/Residue_Variance_Data_All.bin',targetDir))
  ResUBnd <- readRDS(file = sprintf("%s/Datasets/Residue_Data_UnBnd_median.bin",targetDir))
  ResBnd <- readRDS(file = sprintf("%s/Datasets/Residue_Data_Bnd_median.bin",targetDir))
  varClades <- readRDS(file=sprintf('%s/Datasets/Residue_Variance_Data_Clades.bin',targetDir))
  varClasses <- readRDS(file=sprintf('%s/Datasets/Residue_Variance_Data_Classes.bin',targetDir))
  BEallBndr <- readRDS(file = sprintf("%s/Datasets/BE_all_Bnd_raw.bin",targetDir))
  BEallUBndr <- readRDS(file = sprintf("%s/Datasets/BE_all_UBnd_raw.bin",targetDir))
  BEallBndm <- readRDS(file = sprintf("%s/Datasets/BE_all_Bnd_means.bin",targetDir))
  BEallUBndm <- readRDS(file = sprintf("%s/Datasets/BE_all_UBnd_means.bin",targetDir))
  cmplxTypes <- readRDS(file = sprintf("%s/Datasets/ComplexTypes.bin",targetDir))
  scores <- read.table(sprintf("%s/Datasets/BESI_scores.txt",targetDir),stringsAsFactors = FALSE)
  
  assign("all_gp120_PCA",gp120sPCA , envir = .GlobalEnv)
  assign("all_gp120_Density",gp120sDns , envir = .GlobalEnv)
  assign("Residue_Data_Bnd_raw",ResBnd, envir = .GlobalEnv)
  assign("Residue_Data_UnBnd_raw",ResUBnd, envir = .GlobalEnv)
  assign("SeqAlignmentResidue", mySeqAli, envir = .GlobalEnv)
  assign("Residue_Variance_Data_Clades",varClades, envir = .GlobalEnv)
  assign("Residue_Variance_Data_Classes",varClasses, envir = .GlobalEnv)
  assign("Residue_Variance_Data_All",varAll, envir = .GlobalEnv)
  assign("Residue_Data_Bnd_median",ResBnd, envir = .GlobalEnv)
  assign("Residue_Data_UnBnd_median",ResUBnd, envir = .GlobalEnv)
  assign("BEallBnd_raw",BEallBndr,envir = .GlobalEnv)
  assign("BEallUBnd_raw",BEallUBndr,envir = .GlobalEnv)
  assign("BEallBnd_means",BEallBndm,envir = .GlobalEnv)
  assign("BEallUBnd_means",BEallUBndm,envir = .GlobalEnv)
  assign("BESI_Scores",scores,envir = .GlobalEnv)
  assign("ComplexTypes",cmplxTypes,envir = .GlobalEnv)
}
loadResidueData <- function(){
  ResBnd <- get("Residue_Data_Bnd_raw", envir = .GlobalEnv)
  ResUBnd <- get("Residue_Data_UnBnd_raw", envir = .GlobalEnv)
  for (D in 1:length(Dirs)) {
    setwd(Dirs[[D]])
    options(show.error.messages = TRUE)
    seqList <- list.files(path = "seqfiles/.", pattern = "^.*\\.seq$")
    for (i in 1:length(seqList)) {
      seqID = strsplit(seqList[i],'.seq')[[1]]
      cat(sprintf("Loading residue data for %s\n", seqID))
      tmp <- strsplit(SeqAlignmentResidue[[seqID]],'')
      resCounter <- 1
      for (j in 1:length(tmp[[1]])) {
        if(tmp[[1]][j] != '-') {
          resName = Res1[[tmp[[1]][j]]]
          ResBnd[[seqID]][[j]] <- read.table(sprintf('Structures/%s/results/residues/%s_%s_bound-total.mean',
                                                           seqID,resCounter,resName))
          ResUBnd[[seqID]][[j]] <- read.table(sprintf('Structures/%s/results/residues/%s_%s_unbound-total.mean',
                                                            seqID,resCounter,resName))
          names(ResBnd[[seqID]][[j]]) <- pH
          names(ResUBnd[[seqID]][[j]]) <- pH
          resCounter <- resCounter + 1
        }
        else {
          ResBnd[[seqID]][[j]] <- rep(NaN,61)
          ResUBnd[[seqID]][[j]] <- rep(NaN,61)
          names(ResBnd[[seqID]][[j]]) <- pH
          names(ResUBnd[[seqID]][[j]]) <- pH
        }
      }
    }
  }
  cat(sprintf("Residue Data Loaded\n"))
  setwd(Dirs[[1]])
  assign("Residue_Data_Bnd_raw",ResBnd, envir = .GlobalEnv)
  assign("Residue_Data_UnBnd_raw",ResUBnd, envir = .GlobalEnv)
  saveRDS(ResBnd,file=sprintf('%s/Datasets/Residue_Data_Bnd_raw.bin',targetDir))
  saveRDS(ResUBnd,file=sprintf('%s/Datasets/Residue_Data_UnBnd_raw.bin',targetDir))
}
loadSequenceAlignments_Res <- function(){
  mySeqAli <- read.fasta(sprintf('%s/seqfiles/all_seqs.ali',targetDir),
                             as.string = TRUE,seqtype = "AA",set.attributes = FALSE)
  assign("SeqAlignmentResidue", mySeqAli, envir = .GlobalEnv)
  saveRDS(mySeqAli,file=sprintf('%s/Datasets/SeqAlignmentResidue.bin',targetDir))
}
PerformBESIAnalysis <- function(){
  control.list <- list()
  new_cs <- matrix(nrow=num_Sequences,ncol = 2)
  cmean.fixed_2 <- matrix(nrow=num_Sequences,ncol = 1)
  labs <- labels(all_gp120_PCA["unbound"][[1]])
  # calculate PCA for control
  control <- prcomp(template.U,retx=TRUE, center=TRUE, scale=TRUE)
  # calculate cosine sim data
  for(i in 1:length(all_gp120_PCA[[3]])){
    PCACS <- all_gp120_PCA[[3]][[i]]
    if(labels(all_gp120_PCA[[3]])[i] != CONTROL_NAME)
    {
      for (x in 1:2) {
        new_cs[i,x] <- doCSA(PCACS$rotation[,x],control$rotation[,x])
        control.list[i] <- labels(all_gp120_PCA[[3]])[i]
      }
    }
    else{
      control.list[i] <- CONTROL_NAME
      new_cs[i,1:2] <- 1
    }
  }
  # Generate the cosine sim
  for(i in 1:length(all_gp120_PCA[[3]])){
    cmean.fixed_2[i] <- mean(abs(new_cs[i,1:2]))
  }
  scores <- cbind(control.list,cmean.fixed_2)
  assign("BESI_Scores",scores,envir = .GlobalEnv)
  # Plot the cosine similarity graph
  plot_BESI()
  # Generate the phylo correlation graph
  if(length(myGroupings) > 1){
    for(i in SeriesFix(1,length(myGroupings))){
      plot_Phylogeny_vs_BESI(myGroupings[i],
                             theTree = read.tree(sprintf("%s/seqfiles/RAxML_bestTree.%s",targetDir,myGroupings[i])))
    }
  }else{
    plot_Phylogeny_vs_BESI('ALL',theTree = read.tree(sprintf("%s/seqfiles/%s",targetDir,"RAxML_bestTree.BESI")))
  }
  write_BESI_Scores()
}
PerformBEAnalysis <- function(){
  prepBEComparitiveFiles()
  loadBEData_BE_CA()
  for (x in 1:length(myClades)) {
    cladeData1 <- ComplexTypes[which(ComplexTypes[,2] == myClades[x]),]
    # if we have TF and CC subclasses make the comparison
    if(length(mySubClasses) == 2){
      selsData1x <- cladeData1[which(cladeData1[,3] == "TF"),]
      selsData2x <- cladeData1[which(cladeData1[,3] == "CC"),]
      TFDataBndx <- colMeans(do.call(rbind,BEallBnd_means[selsData1x[,1]]))
      TFDataUBndx <- colMeans(do.call(rbind,BEallUBnd_means[selsData1x[,1]]))
      CCDataBndx <- colMeans(do.call(rbind,BEallBnd_means[selsData2x[,1]]))
      CCDataUBndx <- colMeans(do.call(rbind,BEallUBnd_means[selsData2x[,1]]))
      thisTitle <- sprintf("Clade %s - Transmitted Founder (TF) & Chronic Control (CC)",myClades[x])
      thisLegend <- c('Subclass TF','Subclass CC')
      # Plot bound
      thisPlot <- sprintf("%s/BE/Clade_%s_TF_vs_CC_BE_Bound.eps",targetDir,myClades[x])
      thisSubTitle <- sprintf("Subclasses of Clade %s Bound Conformation",myClades[x])
      plotBindingEnergiesCompare(TFDataBndx,CCDataBndx,thisPlot,thisTitle,thisSubTitle,thisLegend)
      # Plot unbound
      thisPlot <- sprintf("%s/BE/Clade_%s_TF_vs_CC_BE_UnBound.eps",targetDir,myClades[x])
      thisSubTitle <- sprintf("Subclasses of Clade %s Unbound Conformation",myClades[x])
      plotBindingEnergiesCompare(TFDataBndx,CCDataBndx,thisPlot,thisTitle,thisSubTitle,thisLegend)
    }
    # if we have more than 1 Clade available make the comparison
    for(y in SeriesFix2(x+1,length(myClades))){
      cladeData2 <- ComplexTypes[which(ComplexTypes[,2] == myClades[y]),]
      
      DataBndx <- colMeans(do.call(rbind,BEallBnd_means[cladeData1[,1]]))
      DataUBndx <- colMeans(do.call(rbind,BEallUBnd_means[cladeData1[,1]]))
      DataBndy <- colMeans(do.call(rbind,BEallBnd_means[cladeData2[,1]]))
      DataUBndy <- colMeans(do.call(rbind,BEallUBnd_means[cladeData2[,1]]))
      
      thisTitle <- sprintf("Clades %s & %s",myClades[x],myClades[y])
      thisLegend <- c(sprintf("Clade %s",myClades[x]),sprintf("Clade %s",myClades[y]))
      # Plot bound
      thisPlot <- sprintf("%s/BE/Clades_%s_vs_Clade_%s_Bound.eps",targetDir,myClades[x],myClades[y])
      thisSubTitle <- sprintf("Clades %s & %s Bound Conformation",myClades[x],myClades[y])
      plotBindingEnergiesCompare(DataBndx,DataBndy,thisPlot,thisTitle,thisSubTitle,thisLegend)
      # Plot unbound
      thisPlot <- sprintf("%s/BE/Clades_%s_vs_Clade_%s_UnBound.eps",targetDir,myClades[x],myClades[y])
      thisSubTitle <- sprintf("Clades %s & %s Unbound Conformation",myClades[x],myClades[y])
      plotBindingEnergiesCompare(DataUBndx,DataUBndy,thisPlot,thisTitle,thisSubTitle,thisLegend)
    }
  }
}
PerformEVMAnalysis <- function(){
  loadSequenceAlignments_Res()
  loadResidueData()
  prepResidueData()
  resLength <- length(Residue_Data_UnBnd_raw[[1]])
  pHLength <- length(pH)
  varAll <- array(dim = resLength)
  varClades <- list()
  varClasses <- list()
  for(i in 1:resLength){
    tmp = array()
    for(j in 1:length(Residue_Data_UnBnd_median)){
      thisSet <- Residue_Data_UnBnd_median[[myTypes[[1]][j]]][i,]
      if(length(tmp) == 1){
        tmp <- thisSet
      }
      else {
        tmp <- rbind(tmp,thisSet)
      }
    }
    tmpArray <- as.matrix(tmp)
    All <- colMeans(tmpArray,na.rm = TRUE)
    dir.create(sprintf("%s/Residues/unbound/Res-%s",targetDir,i), recursive = TRUE)
    # becasue I can't plot a 'NaN' if the row is all NaN, yeah the chances are slim but...!
    if(is.nan(All[1])){
      varAll[i] <- 0
    }
    else{
      varAll[i] <- var(All)
    }
    for(x in SeriesFix(1,length(myClades))){
      varClades[[myClades[x]]][i] <- var(colMeans(tmpArray[which(myTypes[,2]==myClades[x]),],na.rm = TRUE))
    }
    for(x in SeriesFix(1,length(mySubClasses))){
      varClasses[[mySubClasses[x]]][i] <- var(colMeans(tmpArray[which(myTypes[,3]==mySubClasses[x]),],na.rm = TRUE))
    }
    plotResidueData(tmpArray,i,resLength)
  }
  
  if(TnS){
    theTitle <- sprintf("Variance All",i)
  } else {
    theTitle <- NULL
  }
  #thisTitle <- sprintf("Variance All",i)
  thisPlot <- sprintf("%s/Residues/Total_Variance_All_UnBnd.eps",targetDir)
  setEPS()
  postscript(thisPlot)
  theseOpts <- par(no.readonly = TRUE)
  par(mar = c(5,5,4,2) + 0.1)
  plot(varAll,type="l",
       ylab=expression(paste("Variance")),xlab="Residue", main=theTitle,cex.lab = 1.4, cex.main = 1.5)
  dev.off()
  par(theseOpts)
  
  assign("Residue_Variance_Data_Clades",varClades, envir = .GlobalEnv)
  assign("Residue_Variance_Data_Classes",varClasses, envir = .GlobalEnv)
  assign("Residue_Variance_Data_All",varAll, envir = .GlobalEnv)
  saveRDS(varClades,file=sprintf('%s/Datasets/Residue_Variance_Data_Clades.bin',targetDir))
  saveRDS(varClasses,file=sprintf('%s/Datasets/Residue_Variance_Data_Classes.bin',targetDir))
  saveRDS(varAll,file=sprintf('%s/Datasets/Residue_Variance_Data_All.bin',targetDir))
  # Generate the statistics of the variance data
  generateVarianceStats()
  prepEVMSelectionData()
  plotVarianceScreePlot()
  prepLogosData()
  for (D in 1:length(Dirs)) {
    cmdLine <- sprintf("vmd -dispdev none -e %s/tools/VMD_gen_EVM_images.tcl -args %s %s",Dirs[D],targetBase,Dirs[D])
    system(cmdLine)
  }
}
PerformEVIAnalysis <- function(){
}
PerformLoopAnalysis <- function(){
  for (D in 1:length(Dirs)) {
    files <- dir(path = sprintf("%s/seqfiles",Dirs[[D]]),pattern="^.*\\.seq$",full.names = TRUE)
    for(z in 1:length(files)){
      tmp <- strsplit(files[z],'/')
      tmpSeq <- readLines(files[z], n=1,warn = FALSE)
      Sqnc <- tmp[[1]][length(tmp[[1]])]
      seqName <- strsplit(Sqnc[1],".seq")
      outFile <- sprintf("%s/HXB2/HXB2-%s.fasta",targetDir,seqName)
      outFileAligned <- sprintf("%s/HXB2/HXB2-%s.ali",targetDir,seqName)
      output <- sprintf("%s\n>%s\n%s",HXB2CG,seqName,tmpSeq)
      cat(output,file = outFile)
      cmdLine <- sprintf("mafft-einsi --quiet --op 2.0 %s > %s",outFile,outFileAligned)
      system(cmdLine)
      ALI <- read.fasta(outFileAligned, as.string = TRUE,seqtype = "AA",set.attributes = FALSE)
      HXB2 <- ALI[[1]]
      thisSeq <- ALI[[2]]
      #cat(resIDHXB,file= sprintf("%s/HXB2/%s.hxb2cg",targetDir,seqName))
      #HXBConversion[[seqName[[1]]]] <- resIDHXB[as.integer(mysels[[1]])]
      WriteVLoopSeq(seqName[[1]],HXB2,thisSeq,tmpSeq)
    }
  }
  plot_VloopData()
}

PerformPCAReductionDns <- function(){
  gp120sPCA <- get("all_gp120_PCA",envir = .GlobalEnv)
  tDir <- sprintf("%s/PCA",targetDir)
  listN <- list("bound","unbound")
  # Loads charge data, generates EFP graphs, PCA graphs and stores PCA data
  for(i in 1:length(all_gp120_PCA[[3]])){
    PCAu <- all_gp120_PCA[[3]][[i]]
    PCAb <- all_gp120_PCA[[2]][[i]]
    thisSeq <- labels(all_gp120_PCA[[3]])[i]
    # Trim the sequence filename
    # file_path_sans_ext(seqList.current[i])
    unbound <- t(t(PCAu$x[,1:2] %*% t(PCAu$rotation[,1:2])) * PCAu$scale + PCAu$center)
    bound <- t(t(PCAb$x[,1:2] %*% t(PCAb$rotation[,1:2])) * PCAb$scale + PCAb$center)
    plot_Fingerprint(tDir,thisSeq,unbound,sprintf("%s_PCA",strsplit(densityFiles[[3]],'.mean')[[1]]),3)
    plot_Fingerprint(tDir,thisSeq,bound,sprintf("%s_PCA",strsplit(densityFiles[[2]],'.mean')[[1]]),2)
    plot_Fingerprint(tDir,thisSeq,bound - unbound,sprintf("%s_PCA",strsplit(densityFiles[[1]],'.mean')[[1]]),1)
  }
}
prepBEComparitiveFiles <- function(){
  fileNameBnd <- sprintf("%s/BE/_Bound_pH%s_pH%s_Data.csv",targetDir, myCols[1], myCols[2])
  fileNameUBnd <- sprintf("%s/BE/_UBound_pH%s_pH%s_Data.csv",targetDir, myCols[1], myCols[2])
  cat("All units of measure in kJ/mol\n",file = fileNameUBnd,append = FALSE)
  cat("All units of measure in kJ/mol\n",file = fileNameBnd,append = FALSE)
  cat(",",paste(file_path_sans_ext(seqList),collapse=",,"),"\n",file = fileNameUBnd,append = TRUE)
  cat(",",paste(file_path_sans_ext(seqList),collapse=",,"),"\n",file = fileNameBnd,append = TRUE)
  for (x in 1:length(seqList)) {
    cat(sprintf(",pH %s,pH %s",myCols[1],myCols[2]),file = fileNameUBnd,append = TRUE)
    cat(sprintf(",pH %s,pH %s",myCols[1],myCols[2]),file = fileNameBnd,append = TRUE)
  }
}
prepEVMSelectionData <- function(){
  fileName <- sprintf("%s/Datasets/Sequence_EVM_Selections.txt",targetDir)
  cat('',file = fileName,append = FALSE)
  thisSubset <- which(Residue_Variance_Data_All > VarianceCutoff)
  for(i in 1:(length(SeqAlignmentResidue))){
    tmp <- strsplit(SeqAlignmentResidue[[i]],"")
    seqName = labels(SeqAlignmentResidue[i])
    outdata <- tmp[[1]][thisSubset]
    tmpMap <- array(dim = length(tmp[[1]]))
    cnt <- 1 # Indexed according to VMD resID and PDB indexing schemes
    for (j in 1:(length(tmp[[1]]))){
      if(tmp[[1]][j] == '-'){
        tmpMap[j] <- '-1'
      }
      else{
        tmpMap[j] <- cnt
        cnt <- cnt + 1
      }
    }
    mapSelection = tmpMap[thisSubset]
    cat(seqName,file = fileName,append = TRUE)
    cat('\n',file = fileName,append = TRUE)
    cat(mapSelection,file = fileName,append = TRUE)
    cat('\n\n',file = fileName,append = TRUE)
  }
}
prepLogosData <- function(){
  fileName <- sprintf("%s/Logos/ubnd_logoalign_all_sequences.txt",targetDir)
  cat('',file = fileName,append = FALSE)
  thisSubset <- which(Residue_Variance_Data_All > VarianceCutoff)
  for(i in 1:length(SeqAlignmentResidue)){
    tmp <- strsplit(SeqAlignmentResidue[[i]],"")
    outdata <- tmp[[1]][thisSubset]
    cat(outdata,file = fileName,sep = '',append = TRUE)
    cat('\n',file = fileName,append = TRUE)
  }
  outFileName <- sprintf("%s/Logos/ubnd_logoalign_all_sequences.pdf",targetDir)
  if(TnS){
    thisTitle <- 'Distribution of Selected Residues (Unbound) (All)'
    cmdLine <- sprintf("weblogo -f %s -o %s -F pdf -t '%s' -X NO --errorbars NO %s %s -P '' -S 6 %s %s",fileName,outFileName,thisTitle,AR,SPL,TFS,NFS)
  } else{
    cmdLine <- sprintf("weblogo -f %s -o %s -F pdf -X NO --errorbars NO %s %s -P '' -S 6 %s %s",fileName,outFileName,AR,SPL,TFS,NFS)
  }
  system(cmdLine)
  
  # This sections has unanswered questions. As the number of involved sequences changes
  #   the number of aligned residues that match the criteria can and will change.
  # Do we express or surpress this information?
  # How does it change the selection map and would the starting cutoff point need adjusted
  #   which is currently 1/2 standard deviation?
  # Clearly this changes the weblogos expressions but what exactly are we trying to express?
  #   Ans - The conservation of residues across the set - Then use the subset from all sequences
  #   Ans - What is different among the set by clade or subclass? - Then use selections from those subsets
  #   Ans - What is different between each individual sequence? - Then each one needs to be analyzed individually
  for(x in SeriesFix(1,length(mySubClasses))){
    subclassSet <- SeqAlignmentResidue[myTypes[which(myTypes[[3]]==mySubClasses[x]),1]]
    fileName <- sprintf("%s/Logos/ubnd_logoalign_subclass_%s.txt",targetDir,mySubClasses[x])
    cat('',file = fileName,append = FALSE)
    #thisSubset <- which(Residue_Variance_Data_Classes[[mySubClasses[x]]] > VarianceCutoff)
    for(i in 1:length(subclassSet)){
      tmp <- strsplit(subclassSet[[i]],"")
      outdata <- tmp[[1]][thisSubset]
      cat(outdata,file = fileName,sep = '',append = TRUE)
      cat('\n',file = fileName,append = TRUE)
    }
    outFileName <- sprintf("%s/Logos/ubnd_logoalign_subclass_%s.pdf",targetDir,mySubClasses[x])
    if(TnS){
      thisTitle <- sprintf("Distribution of Selected Residues (Unbound) (Subclass %s)",mySubClasses[x])
      cmdLine <- sprintf("weblogo -f %s -o %s -F pdf -t '%s' -X NO --errorbars NO %s %s -P '' -S 6 %s %s",fileName,outFileName,thisTitle,AR,SPL,TFS,NFS)
    } else{
      cmdLine <- sprintf("weblogo -f %s -o %s -F pdf -X NO --errorbars NO %s %s -P '' -S 6 %s %s",fileName,outFileName,AR,SPL,TFS,NFS)
    }
    system(cmdLine)
  }
  for (x in SeriesFix(1,length(myClades))) {
    cladeSet <- SeqAlignmentResidue[myTypes[which(myTypes[[2]]==myClades[x]),1]]
    fileName <- sprintf("%s/Logos/ubnd_logoalign_clade_%s.txt",targetDir,myClades[x])
    cat('',file = fileName,append = FALSE)
    #thisSubset <- which(Residue_Variance_Data_Clades[[myClades[x]]] > VarianceCutoff)
    for(i in 1:length(cladeSet)){
      tmp <- strsplit(cladeSet[[i]],"")
      outdata <- tmp[[1]][thisSubset]
      cat(outdata,file = fileName,sep = '',append = TRUE)
      cat('\n',file = fileName,append = TRUE)
      }
    outFileName <- sprintf("%s/Logos/ubnd_logoalign_clade_%s.pdf",targetDir,myClades[x])
    if(TnS){
      thisTitle <- sprintf("Distribution of Selected Residues (Unbound) (Clade %s)",myClades[x])
      cmdLine <- sprintf("weblogo -f %s -o %s -F pdf -t '%s' -X NO --errorbars NO %s %s -P '' -S 6 %s %s",fileName,outFileName,thisTitle,AR,SPL,TFS,NFS)
    } else {
      cmdLine <- sprintf("weblogo -f %s -o %s -F pdf -X NO --errorbars NO %s %s -P '' -S 6 %s %s",fileName,outFileName,AR,SPL,TFS,NFS)
    }
    system(cmdLine)
  }
}
prepResidueData <- function(){
  ResBnd <- get("Residue_Data_Bnd_median", envir = .GlobalEnv)
  ResUBnd <- get("Residue_Data_UnBnd_median", envir = .GlobalEnv)
  resLength <- length(Residue_Data_UnBnd_raw[[1]])
  pHLength <- length(pH)
  bndData <- matrix(nrow = resLength,ncol = pHLength)
  unbndData <- matrix(nrow = resLength,ncol = pHLength)
  cat("Preparing residue data\n")
  for(i in 1:length(Residue_Data_UnBnd_raw)){
    for(j in 1:resLength){
      if(as.matrix(Residue_Data_UnBnd_raw[[i]][[j]])[1] == "NaN"){
        bndData[j,] <- as.matrix(Residue_Data_Bnd_raw[[i]][[j]])
        unbndData[j,] <- as.matrix(Residue_Data_UnBnd_raw[[i]][[j]])
      }
      else{
        for(p in 1:length(pH)){
          bndData[j,p] <- median(as.matrix(Residue_Data_Bnd_raw[[i]][[j]])[,p],na.rm = TRUE)
          unbndData[j,p] <- median(as.matrix(Residue_Data_UnBnd_raw[[i]][[j]])[,p],na.rm = TRUE)
        }
      }
    }
    ResBnd[[myTypes[[1]][i]]] <- bndData
    ResUBnd[[myTypes[[1]][i]]] <- unbndData
    cat(".")
  }
  cat("\nResidue data is prepared\n")
  assign("Residue_Data_Bnd_median",ResBnd, envir = .GlobalEnv)
  assign("Residue_Data_UnBnd_median",ResUBnd, envir = .GlobalEnv)
  # save the compiled data for later use
  saveRDS(ResUBnd,file = sprintf("%s/Datasets/Residue_Data_UnBnd_median.bin",targetDir))
  saveRDS(ResBnd,file = sprintf("%s/Datasets/Residue_Data_Bnd_median.bin",targetDir))
}
PrepareSequenceData <- function(){
  for(i in 1:length(myGroupings)){
    # future functionality??? must of had something in mind  
  }
  # Step 1 build an all inclusive fasta file
  buildFastas()
  # step 2 align with MAFFT
  AlignMAFFT()
  # Step 3 infer a phylogeny tree with RaXML
  inferTree()
  
}
plotBindingEnergies <- function(thisData,thisPlot,thisTitle,thisSubTitle, ul){
  setEPS()
  postscript(thisPlot,height=8.5,width=16)
  theseOpts <- par(no.readonly = TRUE)
  reg <- lm(thisData ~ seq(1:61))
  par(mar = c(8,6,4,2) + 0.1, mgp=c(3,1,0))
  if(TnS){
    theTitle <- thisTitle
    theSubTitle <- thisSubTitle
  } else {
    theTitle <- NULL
    theSubTitle <- NULL
  }
  plot(thisData,type = "l",col=c("blue"),cex=2, ylim = ul,
       cex.lab = 2.5, xlab="pH",xaxt = "n",ylab="Binding Energy (kJ/mol)",
       cex.main = 3, main=theTitle, cex.axis=1.5)
  if(BELineZero){
    abline(h=0,col="red")
  }
  if(BELineTrend){
    abline(reg,col="green")
  }
  axis(1,at=seq(1,61),labels = pH)
  mtext(side=1,theSubTitle,line=6,cex=2.5)
  dev.off()
  par(theseOpts)
}
plotBindingEnergiesCompare <- function(thisData1,thisData2,thisPlot,thisTitle,thisSubTitle,thisLegend){
  BECdensityRanges <- matrix(c(min(thisData1),max(thisData1),min(thisData2),max(thisData2)),nrow=2,byrow = TRUE)
  setEPS()
  postscript(thisPlot,height=8.5,width=16)
  theseOpts <- par(no.readonly = TRUE)
  par(mar = c(8,6,4,2) + 0.1, mgp=c(3,1,0))
  if(TnS){
    theTitle <- thisTitle
    theSubTitle <- thisSubTitle
  } else {
    theTitle <- NULL
    theSubTitle <- NULL
  }
  plot(thisData1,type = "l",col=c("red"),cex=2,
       cex.lab = 2.5, xlab="pH",xaxt = "n",ylab="Binding Energy (kJ/mol)",
       cex.main = 3, main=theTitle, cex.axis=1.5, ylim = c(min(BECdensityRanges[,1]),max(BECdensityRanges[,2])))
  lines(thisData2,type = "l",col = "green")
  axis(1,at=seq(1,61),labels = pH)
  legend("bottomright",
         legend = thisLegend,
         inset = 0.05,
         col = c("red",'green'),
         lty=c(1,1))
  mtext(side=1,theSubTitle,line=6,cex=2.5)
  dev.off()
  par(theseOpts)
}
plotVarianceScreePlot <- function(){
  thisPlot <- sprintf("%s/Additional_Plots/Residue_Variance_Screeplot_Ubnd.eps",targetDir)
  thisLegend = c(sprintf("Variance = %3.0f",VarianceCutoff))
  if(TnS){
    theTitle <- "Residue Variance - All"
  } else {
    theTitle <- NULL
  }
  setEPS()
  postscript(thisPlot,width = 11,height = 11)
  theseOpts <- par(no.readonly = TRUE)
  par(mar = c(5,5,5,3) +.1)
  #par(yaxp = c(0,1750,250))
  par(mgp = c(3,1,0.2))
  #par(yaxs = "i")
  centers <- barplot(sort(Residue_Variance_Data_All,decreasing=TRUE),col="lightblue",cex.lab = 2, cex.main = 3.5,axes=FALSE,
                     main=theTitle, ylim=c(0,1600),ylab="Total Variance", xlab = "Residues",width=1.1)
  abline(h=VarianceCutoff, col = "red")
  legend("topright",legend=thisLegend,col = c("red"),lty = 1,cex=1.8)
  usr <- par("usr")
  par(usr=c(usr[1:2], 0, 1600))
  axis(side = 2, at=seq(0,1600,100),las=1, tick = TRUE, outer = FALSE)
  dev.off()
  par(theseOpts)  
}
plotResidue <- function(thisData,thisPlot,thisTitle,thisSubTitle){
  if(TnS){
    theTitle <- thisTitle
    theSubTitle <- thisSubTitle
  } else {
    theTitle <- NULL
    theSubTitle <- NULL
  }
  setEPS()
  postscript(thisPlot)
  theseOpts <- par(no.readonly = TRUE)
  par(mar = c(5,5,4,2) + 0.1)
  boxplot(thisData,cex=1,notch=TRUE,xlab="pH",names=pH,
          ylab=expression(paste("ESP (kT/e)")),cex.lab = 1.4, cex.main = 1.5, main=theTitle)
  dev.off()
  par(theseOpts)
}
plotResidueCompare <- function(thisData1,thisData2,thisPlot,thisTitle,thisSubTitle){
  if(TnS){
    theTitle <- thisTitle
    theSubTitle <- thisSubTitle
  } else {
    theTitle <- NULL
    theSubTitle <- NULL
  }
  setEPS()
  postscript(thisPlot)
  theseOpts <- par(no.readonly = TRUE)
  par(mar = c(5,5,4,2) + 0.1)
  boxplot(list(thisData1,thisData2),col=c("blue","green"),cex=2,notch=TRUE,
          ylab=expression(paste("ESP (kT/e)")),cex.lab = 1.4, cex.main = 1.5,
          main=theTitle,names=c("TF","CC"),cex.sub = 1.2)
  dev.off()
  par(theseOpts)
}
plotResidueData <- function(tmpArray,resNum,resLength){
  # common subtitle
  thisSubTitle <- sprintf("Residue %s of %s (Aligned Length) ",resNum,resLength)
  # if we have TF and CC variants do the comparison
  if(length(mySubClasses) == 2){
    # Plot the TF seqeunces
    thisPlot <- sprintf("%s/Residues/unbound/Res-%s/Res-%s_TF_vs_CC_All.eps",targetDir,resNum,resNum)
    thisTitle <- "Variance of Transmitted Founder & Cronic Control"
    tmpData1 <- tmpArray[which(myTypes[,3]=="TF"),]
    tmpData2 <- tmpArray[which(myTypes[,3]=="CC"),]
    if(!length(tmpData1[1]) || !length(tmpData2[1])){
      next
    }
    tmpData1 <- checkAllNaN(tmpData1)
    tmpData2 <- checkAllNaN(tmpData2)
    diffData <- checkAllNaN(tmpData1 - tmpData2)
    # Plot the comparison
    plotResidueCompare(tmpData1,tmpData2,thisPlot,thisTitle,thisSubTitle)
    # Plot the difference
    thisPlot <- sprintf("%s/Residues/unbound/Res-%s/Res-%s_Motif_TF-CC_All.eps",targetDir,resNum,resNum)
    thisTitle <- "Variance Transmitted Founder Minus Cronic Control"
    plotResidue(diffData,thisPlot,thisTitle,thisSubTitle)
  }
  for (x in 1:length(myClades)) {
    # Plot individual Clades, Subclasses at this level become a little more complex and provide little value
    thisPlot <- sprintf("%s/Residues/unbound/Res-%s/Res-%s_Motif_Clade_%s.eps",targetDir,resNum,resNum,myClades[x])
    thisTitle <- sprintf("Motif of Residue %s From Clade %s",resNum,myClades[x])
    cladeData1 <- tmpArray[which(myTypes[,2]==myClades[x]),]
    if(!length(cladeData1)){
      next
    }
    cladeData1 <- checkAllNaN(cladeData1)
    plotResidue(cladeData1,thisPlot,thisTitle,thisSubTitle)
    for(y in SeriesFix(x+1,length(myClades))){
      thisPlot <- sprintf("%s/Residues/unbound/Res-%s/Res-%s_Motif_Clade_%s_vs_%s.eps",targetDir,resNum,resNum,myClades[x],myClades[y])
      thisTitle <- sprintf("Motif of Residue %s in Clade %s vs %s",resNum,myClades[x],myClades[y])
      cladeData2 <- myTypes[which(myTypes[,2] == myClades[y]),]
      cladeData2 <- checkAllNaN(cladeData2)
      plotResidueCompare(cladeData1,cladeData2,thisPlot,thisTitle,thisSubTitle)
    }
  }
  # Plot All Clades, Subclasses if more than one Clade present
  if(length(tmpArray) && length(myClades) > 1){
    thisPlot <- sprintf("%s/Residues/unbound/Res-%s/Res-%s_Motif_ALL_%s.eps",targetDir,resNum,resNum,myClades[x])
    thisTitle <- sprintf("Motif of Residue %s From All Clades",resNum)
    tmpArray <- checkAllNaN(tmpArray)
    plotResidue(tmpArray,thisPlot,thisTitle,thisSubTitle)
  }
}
plot_BESI <- function(){
  # Setup the cosine similarity graph
  thisPlot <- sprintf("%s/BESI/Sequence-BESI-Score_PC_columns.eps",targetDir)
  if(TnS){
    thisTitle <- sprintf("%s - BESI Scores",studyName)
    if(doSens){
      theSubTitle <- "Red indicates experimental agreement"
    } else {
      theSubTitle <- "Red indicates score >= 0.8"
    }
  } else {
    theTitle <- NULL
    theSubTitle <- NULL
  }
  setEPS()
  postscript(thisPlot,width = 12,height = 9)
  theseOpts <- par(no.readonly = TRUE)
  par(mar = c(5,5,4,2) + 0.1)
  if(doSens){
    centers <-barplot(as.double(BESI_Scores[,2]),
                      col = ifelse(mySensitivity[,2],"red","lightblue"),
                      ylab="Score", sub=theSubTitle,
                      xlab="Sequence", cex.lab = 2.5, cex.main = 3, main=theTitle,
                      cex.axis = 1.5, axisnames = FALSE, ylim = c(0,1.05),las=2)
    
  } else {
  centers <-barplot(as.double(BESI_Scores[,2]),col = ifelse(BESI_Scores[,2] >= 0.8,"red","lightblue"),
                    ylab="Score", sub=theSubTitle,
                    xlab="Sequence", cex.lab = 2.5, cex.main = 3, main=theTitle,
                    cex.axis = 1.5, axisnames = FALSE, ylim = c(0,1.05),las=2)
  }
  par(srt=90)
  if(ISN){
    text(centers,0.01,labels = BESI_Scores[,1],adj = c(0,0.5),cex = 2)
  }
  par(srt=0)
  #axis(side = 1, at=centers[,1], labels = as.character(1:num_Sequences), las=1, lwd = 0.1)
  dev.off()
  par(theseOpts)
}
plot_Phylogeny_vs_BESI <- function(theGroup,theTree){
  mySensitivity <- get("mySensitivity",envir = .GlobalEnv)
  # thisTitle <- sprintf("", studyName,theGroup)
  # thisTitle <- sprintf("%s-%s (Phylogenetic Tree vs BESI)", studyName,theGroup)
  # theTree <- read.tree(sprintf("%s/seqfiles/%s",targetDir,"RAxML_bestTree.BESI"))

  # Determine a more suitable x.lim for the phylo plot
  L <- array()
  xLim <- NULL
  for(i in 1:length(theTree$tip.label)){
    L[i] <- stringi::stri_length(theTree$tip.label[i])
  }
  if(max(L) > charLimit){
    # Need to get some dimensions of the plot to adjust for long tip lables
    setEPS()
      postscript("/dev/null",width = 24,height = 24)
    theseOpts <- par(no.readonly = TRUE)
    nf <- layout(matrix(c(1,2), 1, 2, byrow = TRUE), widths=c(0.6,3.2), heights=c(0.1,1.9))
    #nf <- layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE), widths=c(0.6,3.2), heights=c(0.1,1.9))
    #par(mar = c(0.0,0.5,4.0,0.5) + 0.1)
    #plot(1, ylim=c(0,0), xlim=c(0,3), type="n", axes=F, xlab="", ylab="",main=thisTitle,cex.main=4.5)
    par(mar = c(1,1,1,1) + 0.1)
    plot(1, ylim=c(1,10), xlim=c(0,3), type="n", axes=F, xlab="", ylab="")
    par(mar = c(1.5,0,0,35) + 0.1,xpd=NA)
    pinfo <- plot(theTree, underscore = TRUE, tip.color = "white",
                  use.edge.length = TRUE, root.edge = TRUE,no.margin = FALSE,
                  cex.lab = 3,cex = 3.5,x.lim = xLim)
    par(theseOpts)
    dev.off()
    par(theseOpts)
    
    xLim <- pinfo$x.lim
    #xLim <- c(0.0,((max(L) - charLimit) * xLimFactor + pinfo$x.lim[2])*0.37)
  }
  
  R <- rev(heat.colors(100))
  D <- rev(DonorColors(100))
  myColors <- list(R,D)
  names(myColors) <- list("R", "D")
  twoBars <- TRUE
  if(length(myDnR) == 1){
    twoBars <- FALSE
  }
  
  thisPlot <- sprintf("%s/BESI/%s_%s_Phylo_vs_BESI.eps",targetDir,studyName,theGroup)
  setEPS()
  postscript(thisPlot,width = 24,height = 24)
  theseOpts <- par(no.readonly = TRUE)
  nf <- layout(matrix(c(1,2), 1, 2, byrow = TRUE), widths=c(0.6,3.2), heights=c(0.1,1.9))
  #nf <- layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE), widths=c(0.6,3.2), heights=c(0.1,1.9))
  #par(mar = c(0.0,0.5,4.0,0.5) + 0.1)
  #plot(1, ylim=c(0,0), xlim=c(0,3), type="n", axes=F, xlab="", ylab="",main=thisTitle,cex.main=4.5)
  par(mar = c(1,1,1,1) + 0.1)
  plot(1, ylim=c(1,10), xlim=c(0,3), type="n", axes=F, xlab="", ylab="")
  par(srt=90,xpd=FALSE)
  if(twoBars){
    gradient.rect(1,1,1.5,6,col = R,gradient = "vertical")
    gradient.rect(1.5,1,2,6,col = D,gradient = "vertical")
    par(srt=90)
    text(x=1.7,y=1.1,pos=4,labels = "Donor",cex = 3.5)
    text(x=1.2,y=1.1,pos=4,labels = "Recipient",cex = 3.5)
  } else {
    gradient.rect(1,1,2.0,6,col = R,gradient = "vertical")
  }
  text(x=0.8,y=3,pos = 3,labels = "Similarity to Control",cex = 3.5)
  par(srt=0,xpd=FALSE)
  text(x=0.6,y=6.0,pos=1,labels = "1.0",cex=3)
  text(x=0.6,y=1.0,pos=3,labels = "0.0",cex=3)
  text(x=2.5,y=6.0,pos=1,labels = "High",cex=3)
  text(x=2.5,y=1.0,pos=3,labels = "Low",cex=3)
  
  par(mar = c(1.5,0,0,35) + 0.1,xpd=NA)
  pinfo <- plot(theTree, underscore = TRUE, tip.color = "white",
                use.edge.length = TRUE, root.edge = TRUE,no.margin = FALSE,
                cex.lab = 3,cex = 3.5,x.lim = xLim
  )
  add.scale.bar(0,0.8,length = NULL, ask = FALSE, cex = 2, font = 3.5, col = "red")
  for(p in 1:length(theTree$tip.label)){
    val <- BESI_Scores[[which(BESI_Scores[,1]==theTree$tip.label[p]),2]]
    tiplabels(theTree$tip.label[p],p,
              bg=myColors[[myTypes[which(myTypes[,1] == theTree$tip.label[p]),4]]][as.integer(as.double(val)*100)],
              adj=0,
              col = "black",
              cex=3.5)
    #tiplabels(pch = 8,p,col = "red",adj=1.4,cex=2)
    if(doSens){
      if (mySensitivity[grep(theTree$tip.label[p],mySensitivity[,1]),2]){
        #tiplabels(pch = 19,tip=p,col = "red",adj=1.4,cex=2)
        tiplabels(NULL,tip = p,adj=0.48,pch=16,col = "red",cex=7)
      }
      legend("topright", c("Experimental Agreement"),
             pch = 16,col = "red", 
             bg = "gray95",cex = 3)
      
    }  
  }
  par(theseOpts)
  dev.off()
}
plot_Fingerprint <- function(tDir,thisSeq,thisData,thisDenseFile,j){
  if(TnS){
    theTitle <- densityTitles[j]
    theSubTitle <- sprintf("%s",thisSeq)
  } else {
    theTitle <- NULL
    theSubTitle <- NULL
  }
  thisPlot <- sprintf("%s/FingerPrints/%s_%s.eps",tDir,thisSeq,thisDenseFile)
  setEPS()
  postscript(thisPlot,height=8.5,width=16)
  theseOpts <- par(no.readonly = TRUE)
  par(mar = c(8,6,4,2) + 0.1, mgp=c(3,1,0))
  
  boxplot(thisData,col=c("lightblue"),cex=2,notch=TRUE,
          cex.lab = 2.5, xlab="pH",ylab=ylabels[j],names=pH,
          cex.main = 3, main=theTitle, cex.axis=1.5,
          ylim = densityRanges[j,],las=1, lheight = 3)
  mtext(side=1,theSubTitle,line=6,cex=2.5)
  dev.off()
  par(theseOpts)
}
plot_VloopData <- function(){
  for (f in 1:5) {
    if(TnS){
      theTitle <- sprintf("Variable Loop %s vs BESI Score",f)
      theSubTitle <- vloopSubTitle
    } else {
      theTitle <- NULL
      theSubTitle <- NULL
    }
    inFile <- sprintf("%s/Datasets/V%s-length_data.dat",targetDir,f)
    X <- read.table(inFile)
    
    #thisTitle <- sprintf("Variable Loop %s vs BESI Score",f)
    thisPlot <- sprintf("%s/VLoops/vloop-%s_vs_BESI.eps",targetDir,f)
    setEPS()
    postscript(thisPlot,width = 17,height = 11)
    theseOpts <- par(no.readonly = TRUE)
    par(mar = c(6,7,3,2) + 0.1)
    
    plot(X[2:3],col = ifelse(X[3][,1] == 1.0,"red","lightblue"),
         ylab="BESI Score",sub=theSubTitle,
         xlab="Number of Residues", cex.sub = 2.5,cex.lab = 3.0, cex.main = 3, main=theTitle,
         cex.axis = 2.0, ylim = c(0,1.05),las=1,pch = 19,cex = 3)
    dev.off()
    par(theseOpts)
    
    #thisTitle <- sprintf("Variable Loop %s vs BESI Score",f)
    thisPlot <- sprintf("%s/VLoops/vloop-%s_vs_BESI_quadrants.eps",targetDir,f)
    setEPS()
    postscript(thisPlot,width = 17,height = 11)
    theseOpts <- par(no.readonly = TRUE)
    par(mar = c(6,7,3,2) + 0.1)
    
    plot(X[2:3],col = ifelse(X[3][,1] == 1.0,"red","lightblue"),ylab = '',
         xlab="Number of Residues", cex.sub = 2.5,cex.lab = 3.0, cex.main = 3, main=theTitle,
         cex.axis = 2.5, ylim = c(0,1.05),las=1,pch = 19,cex = 3.5)
    title(ylab="BESI Score", cex.lab = 3.0, line = 4)
    title(sub=theSubTitle,cex.sub = 2.5, line = 5)
    abline(h=0.5, v=(min(X[2])+((max(X[2])- min(X[2]))*.5)), col = "black")
    dev.off()
    par(theseOpts)
  }
}
SeriesFix <- function(lbnd, ubnd, by=1) {
  # Returns a NULL value when the lower boundary is greater than the upper
  mySeries <- c()
  if(lbnd < ubnd){
    mySeries <- seq(lbnd,ubnd, by=by)
  }
  return(mySeries)
}
SeriesFix2 <- function(lbnd, ubnd, by=1) {
  # Returns a NULL value when the lower boundary is greater than the upper
  mySeries <- c()
  if(lbnd <= ubnd){
    mySeries <- seq(lbnd,ubnd, by=by)
  }
  return(mySeries)
}
write_BESI_Scores <- function() {
  fileName = sprintf("%s/Datasets/BESI_scores.txt",targetDir)
  cat('',file = fileName,append = FALSE)
  info <- sprintf("%-50s %s","Sequence","Score")
  cat(info,file = fileName,sep = '',append = TRUE)
  cat('\n',file = fileName,append = TRUE)
  for (i in 1:length(BESI_Scores[,1]))
  {
    info <- sprintf("%-50s %s",BESI_Scores[i,1],BESI_Scores[i,2])
    cat(info,file = fileName,sep = '',append = TRUE)
    cat('\n',file = fileName,append = TRUE)
  }
}
WriteVLoopSeq <- function(seqName,HXB2,thisSeq,tmpSeq) { 
  # Determines and outputs V loops, length and score comparison data
  seqLength <- stringi::stri_length(HXB2[[1]])
  SeqResNums <- array(dim = seqLength)
  HXB2ResNums <- array(dim = seqLength)
  HXB2ResNum <- 1
  ResNum <- 1
  for (i in 1:seqLength) {
    if (stringi::stri_sub(HXB2,i,i) == "-") {
      SeqResNums[i] <- ResNum
      ResNum <- ResNum + 1
      HXB2ResNums[i] <- -1
    } else if (stringi::stri_sub(thisSeq,i,i) == "-") {
      HXB2ResNums[i] = HXB2ResNum
      HXB2ResNum <- HXB2ResNum + 1
      SeqResNums[i] <- -1
    } else {
      HXB2ResNums[i] = HXB2ResNum
      HXB2ResNum <- HXB2ResNum + 1
      SeqResNums[i] <- ResNum
      ResNum <- ResNum + 1
    }
  }
  for ( x in 1:length(VLOOPS)) {
    thisLoop <- list()
    tmp <- SeqResNums[which(HXB2ResNums ==VLOOPS[[x]][1]):which(HXB2ResNums ==VLOOPS[[x]][2])]
    sels <- tmp[which(tmp != -1)]
    loopLength <- length(sels)
    fileName <- sprintf("%s/VLoops/%s-V%i.seq",targetDir,seqName,x)
    fileLoopData <- sprintf("%s/Datasets/V%s-length_data.dat",targetDir,x)
    BESIScore <- BESI_Scores[which(BESI_Scores[,1] == seqName),2]
    loopData <- sprintf("%-60s %-10s %-15s\n",seqName,loopLength,BESIScore)
    cat(loopData,file = fileLoopData,append = TRUE)
    data <- sprintf(">%s\n",seqName)
    for(y in 1:loopLength) {
      thisLoop[y] <- stringi::stri_sub(tmpSeq,sels[y],sels[y])
    }
    cat(sprintf(">%s\n",seqName),file = fileName,append = FALSE)
    cat(as.character(thisLoop),sep = "",file = fileName,append = TRUE)
  }
}
