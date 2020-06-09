#!/usr/bin/env Rscript

#==============================================================================
#     BESI_EVM_analysis.r - A driver for libBESI_EVM.r
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

# BESI/EVM analysis tool 
#   author: Scott P. Morton
#   Center for Computational Science
#   Middle Tennessee State University
#   12/24/2018   

#### required external software ####
#
# MAFFT script is coded against version 7.273 from
# http://mafft.cbrc.jp/alignment/software/
#
# RAxml script is coded against version 8.2.11 from
# https://github.com/stamatak/standard-RAxML
#
# WebLogo 3.6.0 from
# https://github.com/WebLogo/weblogo
#
# ImageMagick binaries
#
# VMD
#
#### ####

# clear the environment first
remove(list = ls())

# set the working directory, where libBESI_EVM.r is located
setwd("~/Studies/Boeras")


# set the raxml binary to use
raxmlBin <- "raxmlHPC-PTHREADS-AVX"

# Capture command line args and flag env for cli options
# cliMode implies automation and assumes src and dest dirs are the same
args = commandArgs(trailingOnly=TRUE)
if (length(args) == 2){
  studyName <- args[1] # needs to be meaningful as this will be used in some graphics
  targetBase <- args[2]
  Dirs <- list(args[2])
} else {
  # Adjust these 3 parameters for manual mode operation
  studyName = 'Binding Energies' # needs to be meaningful as this will be used in some graphics
  # Creates a subfolder here of variable "studyName"
  targetBase <- "/home/scott/mounts/biosim1/Studies/Test"
  # Dirs is a list to allow manual runs from multiple source folders into a single destination folder
  Dirs <- list("/home/scott/mounts/biosim1/Studies/Test")
  
  # !!!!!!!!!!!!!!!!!!!!!!! Warning !!!!!!!!!!!!!!!!!!!!!!!
  #
  cat("Automatic mode requires exactly two (2) command line arguments:\n")
  cat("\t'The_Study_Name' and explicit src eg '/electrostatic/pipline/folder'\n")
  cat(sprintf("Current settings are:\n\t Study Name - %s\n\t Output Folder - %s\n\t Source Folder - %s\n",
      studyName,targetBase,Dirs))
  cat("Continue in manual mode y/N: ")
  if(askYesNo("Continue in manual mode?")){
    cat("continuing...\n")
  } else{quit(save = "ask")}
  #
  # !!!!!!!!!!!!!!!!!!!!!!! Warning !!!!!!!!!!!!!!!!!!!!!!!
}

doSens <- TRUE  # Tag phylo-tree with pH sensitivity mark and sort barplot by sensitivity
doBESI <- TRUE  # Perform Biomolecular Electro-Static Indexing, infer tree with gradient overlay
doBE <- TRUE    # Perform Binding Energies Analysis for complex models
uniformBEScale <- TRUE # Create BE graphs with a uniform scale across all graphs
doEVM <- TRUE   # Perform Electrostatic Variance Masking on residue data
doEVI <- FALSE   # (NOT IMPLEMENTED)Perform Electrostatic Variance Indexing on residue data (name subject to change)
doHXB2 <- TRUE  # Prepare HXB2 comparison data, plot residue loop data
doLoopAnalysis <- TRUE # Perform variable  loop comparison to BESI control

#### DO NOT EDIT ####
thisSrcDir <- getwd()
source(sprintf("%s/libBESI_EVM.r",thisSrcDir))
# findLineNum(sprintf("%s/libBESI_EVM.r#1044",thisSrcDir))
# setBreakpoint(sprintf("%s/libBESI_EVM.r#1044",thisSrcDir), envir = PerformBESIAnalysis)
#### DO NOT EDIT ###

#### These settings must come after the library is loaded ####
#### Init the library ###
### Titles, SubTitles and the Such

# Print Titles and Subtitles? TRUE / FALSE 
TnS <- FALSE

# Do inline sequence names in BESI Score Bar Graph? TRUE /FALSE
# Does not work well with large lists of sequences or sequences with long names
ISN <- FALSE

# Effects scales on certain graphs (fingerprints) Adjust to your liking, requires observation of resulting graphs
# This is meant to give a uniform graph scale across all EFP graphs
densityRanges <- matrix(c(-3.75,1.75,-2,10,-2,10),nrow=3,byrow = TRUE)

# Decide to print red horizontal line at 0 in binding energies graph TRUE/FALSE
BELineZero = FALSE
# Decide to print green trend line in binding energies graph TRUE/FALSE
BELineTrend = FALSE

# The variant used as the control in the study - trapped during processing, make sure this name matches the control name
CONTROL_NAME <- "Z242MPL25JAN03PCR23ENV1.1-_Donor_Transmitted"

# The control unbound ESSC data
template.U <- read.table(sprintf("%s/Datasets/UB_mean_template.data",targetDir))

# BESI Phylogeny plot x.lim factor, adjust to best fit long tip labels (greater than 'charLimit' characters)
charLimit <- 20
xLimFactor <- 0.035

# The columns to use for comparitive analysis of BE data (pH values) ie the Josh Sheet
myCols <- c("5.5","7.4")

# The subtitle for the variable loop vs BESI graphs
# Maybe this study does not have a control variant for this part, then set to " "
vloopSubTitle <- "Red indicates 'control'" 

# Weblogos settings
AR <- "--aspect-ratio 8"        # Ratio of stack height to width (default: 5)
SPL <- "--stacks-per-line 64"   # Maximum number of logo stacks per logo line.(default: 40)
                                # gets updated by the function 'generateVarianceStats()' to the number of
                                # residues selected by EVM
TFS <- "--title-fontsize 16"    # Title text font size in points (default: 12)
NFS <- "--number-fontsize 10"   # Axis numbers font size in points (default: 8)

# The last step is to prepare the EVM imagery
#   Copy from the vmdScripts repo folder 'VMD_gen_EVM_images.tcl' to the EVM folder 
#   You will need to check the orientation of the x,y,z axis along with the scale before executing
#   This assumes that all pdb files have been oriented to a single reference
#   This action will create a '.tga' and a '.pdf' file and requires the ImageMagick binaries to be installed on your machine

#### PROGRAM ####

# Verify directory structure exists
if(!(dir.exists(targetDir))){
  buildDirectory()
}

# prints EFP graphs, adjust 'densityRanges' above and rerun this section,
# this will overwrite EFP graphs and PCA data.
LoadChargeData_EFP_PCA()

if(doBESI){
  #### Basic required operations ####
  PrepareSequenceData()
  PerformBESIAnalysis()
}
 if(doBE){
   PerformBEAnalysis()
 }
 if(doEVM){
   PerformEVMAnalysis()
 }
 if(doEVI){
   PerformEVIAnalysis() # (NOT IMPLEMENTED) Future function (TBD)
 }
 if(doHXB2){
   AlignHXB2()
 }
if(doLoopAnalysis){
  PerformLoopAnalysis()
}