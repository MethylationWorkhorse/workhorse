
rm(list=ls(all=TRUE))

suppressWarnings(suppressPackageStartupMessages(require("tidyverse")) )
suppressWarnings(suppressPackageStartupMessages(require("stringr")) )
suppressWarnings(suppressPackageStartupMessages(require("glue")) )
suppressWarnings(suppressPackageStartupMessages(require("scales")) )
suppressWarnings(suppressPackageStartupMessages(require("matrixStats")) )

COM <- ","
TAB <- "\t"
RET <- "\n"

# Scratch space for workspace investigation

verbose <- 0

grn_idat <- '/Users/bbarnes/Documents/Projects/workhorse/idats/ReferenceBETA/201502830033/201502830033_R01C01_Grn.idat.gz'
red_idat <- '/Users/bbarnes/Documents/Projects/workhorse/idats/ReferenceBETA/201502830033/201502830033_R01C01_Red.idat.gz'
grn_dat  <- illuminaio::readIDAT(grn_idat)
red_dat  <- illuminaio::readIDAT(red_idat)



# End of file
