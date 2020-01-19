

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                             Brief Description::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# This script is designed to pick 10k probes for Dan's 48x1 chip format experiments.
#

rm(list=ls(all=TRUE))

suppressWarnings(suppressPackageStartupMessages(require("optparse",quietly=TRUE)))

suppressWarnings(suppressPackageStartupMessages(require("tidyverse")) )
suppressWarnings(suppressPackageStartupMessages(require("stringr")) )
suppressWarnings(suppressPackageStartupMessages(require("glue")) )
suppressWarnings(suppressPackageStartupMessages(require("scales")) )
suppressWarnings(suppressPackageStartupMessages(require("matrixStats")) )

# Parallel Computing Packages
suppressWarnings(suppressPackageStartupMessages(require("doParallel")) )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                             Default Params::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

doParallel::registerDoParallel()
num_cores   <- detectCores()
num_workers <- getDoParWorkers()

opt <- NULL
opt$prgmTag  <- 'selectMaxGroups'
opt$manDir   <- '/Users/bbarnes/Documents/Projects/manifests/methylation'
# opt$topDir   <- '/Users/bbarnes/Documents/Projects/workhorse'
opt$topDir   <- '/Users/bbarnes/Documents/Projects/workhorse/git/workhorse'
opt$sesDir   <- '/Users/bbarnes/Documents/Projects/darkmatter/Projects/sesamize'
opt$ssSrcDir <- file.path(opt$sesDir, 'sampleSheets/official')
opt$datDir   <- file.path(opt$topDir, 'dat')
opt$srcDir   <- file.path(opt$topDir, 'scripts/R')
opt$funSrc   <- file.path(opt$srcDir, 'probeSelection/selection_functions.R')
opt$wrkSrc   <- file.path(opt$srcDir, 'workhorse_functions.R')

opt$platform  <- 'EPIC-B4'
opt$build     <- 'hg19'

opt$passPercMinDetpNegs <- 97

opt$verbosity <- 3

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                           Command Line Options::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# Add command line later...


# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                             Preprocessing::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
stopifnot(file.exists(opt$funSrc))
stopifnot(file.exists(opt$wrkSrc))
source(opt$funSrc)
source(opt$wrkSrc)


# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#               Load Manifest Annotation to Select from::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

man.ann.tib <- getManifestAnnotation(verbose=opt$verbosity)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                  Load Probe Design Database to Select from::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# Full Manifest Design File::EPIC-B4???
man.des.csv <- file.path(opt$manDir, 'MethylationEPIC_v-1-0_B4.core.cpg-only.table.csv.gz')
man.des.tib <- suppressMessages(suppressWarnings(readr::read_csv(man.des.csv) )) %>%
  dplyr::arrange(IlmnID) %>% 
  dplyr::rename(Full_ID=IlmnID, Probe_ID=Name, U=AddressA_ID, M=AddressB_ID) %>%
  dplyr::mutate(U=str_remove(U, '^0+'), M=str_remove(M, '^0+'))


min.des.tib <- man.des.tib %>% dplyr::select(Probe_ID, Full_ID, CpgCnt, ScoreMin)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                  Load Genomic Coordinate Grouping Manifest::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

man.grp.rds <- file.path(opt$datDir, 'manifest', paste0(opt$platform,'.manifest.sesame.grouped.cpg-sorted.rds'))
man.grp.csv <- file.path(opt$datDir, 'manifest', paste0(opt$platform,'.manifest.sesame.grouped.cpg-sorted.csv.gz'))
if (file.exists(man.grp.rds)) {
  cat(glue::glue("[{opt$prgmTag}]: Loading Sesame Grouped Manifest({opt$platform}) RDS(CPG)={man.grp.rds}"),"\n", sep='')
  man.grp.tib <- readr::read_rds(man.grp.rds)
} else {
  stop("\n", glue::glue("[{opt$prgmTag}]: Grouping Manifest does not exist; RDS={man.grp.rds}"),"\n\n", sep='')
  q()
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                Build Human Currated Sample Sheet to Date::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

ss.curated.tib <- getCurratedHumanSampleSheet(opt=opt, verbose=opt$verbosity)

ss_curated_csv <- '/Users/bbarnes/Documents/Projects/workhorse/sampleSheets/EPIC-B4.SampleSheets.n4331.Dec1-2019.csv.gz'
ss_curated_tib <- readr::read_csv(ss_curated_csv)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                      Load all Auto-SampleSheets and 
#                         Filter Out Poor Samples::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# opt$ssAutoList <- list.files(file.path(opt$sesDir, 'idats_refs'), 
#                           pattern='EPIC-B0_both.SampleSheet.csv.gz', recursive=TRUE, full.names=TRUE)
# 
# opt$ssAutoList <- list.files(file.path(opt$sesDir, 'idats_refs'), 
#                           pattern='EPIC-B2.sset_core.Auto_SampleSheet.csv.gz', recursive=TRUE, full.names=TRUE)

opt$ssAutoList <- list.files(file.path(opt$topDir, 'builds/workhorse/References'), 
                             pattern='EPIC-B4_both.SampleSheet.csv.gz', recursive=TRUE, full.names=TRUE)

ss.auto.raw.tib <- suppressMessages(suppressWarnings(lapply(opt$ssAutoList, readr::read_csv))) %>%
  dplyr::bind_rows() %>%
  dplyr::mutate(Auto_Sample_DB_Key=stringr::str_remove(Auto_Sample_DB_Key, '_Beta$'),
                Auto_Sample_R2_Key=stringr::str_remove(Auto_Sample_R2_Key, '_Beta$'))

ss.auto.sam.tib <- ss.auto.raw.tib %>%
  dplyr::filter(PassDetpNegs_Percent_CG>=opt$passPercMinDetpNegs) %>%
  dplyr::arrange(Sentrix_Name) %>%
  split(.$Auto_Sample_R2_Key)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                      Load Tibs from Auto SampleSheets 
#                        & Filter Out Poor DetP Betas::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# man.ref.tib <- foreach (sample=names(ss.auto.sam.tib), .inorder=T, .final = function(x) setNames(x, names(ss.auto.sam.tib))) %dopar% {
sam_tib <- NULL
for (sample in names(ss.auto.sam.tib)) {
  # sample <- 'U'
  ss.auto.tib <- ss.auto.sam.tib[[sample]]
  ss.cur.tib  <- ss.curated.tib %>% dplyr::inner_join(ss.auto.tib, by="Sentrix_Name")
  # print(ss.cur.tib)
  
  ss.tib.locs <- file.path(opt$topDir, 'builds/workhorse/References', ss.cur.tib$Barcode, opt$platform)
  ss.tib.list <- paste0(ss.tib.locs,'/',ss.cur.tib$Sentrix_Name, '_EPIC-B4_both.Probes.tib.rds')
  # cat(ss.tib.locs,"\t",ss.tib.list, sep=' ')
  tibs <- lapply(ss.tib.list, loadBetaTib) %>%
    dplyr::bind_cols() %>% 
    dplyr::select(Probe_ID, starts_with('Beta'))
  
  sam_tib[[sample]] <- tibs
  # sam_tib <- sam_tib %>% dplyr::bind_rows(.id=sample)
}
# Quick bind_cols verification on via size::
#   lapply(sam_tib, base::nrow) %>% unlist() %>% unique


# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                              Main Workflow::
#
#   1. Proof of Concept: Show Group Variance for Largest Groups
#      - Filter/Group man.grp.tib
#      - Join man.grp.tib with each Beta
#      - Plot Group Stats by: size vs. sample type
#
#   2. Reselect and Expand EPIC designs
#      - Create Search Space Buckets to tally probes
#      - Traversing the largest groups selecting probes based on:
#        A. Is EPIC::
#           a. Is Bin Available::
#           b. Select as many same locus sites as allowed
#        B. NOT EPIC::
#           a. Is Poor Performing Probe and Bin Available::
#
#   Pre-Data::
#
#   1. man.des.tib:: Design Manifest (Scores, #CpGs, Strands, ...)
#      min.des.tib:: Subset of data
#   2. man.grp.tib:: Groupd Manifest (Grouping Data and Genomic Coordinates)
#   3. sam_tib[s] :: Beta Values already filtered by DetP > [0.05/0.02]
#   4. ss.auto.sam.tib:: Sample Sheet Annotation (Human/Auto)
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

#   1. Proof of Concept: Show Group Variance for Largest Groups
#      - Filter/Group man.grp.tib
man.top.tib <- man.grp.tib %>% dplyr::ungroup() %>% 
  dplyr::filter(Local_Group_Cnt>4) %>% 
  dplyr::group_by(seqnames, Local_Group_Key, Local_Group_Cnt)

# Top Groups Sorted by Probe_ID
cpg.top.tib <- man.top.tib %>% dplyr::arrange(Probe_ID)

# On-target Summary::
sum.25k.tib <- cpg.top.tib %>% dplyr::left_join(min.des.tib, by="Probe_ID") %>% 
  dplyr::group_by(CpgCnt, as.integer(ScoreMin*10)) %>% 
  dplyr::summarise(Count=n()) %>% 
  purrr::set_names(c('CpgCnt', 'Score', 'Count')) %>% 
  dplyr::filter(!is.na(CpgCnt), !is.na(Score)) %>%
  dplyr::ungroup() %>% 
  dplyr::mutate(CpgCnt=as.integer(CpgCnt))

# gg <- ggplot(data=sum.25k.tib, aes(x=CpgCnt, y=Score, z=Count, fill=Sample)) +

library(hexbin)
library(plotly)
library(gridExtra)
library(grid)
library(akima)
library(rgl)
library(latticeExtra)
library(png)

while (dev.cur()>1) dev.off()

cc <- cloud(sum.25k.tib %>% as.matrix, panel.3d.cloud=panel.3dbars, col.facet='grey', 
            xbase=0.4, ybase=0.4, scales=list(arrows=FALSE, col=1), 
            par.settings = list(axis.line = list(col = "transparent")))



gg <- ggplot(data=sum.25k.tib, aes(x=CpgCnt, y=Score, fill=Count)) +
  ggplot2::geom_density_2d(alpha=0.5)

cc <- cloud(y~x+z, sum.25k.tib, panel.3d.cloud=panel.3dbars, col.facet='grey', 
      xbase=0.4, ybase=0.4, scales=list(arrows=FALSE, col=1), 
      par.settings = list(axis.line = list(col = "transparent")))


grp.stats <- NULL
for (sample in names(sam_tib)) {
  #      - Join man.grp.tib with each Beta
  # sample <- 'U'
  
  grp.stats[[sample]] <- cpg.top.tib %>% dplyr::inner_join(sam_tib[[sample]], by="Probe_ID") %>%
    dplyr::arrange(seqnames, Local_Group_Key, Local_Group_Cnt) %>%
    dplyr::group_by(seqnames, Local_Group_Key, Local_Group_Cnt) %>%
    dplyr::summarise(Beta_Avg=mean(Beta, na.rm=TRUE), 
                     Beta_Med=median(Beta, na.rm=TRUE), 
                     Beta_SD=sd(Beta, na.rm=TRUE))

  # gg <- ggplot(data=grp.stats) +
  #   ggplot2::geom_density(aes(x=Beta_Avg), fill="Red") +
  #   ggplot2::geom_density(aes(x=Beta_Med), fill="Blue")
}
sam.stats <- grp.stats %>% dplyr::bind_rows(.id='Sample')

gg <- ggplot(data=sam.stats, aes(x=Beta_Avg, fill=Sample)) +
  ggplot2::geom_density(alpha=0.5)

# man.des.tib


ctime <- Sys.time()
cat(glue::glue("[{opt$prgmTag}]: Finished(time={ctime})"),"\n\n",sep='')

# End of file
