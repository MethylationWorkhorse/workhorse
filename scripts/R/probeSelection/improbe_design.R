
rm(list=ls(all=TRUE))

suppressWarnings(suppressPackageStartupMessages(require("optparse",quietly=TRUE)))

suppressWarnings(suppressPackageStartupMessages(require("tidyverse")) )
suppressWarnings(suppressPackageStartupMessages(require("stringr")) )
suppressWarnings(suppressPackageStartupMessages(require("glue")) )
suppressWarnings(suppressPackageStartupMessages(require("scales")) )
suppressWarnings(suppressPackageStartupMessages(require("matrixStats")) )
suppressWarnings(suppressPackageStartupMessages(require("grid")) )

suppressWarnings(suppressPackageStartupMessages(require("Biostrings")) )

# Parallel Computing Packages
suppressWarnings(suppressPackageStartupMessages(require("doParallel")) )

doParallel::registerDoParallel()
num_cores   <- detectCores()
num_workers <- getDoParWorkers()

COM <- ","
TAB <- "\t"
RET <- "\n"
BNG <- "|"

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                             Default Params::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

opt <- NULL
opt$prgmTag <- 'improbe_design'

opt$single   <- FALSE
opt$parallel <- FALSE

opt$desDir <- NULL
opt$manDir <- NULL
opt$outDir <- NULL

opt$name   <- NULL

opt$prbSrc <- NULL
opt$format <- NULL
opt$subCnt <- NULL
opt$blockSize <- 10000

opt$exactMatch <- FALSE
opt$impExecute <- FALSE

opt$verbosity <- 3

improbeHeader <- c('Seq_ID', 'Sequence', 'Genome_Build', 'Chromosome', 'Coordinate', 'CpG_Island')

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                       System (Linux/Mac) Directories::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
opt$mac_dir <- '/Users/bbarnes/Documents/Projects'
opt$lix_dir <- '/illumina/scratch/darkmatter/data'
opt$sys_key <- NULL

if (dir.exists(opt$mac_dir)) {
  opt$sysName <- 'mac'
  opt$prgmTag  <- paste(opt$prgmTag, 'RStudio', sep='_')
  opt$prgDir   <- opt$mac_dir
  opt$manDir   <- '/Users/bbarnes/Documents/Projects/manifests/methylation'
  opt$desDir   <- '/Users/bbarnes/Documents/Projects/methylation/NZT_Limitless/data/imDesignOutput'
  opt$topDir   <- file.path(opt$mac_dir, 'workhorse')
  opt$gitDir   <- file.path(opt$topDir, 'git/workhorse')

  # NOT REAL ON MAC::
  opt$improbe  <- file.path(opt$lix_dir, 'bin/improbe')
  opt$imp13mer <- file.path(opt$lix_dir, 'dat/improbe/human-36.1-methyl-and-unmethyl-13mer-s3-for-infinium-methylation.dat')
  opt$impTango <- file.path(opt$lix_dir, 'dat/improbe/Tango_A_or_B_11mer_s1.dat')
} else if (dir.exists(opt$lix_dir)) {
  opt$sysName <- 'linux'
  opt$prgDir   <- opt$lix_dir
  opt$manDir   <- '/illumina/scratch/darkmatter/manifests/data/methylation'
  opt$desDir   <- NULL
  opt$topDir   <- file.path(opt$lix_dir)
  opt$gitDir   <- file.path(opt$topDir, 'git/workhorse')
  
  opt$improbe  <- file.path(opt$lix_dir, 'bin/improbe')
  opt$imp13mer <- file.path(opt$lix_dir, 'dat/improbe/human-36.1-methyl-and-unmethyl-13mer-s3-for-infinium-methylation.dat')
  opt$impTango <- file.path(opt$lix_dir, 'dat/improbe/Tango_A_or_B_11mer_s1.dat')
} else {
  stop(glue::glue("[ERROR]: Unrecognized prgDir=({opt$mac_dir}, {opt$lix_dir})"),"\n", sep='')
  q()
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                           Command Line Options::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
args.dat <- commandArgs(trailingOnly = FALSE)
opt_parser <- NULL
if (args.dat[1]=='RStudio') {
  opt$srcDir   <- file.path(opt$gitDir, 'scripts')
  opt$prgmPath <- file.path(opt$srcDir, 'R/plotsByExperiment.R')
  
  opt$format <- 'des'
  opt$format <- 'man'
  
  opt$prbSrc <- 'cg'
  opt$prbSrc <- 'rs'
  opt$prbSrc <- 'ch'
  
  opt$subCnt <- 20000
  
  opt$name   <- 'EPIC-B2'
  
  opt$exactMatch <- FALSE
  opt$impExecute <- FALSE
  
  opt$outDir <- file.path(opt$topDir, 'workspace/manifest')
  
} else {
  opt$prgmPath <- substring(args.dat[grep("--file=", args.dat)], 8)
  opt$runMode  <- 'CommandLine'
  opt$prgmTag  <- sub('\\.R$', '', basename(substring(args.dat[grep("--file=", args.dat)], 8)))
  # opt$srcDir   <- dirname(normalizePath(dirname(substring(args.dat[grep("--file=", args.dat)], 8)) ))
  opt$srcDir   <- dirname(dirname(normalizePath(dirname(substring(args.dat[grep("--file=", args.dat)], 8)) )))
  opt$gitDir   <- file.path(opt$srcDir, '../../..')
  
  args.dat <- commandArgs(trailingOnly = TRUE)
  option_list = list(
    make_option(c("--Rscript"), type="character", default=opt$Rscript, 
                help="Rscript path [default= %default]", metavar="character"),
    make_option(c("--prgmTag"), type="character", default=opt$prgmTag, 
                help="Program name [default= %default]", metavar="character"),
    make_option(c("--sysName"), type="character", default=opt$sysName, 
                help="System Name. Should be auto-detected [linux, mac] [default= %default]", metavar="character"),
    
    # Directories/Files::
    make_option(c("--topDir"), type="character", default=opt$topDir, 
                help="Top data directory location [default= %default]", metavar="character"),
    make_option(c("--srcDir"), type="character", default=opt$srcDir, 
                help="Source data directory location [default= %default]", metavar="character"),
    make_option(c("--manDir"), type="character", default=opt$manDir, 
                help="Manifest directory location [default= %default]", metavar="character"),
    make_option(c("--desDir"), type="character", default=opt$desDir, 
                help="Improbe Design data directory location [default= %default]", metavar="character"),
    make_option(c("--outDir"), type="character", default=opt$outDir, 
                help="Output directory location [default= %default]", metavar="character"),
    make_option(c("--name"), type="character", default=opt$name, 
                help="Output name [default= %default]", metavar="character"),
    
    # Processing Parameters::
    make_option(c("--single"), action="store_true", default=opt$single, 
                help="Boolean variable to run a single sample on a single-core [default= %default]", metavar="boolean"),
    make_option(c("--parallel"), action="store_true", default=opt$parallel, 
                help="Boolean variable to run parallel on multi-core [default= %default]", metavar="boolean"),
    
    # Run Parameters::
    make_option(c("--prbSrc"), type="character", default=opt$prbSrc,
                help="Probe Type (cg,ch,rs) [default= %default]", metavar="character"),
    make_option(c("--format"), type="character", default=opt$format,
                help="Input file format  [default= %default]", metavar="character"),
    make_option(c("--subCnt"), type="double", default=opt$subCnt, 
                help="Sub Sample Input Count [default= %default]", metavar="double"),
    make_option(c("--blockSize"), type="double", default=opt$blockSize, 
                help="Loci (x4) block size (mostly for cluster and memory partitioning) [default= %default]", metavar="double"),
    make_option(c("--exactMatch"), action="store_true", default=opt$exactMatch, 
                help="Boolean variable for exact probe matching vs. mismatch counting [default= %default]", metavar="boolean"),
    
    # improbe Parameters::
    make_option(c("--impExecute"), action="store_true", default=opt$impExecute, 
                help="Boolean variable for executing improbe [default= %default]", metavar="boolean"),
    
    make_option(c("--improbe"), type="character", default=opt$improbe,
                help="improbe executable [default= %default]", metavar="character"),
    make_option(c("--impTango"), type="character", default=opt$impTango,
                help="improbe tangos file [default= %default]", metavar="character"),
    make_option(c("--imp13mer"), type="character", default=opt$imp13mer,
                help="improbe 13mer file [default= %default]", metavar="character"),
    
    # Verbosity::
    make_option(c("-v", "--verbosity"), type="double", default=opt$verbosity, 
                help="0-5 (5 is very verbosity) [default= %default]", metavar="double")
  )
  opt_parser = OptionParser(option_list=option_list)
  opt = parse_args(opt_parser)
}

if (is.null(opt$prgmTag)  || 
    is.null(opt$prbSrc)   ||
    is.null(opt$format)   ||
    is.null(opt$topDir)   ||
    is.null(opt$manDir)   ||
    is.null(opt$outDir)   ||
    is.null(opt$name)     ||
    is.null(opt$srcDir) ) {
  
  if (!is.null(opt_parser)) print_help(opt_parser)
  print(as.data.frame(tidyr::gather(dplyr::bind_rows(opt),"Options", "Value")))
  stop("Null arguments!\n\n")
  q()
}
if (is.null(opt$format) && is.null(opt$desDir)) {
  if (!is.null(opt_parser)) print_help(opt_parser)
  print(as.data.frame(tidyr::gather(dplyr::bind_rows(opt),"Options", "Value")))
  stop("If using improbe design files for analysis you must provide improbe design directory (--desDir)!\n\n")
  q()
}

# Automatic Params
opt$MinMisMatch <- 0
if (opt$prbSrc=='ch') opt$MinMisMatch <- 1

# R Source Code Files
opt$fIOSrc   <- file.path(opt$srcDir, 'R/workhorse_IO_functions.R')
opt$wrkSrc   <- file.path(opt$srcDir, 'R/workhorse_functions.R')
opt$plotSrc  <- file.path(opt$srcDir, 'R/workhorse_plotting_functions.R')
opt$impSrc   <- file.path(opt$srcDir, 'R/probeSelection/improbe_functions.R')

if (opt$verbosity>0) print(as.data.frame(tidyr::gather(dplyr::bind_rows(opt),"Options", "Value")))

stopifnot(file.exists(opt$fIOSrc))
stopifnot(file.exists(opt$plotSrc))
stopifnot(file.exists(opt$impSrc))
# stopifnot(file.exists(opt$wrkSrc))

source(opt$fIOSrc)
source(opt$plotSrc)
source(opt$impSrc)
# source(opt$wrkSrc)

# Output Directory::
if (!dir.exists(opt$outDir)) dir.create(opt$outDir, recursive=TRUE)
opt$out_man_csv <- file.path(opt$outDir, paste(opt$name,'_',opt$prbSrc,'.manifest.csv.gz', sep=''))

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                             Default Files::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# improbe Design Files::
cpg_des_tsv <- file.path(opt$desDir, 'hs_ref_both_both_RefSeq37_probeDesignOutput.tsv.gz')
cpg_des_tsv <- file.path(opt$desDir, 'hs_ref_both_both_RefSeq37_probeDesignOutput.1k.tsv.gz')
cph_des_tsv <- file.path(opt$desDir, 'ch.improbe-designOutput.tsv.gz')

# 450k Formats::
# cph_man_csv <- file.path(opt$manDir, 'ch-repair/HumanMethylation450_15017482_v.1.2.ch-noAnnotation.csv.gz')
# snp_man_csv <- file.path(opt$manDir, 'rs-repair/HumanMethylation450_15017482_v.1.2.rs.no-FwdSeq.csv.gz')
# snp_swp_tsv <- file.path(opt$manDir, 'rs-repair/hs_ref_RefSeq37_probeDesign_SWAP.tsv.gz')

# EPIC Formats::
# cpg_man_csv <- file.path(opt$manDir, 'MethylationEPIC_v-1-0_B4.core.cpg-only.table.1k.csv.gz')
# cpg_man_csv <- file.path(opt$manDir, 'MethylationEPIC_v-1-0_B4.core.cpg-only.table.csv.gz')
cpg_man_csv <- file.path(opt$manDir, 'MethylationEPIC_v-1-0_B2.csv.gz')
cpg_man_csv <- file.path(opt$manDir, 'MethylationEPIC_v-1-0_B2.1k.csv.gz')

cph_man_csv <- file.path(opt$manDir, 'ch-repair/MethylationEPIC_v-1-0_B4.ch-noAnnotation.csv.gz')
snp_man_csv <- file.path(opt$manDir, 'rs-repair/MethylationEPIC_v-1-0_B4.rs.no-header.csv.gz')
snp_swp_tsv <- file.path(opt$manDir, 'rs-repair/rs_swap_improbe_input.tsv.gz')

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                         Load/Properly-Build SNP Manifest::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# if (opt$format=='des' && opt$prbSrc=='cg') src_man_tib <- loadImprobeDesign(cpg_des_tsv)
# if (opt$format=='man' && opt$prbSrc=='cg') cpg_man_tib <- loadManifestCG(cpg_man_csv)
# if (opt$format=='man' && opt$prbSrc=='ch') cph_man_tib <- loadManifestCH(cph_man_csv, ry='R')
# if (opt$format=='man' && opt$prbSrc=='rs') snp_man_tib <- loadManifestRS(snp_man_csv, snp_swp_tsv)

# MAN - Manifest
# DES - Template Sequence
# PRB - Probes
# IMP - Improbe

cpg_man_tib <- loadManifestCG(cpg_man_csv) %>% dplyr::mutate(PRB_DES='cg', CHR=as.character(CHR))
cph_man_tib <- loadManifestCH(cph_man_csv, ry='R') %>% dplyr::mutate(PRB_DES='ch', CHR=as.character(CHR))
snp_man_tib <- loadManifestRS(snp_man_csv, snp_swp_tsv) %>% dplyr::mutate(PRB_DES='rs', CHR=as.character(CHR))
src_man_tib <- dplyr::bind_rows(cpg_man_tib, cph_man_tib, snp_man_tib)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                                 Sub Sample::
if (!is.null(opt$subCnt) && opt$subCnt!=0)
  src_man_tib <- src_man_tib %>% head(n=opt$subCnt)
src_man_len <- src_man_tib %>% base::nrow()

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                            Format Manifest::
src_man_tib <- src_man_tib %>% select(IlmnID,Forward_Sequence,Genome_Build,CHR,MAPINFO,PRB_DES,diNUC) %>%
  dplyr::rename(Seq_ID=IlmnID,Chromosome=CHR,Coordinate=MAPINFO,IUPAC_Forward_Sequence=Forward_Sequence) %>%
  dplyr::mutate(Sequence=stringr::str_replace(IUPAC_Forward_Sequence, '\\[[A-Za-z][A-Za-z]\\]', '[CG]'),
                CpG_Island='FALSE',PRB_CGN=Seq_ID) %>%
  dplyr::select(improbeHeader, PRB_CGN,PRB_DES,diNUC,IUPAC_Forward_Sequence)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                    Include Full Design Data from improbe::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
imp_inp_tib <- src_man_tib %>% dplyr::select(improbeHeader)
imp_out_key <- paste(opt$name,opt$prbSrc,paste0(1,'-',src_man_len), sep='_')
imp_des_tsv <- writeImprobeInput(tib=imp_inp_tib,name=imp_out_key,dir=opt$outDir,run=opt$impExecute,
                                 exe=opt$improbe,impTango=opt$impTango,imp13mer=opt$imp13mer,
                                 verbose=opt$verbosity)

if (opt$sysName=='mac') {
  imp_cpg_tsv <- file.path(opt$topDir, 'improbe/fromCluster/EPIC-B2_cg_1-1000.improbe-output.tsv.gz')
  imp_cph_tsv <- file.path(opt$topDir, 'improbe/fromCluster/EPIC-B2_ch_1-2932.improbe-output.tsv.gz')
  imp_snp_tsv <- file.path(opt$topDir, 'improbe/fromCluster/EPIC-B2_rs_1-58.improbe-output.tsv.gz')
  
  imp_cpg_tib <- loadImprobeDesign(imp_cpg_tsv) %>% dplyr::mutate(Chromosome=as.character(Chromosome))
  imp_cph_tib <- loadImprobeDesign(imp_cph_tsv) %>% dplyr::mutate(Chromosome=as.character(Chromosome))
  imp_snp_tib <- loadImprobeDesign(imp_snp_tsv) %>% dplyr::mutate(Chromosome=as.character(Chromosome))
  imp_des_tib <- dplyr::bind_rows(imp_cpg_tib,imp_cph_tib,imp_snp_tib) # %>%
  #  dplyr::select(-c(Genome_Build,Chromosome,Coordinate))
  
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                      Calcluate Probes on All Strands::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
des_seq_F_C <- src_man_tib %>% 
  dplyr::mutate(FR=TRUE, CO=TRUE, DesSeqN=shearBrac(IUPAC_Forward_Sequence))

des_seq_R_C <- des_seq_F_C %>% dplyr::mutate(
  FR=!FR,CO=CO, DesSeqN=revCmp(DesSeqN) )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                        Build All Design Strands::
bsc_tibs <- NULL

# BSC-Forward-Converted::
bsc_tibs$FC <- des_seq_F_C %>% dplyr::mutate(
  DesBscU = bscUs(DesSeqN),
  DesBscM = bscMs(DesSeqN),
  DesBscD = bscDs(DesSeqN) ) %>%
  dplyr::mutate(Seq_ID)
# BSC-Foward-Opposite::
bsc_tibs$FO <- bsc_tibs$FC %>% dplyr::mutate(
  FR=FR,CO=!CO, 
  # DesSeqN=revCmp(DesSeqN),
  DesBscU=revCmp(DesBscU),
  DesBscM=revCmp(DesBscM),
  DesBscD=revCmp(DesBscD) ) %>%
  dplyr::mutate(Seq_ID)

# BSC-Reverse-Converted::
bsc_tibs$RC <- des_seq_R_C %>% dplyr::mutate(
  DesBscU = bscUs(DesSeqN),
  DesBscM = bscMs(DesSeqN),
  DesBscD = bscDs(DesSeqN) ) %>%
  dplyr::mutate(Seq_ID)
# BSC-Reverse-Opposite::
bsc_tibs$RO <- bsc_tibs$RC %>% dplyr::mutate(
  FR=FR,CO=!CO, 
  # DesSeqN=revCmp(DesSeqN),
  DesBscU=revCmp(DesBscU),
  DesBscM=revCmp(DesBscM),
  DesBscD=revCmp(DesBscD) ) %>%
  dplyr::mutate(Seq_ID)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                  Compare Manifest with Calculated Probes::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
prb_tibs <- NULL
prb_tibs <- foreach (srd=names(bsc_tibs), .combine=rbind) %dopar% {
  lapply(split(bsc_tibs[[srd]], bsc_tibs[[srd]]$PRB_DES), desAllPrbs) %>% dplyr::bind_rows()
}

prb_des_tibs <- prb_tibs %>% 
  dplyr::inner_join(imp_des_tib, 
                    by=c('Seq_ID','FR','CO','Genome_Build','Chromosome','Coordinate'),
                    suffix=c('_IUP','_IMP'))

# Now we have comparisons::
prb_mat_tibs <- prb_des_tibs %>% 
  cmpInfI_MisMatch('PRB1_U_IUP', 'PRB1_U_IMP', 'PRB1_M_IUP', 'PRB1_M_IMP',
                   verbose=opt$verbosity) # %>% dplyr::select(Man_MisMatch,Man_TarMatch)

plot_ord_tib <- prb_mat_tibs %>% dplyr::arrange(-Man_MisMatch) %>% 
  dplyr::group_by(FR,TB,CO,PRB_DES,Man_MisMatch,Man_TarMatch) %>% 
  dplyr::top_n(1,wt=Seq_ID) %>% 
  dplyr::select(Seq_ID,FR,TB,CO,PRB_DES,Man_MisMatch,Man_TarMatch)

plot_nrows <- plot_ord_tib %>% base::nrow()

for (ii in seq(1,plot_nrows)) {
  tag_tib <- plot_ord_tib[ii,]
  print(tag_tib)
  
  cur_tib <- prb_mat_tibs %>% dplyr::filter(Seq_ID==tag_tib$Seq_ID[1])
  tibs <- cur_tib %>% srdsToBrac()
  
  # Build Printable Strings for reach strand::
  cat(glue::glue("IlmnID={IlmnId}:{RET}"))
  strs <- NULL
  strs$F <- prbsToStr(tibs$F, pr=tag_tib$PRB_DES[1], verbose=opt$verbosity)
  strs$R <- prbsToStr(tibs$R, pr=tag_tib$PRB_DES[1], verbose=opt$verbosity)
  
  # cur_tib <- mat_tib %>% dplyr::filter(IlmnID==IlmnId)
  if (FALSE && cur_tib$Infinium_Design_Type[1]=='II') {
    cur_tib %>% dplyr::select(FR,CO, PRB2_D_IUP, PRB2_D_IMP,Man_MisMatch,Man_TarMatch) %>% print()
  } else {
    cur_tib %>% dplyr::select(FR,CO, PRB1_U_IUP, PRB1_U_IMP,Man_MisMatch,Man_TarMatch) %>% print()
    cur_tib %>% dplyr::select(FR,CO, PRB1_M_IUP, PRB1_M_IMP,Man_MisMatch,Man_TarMatch) %>% print()
  }
  
  break
}



# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                          END OF CURRENT WORK::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# TBD::
#
#  1. RS/CH A/B slot dual improbe designs for scores
#  2. Testing: mm10, Horvath, NZT, 1um-15k
#
#  3. Plot CSS vs. SCR_U/SCR_M and Variability
#     - CSS-Area
#
#         CSS_U  vs. SCR_U
#        _____
#       /    /
#      /    /   CSS_M  vs. SCR_M
#     /    /
#    -----
#

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                          END OF CURRENT WORK::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #





prb_tibs <- NULL
prb_tibs <- foreach (srd=names(bsc_tibs), .combine=rbind) %dopar% {
  bsc_pr_tibs <- bsc_tibs[[srd]] %>% split(.$PRB_DES)
  lapply(bsc_pr_tibs, desAllPrbs) %>% dplyr::bind_rows()
}


prb_tibs <- NULL
# for (srd in names(bsc_tibs)) { # break }
prb_tibs <- foreach (srd=names(bsc_tibs), .combine=rbind) %dopar% {
  fr <- bsc_tibs[[srd]] %>% dplyr::distinct(FR) %>% base::as.logical()
  co <- bsc_tibs[[srd]] %>% dplyr::distinct(CO) %>% base::as.logical()
  bsc_pr_tibs <- bsc_tibs[[srd]] %>% split(.$PRB_DES)
  lapply(bsc_pr_tibs, desAllPrbs, fr,co,pr) %>% dplyr::bind_rows()
}


  # prb_tib <- NULL
  # for (pr in names(bsc_pr_tibs)) {
  #   prb_tib <- prb_tib %>% 
  #     dplyr::bind_rows(bsc_pr_tibs[[pr]] %>% 
  #                        des2prbs(fwd=fr, con=co, pr=pr, mu='U', strand='DesBscU') %>%
  #                        des2prbs(fwd=fr, con=co, pr=pr, mu='M', strand='DesBscM') %>%
  #                        des2prbs(fwd=fr, con=co, pr=pr, mu='D', strand='DesBscD') )
  # }
  # prb_tib
}


    
  prb_tib <- prb_tib %>% dplyr::bind_rows(prb)
  
  # Infinium II Match::
  mat2_nrows <- 0
  if (!is.null(prb_man_tib[['II']])) {
    if (opt$exactMatch) {
      mat2 <- prb %>% dplyr::inner_join(prb_man_tib[['II']], by="IlmnID") %>%
        cmpInfII('PRB2_D', 'AlleleA_ProbeSeq', mu='D', verbose=opt$verbosity)
    } else {
      mat2 <- prb %>% dplyr::inner_join(prb_man_tib[['II']], by="IlmnID") %>%
        cmpInfII_MisMatch('PRB2_D', 'AlleleA_ProbeSeq', mu='D', verbose=opt$verbosity)
    }
    
    mat2_nrows <- mat2 %>% dplyr::filter(Man_TarMatch & Man_MisMatch < 3) %>%
      dplyr::distinct(IlmnID) %>% base::nrow()
    mat_tib <- mat_tib %>% dplyr::bind_rows(mat2)
  }
  
  # Infinium I Match::
  mat1_nrows <- 0
  if (!is.null(prb_man_tib[['I']])) {
    if (opt$exactMatch) {
      mat1 <- prb %>% dplyr::inner_join(prb_man_tib[['I']], by="IlmnID") %>% 
        cmpInfI('PRB1_U', 'AlleleA_ProbeSeq',
                'PRB1_M', 'AlleleB_ProbeSeq', verbose=opt$verbosity)
    } else {
      mat1 <- prb %>% dplyr::inner_join(prb_man_tib[['I']], by="IlmnID") %>% 
        cmpInfI_MisMatch('PRB1_U', 'AlleleA_ProbeSeq',
                         'PRB1_M', 'AlleleB_ProbeSeq', verbose=opt$verbosity)
    }

    mat1_nrows <- mat1 %>% dplyr::filter(Man_TarMatch & Man_MisMatch < 3) %>% 
      dplyr::distinct(IlmnID) %>% base::nrow()
    mat_tib <- mat_tib %>% dplyr::bind_rows(mat1)
  }
  
  mat_nrows <- mat_tib %>% dplyr::distinct(IlmnID) %>% base::nrow()
  cat(glue::glue("[Matching]: fr={fr}, co={co}, I={mat1_nrows}, II={mat2_nrows}, Total={mat_nrows}.{RET}"))
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                          Visualization of Probes::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

vizProbes <- FALSE
vizProbes <- TRUE

if (vizProbes) {
  prb_tibs <- prb_tib %>% split(.$IlmnID)
  mat_tibs <- mat_tib %>% dplyr::filter(Man_TarMatch & Man_MisMatch<=opt$MinMisMatch) %>% split(.$IlmnID)

  IlmnIds <- names(prb_tibs)
  
  for (IlmnId in IlmnIds) { # break }
    # IlmnId <- 'rs10033147'
    # IlmnId <- 'ch.1.101940785F'
    # IlmnId <- 'ch.1.101940785F'
    #  R_O should be shift to one to the right...
    # IlmnId=rs715359 Current RS Check
    #
    # Hovarth CH::
    #  IlmnId <- 'ch.2.30415474F'
    #  IlmnId <- 'ch.3.741789F'
    
    if (length(mat_tibs[[IlmnId]])==0) {
      # Add Bracs & Split On FR Plot Con/Opp::
      tibs <- prb_tibs[[IlmnId]] %>% srdsToBrac()
      
      # Build Printable Strings for reach strand::
      cat(glue::glue("IlmnID={IlmnId}:{RET}"))
      strs <- NULL
      strs$F <- prbsToStr(tibs$F, pr=opt$prbSrc, verbose=opt$verbosity)
      strs$R <- prbsToStr(tibs$R, pr=opt$prbSrc, verbose=opt$verbosity)

      cur_tib <- mat_tib %>% dplyr::filter(IlmnID==IlmnId)
      if (cur_tib$Infinium_Design_Type[1]=='II') {
        cur_tib %>% dplyr::select(FR,CO, PRB2_D, AlleleA_ProbeSeq,Man_MisMatch,Man_TarMatch) %>% print()
      } else {
        cur_tib %>% dplyr::select(FR,CO, PRB1_U, AlleleA_ProbeSeq,Man_MisMatch,Man_TarMatch) %>% print()
        cur_tib %>% dplyr::select(FR,CO, PRB1_M, AlleleB_ProbeSeq,Man_MisMatch,Man_TarMatch) %>% print()
      }
      
      # break
    } else {

      # break
    }
  }
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                 Select Top Ranked Match to Manifest::
#
#  1. Remove PRB1_U==PRB1_M
#  2. Rank Man_MisMatch, Man_TarMatch
mat_tib <- mat_tib %>% dplyr::filter(PRB1_U!=PRB1_M) %>% 
  dplyr::arrange(IlmnID,Man_MisMatch,!Man_TarMatch) %>% 
  dplyr::group_by(IlmnID)  %>% top_n(1, wt=-Man_MisMatch) %>% 
  dplyr::distinct(IlmnID, .keep_all=TRUE) %>% 
  dplyr::ungroup()
  # dplyr::select(IlmnID,FR,CO,Man_MisMatch,Man_TarMatch)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                          Build New Manifest::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
# Calculate New Fields::
#  - FLAG Matched and Failed
#  - Calculate T/B from Forward_Sequence
#    mat_tib$Forward_Sequence %>% str_sub(1,60) ...
#
mat_man_tib <- mat_tib %>% 
  dplyr::mutate(
    Strand_FR=dplyr::case_when( FR ~ 'F', !FR ~ 'R', TRUE ~ NA_character_ ),
    Strand_TB='U',
    Strand_CO=dplyr::case_when( CO ~ 'C', !CO ~ 'O', TRUE ~ NA_character_ ),
    Infinium=case_when( is.na(AlleleB_ProbeSeq) ~ 'II', is.na(AlleleA_ProbeSeq) ~ NA_character_, TRUE ~ 'I'),
    Next_Base_Ref=NXB_D,
    Color_Channel_Ref=case_when(Next_Base_Ref=='A'|Next_Base_Ref=='T'|Next_Base_Ref=='W'|
                                  Next_Base_Ref=='a'|Next_Base_Ref=='t'|Next_Base_Ref=='w' ~ 'Red',
                                Next_Base_Ref=='C'|Next_Base_Ref=='G'|Next_Base_Ref=='S'|
                                  Next_Base_Ref=='c'|Next_Base_Ref=='g'|Next_Base_Ref=='s' ~ 'Grn',
                                TRUE ~ 'Mix'),
    Name=paste(CGN,Strand_FR,Strand_TB,Strand_CO,Infinium, sep='_'),
    AddressA_ID=NA,
    AddressB_ID=NA
  ) %>% 
  dplyr::select(IlmnID,Name,AddressA_ID,AlleleA_ProbeSeq,AddressB_ID,AlleleB_ProbeSeq,
                CGN,Strand_FR,Strand_TB,Strand_CO,Infinium,Next_Base_Ref,Color_Channel_Ref,
                Forward_Sequence,DesSeqN,DesBscU,DesBscM,DesBscD,
                PRB1_U,PRB1_M,PRB2_D,
                NXB_U,NXB_M,NXB_D,
                CPN_U,CPN_M,CPN_D,
                TAR_U,TAR_M,TAR_D,
                BOD_U,BOD_M,BOD_D,
                END_U,END_M,END_D,
                Man_MisMatch,Man_TarMatch)

new_man_tib <- src_man_tib %>% dplyr::full_join(mat_man_tib, by="IlmnID", suffix=c("_SRCX","_SRCY"))

# QC Conflicts & Remove Redundancy::

# QC: Must be equal::
#     - AlleleA_ProbeSeq_SRCX
#     - AlleleB_ProbeSeq_SRCX
#     - Forward_Sequence_SRCX
qc_fail_tib <- new_man_tib %>% dplyr::filter(AlleleA_ProbeSeq_SRCX!=AlleleA_ProbeSeq_SRCY,
                                             AlleleB_ProbeSeq_SRCX!=AlleleB_ProbeSeq_SRCY,
                                             Forward_Sequence_SRCX!=Forward_Sequence_SRCY)
qc_fail_nrows <- base::nrow(qc_fail_tib)
if (qc_fail_nrows!=0) {
  print(qc_fail_tib)
  stop(glue::glue("{RET}[QC-FAILURE]: ERROR: qc_fail_nrows={qc_fail_nrows}.{RET}{RET}"))
  q()
}

# Remove Redundancy::
fin_man_tib <- new_man_tib %>% 
  dplyr::rename(Name=Name_SRCY,
                CGN=CGN_SRCY,
                AddressA_ID=AddressA_ID_SRCX,
                AlleleA_ProbeSeq=AlleleA_ProbeSeq_SRCX,
                AddressB_ID=AddressB_ID_SRCX,
                AlleleB_ProbeSeq=AlleleB_ProbeSeq_SRCX,
                Forward_Sequence=Forward_Sequence_SRCX) %>%
  dplyr::select(-c(Name_SRCX,CGN_SRCX,AddressA_ID_SRCY,AlleleA_ProbeSeq_SRCY,
                   AddressB_ID_SRCY,AlleleB_ProbeSeq_SRCY,
                   Forward_Sequence_SRCY))

# Review Calculation Conflicts::
#  - Infinium_Design_Type
#  - Next_Base
#  - Color_Channel
#
qc_inf_fail_nrows <- fin_man_tib %>% 
  dplyr::filter(Infinium != Infinium_Design_Type) %>% base::nrow()
qc_nxb_fail_nrows <- fin_man_tib %>% 
  dplyr::filter(!is.na(Next_Base_Ref), 
                stringr::str_to_upper(cmpl(Next_Base_Ref)) != Next_Base) %>% base::nrow()
qc_col_fail_nrows <- fin_man_tib %>% 
  dplyr::filter(Color_Channel_Ref != Color_Channel) %>% base::nrow()

cat(glue::glue("[Stats]: qc_inf_fail_nrows={qc_inf_fail_nrows}.{RET}"))
cat(glue::glue("[Stats]: qc_nxb_fail_nrows={qc_nxb_fail_nrows}.{RET}"))
cat(glue::glue("[Stats]: qc_col_fail_nrows={qc_col_fail_nrows}.{RET}"))

# Lastly Add Design QC and Sort Final Manifest::
fin_man_tib <- fin_man_tib %>%
  dplyr::mutate(Design_QC=case_when(Man_TarMatch & Man_MisMatch<=opt$MinMisMatch ~ TRUE, TRUE ~ FALSE)) %>%
  dplyr::select(IlmnID,Name,AddressA_ID,AlleleA_ProbeSeq,AddressB_ID,AlleleB_ProbeSeq,
                Design_QC,CGN,Strand_FR,Strand_TB,Strand_CO,
                Infinium,Next_Base_Ref,Color_Channel_Ref,
                Forward_Sequence,DesSeqN,DesBscU,DesBscM,DesBscD,
                PRB1_U,PRB1_M,PRB2_D,
                NXB_U,NXB_M,NXB_D,
                CPN_U,CPN_M,CPN_D,
                TAR_U,TAR_M,TAR_D,
                BOD_U,BOD_M,BOD_D,
                END_U,END_M,END_D,
                dplyr::everything())

fin_man_tib %>% dplyr::group_by(Design_QC) %>% summarise(n()) %>% print()

readr::write_csv(fin_man_tib, opt$out_man_csv)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                                TBD LIST::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# 0. Investigate RS Infinium I Failures Due to lack of SNP (or C/G base missing)
# 1. [***SCRATCH]:: Kepp Reference!!! RevCmp Whole Sequece and Change Substring Indexes to 
#    - Put All Sequences in Design Probe Sequence Orientation 
# 2. [*DONE]:: improbe-design.tsv -> improbe-design-manifest.csv
#    - Compare Strands TB/FR/CO, etc. 
# 3. [*DONE] Command Line Parameters
# 4. Add Top/Bot Calculation
# 5. Read Mixed Probe Type Manifest (CG,CH,RS, etc.)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                               Testing LIST::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# 1. Plot CH/RS Failures across samples
#    - Use Design_QC flag as standard for plotting

# 2. improbe-design.tsv
# 3. Full EPIC Testing
# 4. NZT Testing
# 5. 1um 15k Testing
# 6. mm10 Wanding Testing

# Core Manifest Fields
# dplyr::select(IlmnID, Name, AddressA_ID, AlleleA_ProbeSeq, AddressB_ID, AlleleB_ProbeSeq,
#               Infinium_Design_Type, Next_Base, Color_Channel, Forward_Sequence, Genome_Build, CHR, MAPINFO)


# End of file
