
rm(list=ls(all=TRUE))

suppressWarnings(suppressPackageStartupMessages(require("optparse",quietly=TRUE)))

suppressWarnings(suppressPackageStartupMessages(require("tidyverse")) )
suppressWarnings(suppressPackageStartupMessages(require("stringr")) )
suppressWarnings(suppressPackageStartupMessages(require("glue")) )
suppressWarnings(suppressPackageStartupMessages(require("scales")) )
suppressWarnings(suppressPackageStartupMessages(require("matrixStats")) )
suppressWarnings(suppressPackageStartupMessages(require("grid")) )

# Plotting Extras
suppressWarnings(suppressPackageStartupMessages(require("GGally")) )
suppressWarnings(suppressPackageStartupMessages(require("hexbin")) )

# Parallel Computing Packages
suppressWarnings(suppressPackageStartupMessages(require("doParallel")) )

# c++ Packages
suppressPackageStartupMessages(require("Rcpp"))

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                              Global Params::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

doParallel::registerDoParallel()
num_cores   <- detectCores()
num_workers <- getDoParWorkers()

COM <- ","
TAB <- "\t"
RET <- "\n"

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                      Define Default Params and Options::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
par <- NULL
opt <- NULL

# Illumina based directories::
par$macDir <- '/Users/bbarnes/Documents/CustomerFacing'
par$lixDir <- '/illumina/scratch/darkmatter/data'

# Executables::
opt$Rscript <- NULL

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#  Experiment Input Variables::
opt$expNums <- NA
opt$expVars <- NULL

opt$filterKeys <- NA
opt$filterVals <- NA

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#  Directories::
opt$outDir <- NULL
opt$datDir <- NULL
opt$samDir <- NULL

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#  Optional Files::
opt$manifestPath <- NULL
opt$addressPath  <- NULL

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#  Platform/Method Parameters::
opt$platform  <- 'EPIC'
opt$manifest  <- 'B4'
opt$genome    <- 'hg19'

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#  Run Parameters::
opt$writeSigs  <- FALSE
opt$writeCall  <- FALSE

opt$fresh      <- FALSE
opt$autoDetect <- FALSE
opt$minPval    <- 0.02
opt$minDelta   <- 0.2

opt$auto_pval_field <- 'NDI_negs_pval'
opt$auto_beta_field <- 'NDI_beta'

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#  Parallel/Cluster Parameters::
opt$single   <- FALSE
opt$parallel <- FALSE
opt$cluster  <- FALSE
opt$clean    <- FALSE

opt$plotControls <- FALSE
opt$chipSummary  <- FALSE

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#  Remaining Plot Parameters::
#
opt$pvalPassKey <- 'PassDetpNegs_Percent_CG'
opt$pvalPassMin <- 96

opt$poobMinPval <- 0.2
opt$negsMinPval <- 0.02

opt$plotMax <- 3
opt$plotSub <- 10000
opt$plotSub <- 5000

opt$plotSpread  <- 'top'
opt$plotSpread  <- 'bot'
opt$plotSpread  <- 'mid'

opt$plotType   <- 'points'
opt$plotType   <- 'density'
opt$plotType   <- 'auto'

opt$plotFormat <- 'pdf'
opt$plotFormat <- 'png'

opt$dpi <- 72
opt$dpi <- 120

opt$verbosity <- 3


# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                               Parse Options::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
args.dat <- commandArgs(trailingOnly = FALSE)
if (args.dat[1]=='RStudio') {
  
  par$runMode    <- args.dat[1]
  par$dirRelPath <- file.path(par$macDir, 'workhorse/scripts')
  par$prgmTag    <- 'lighthoof_plotByExperiments'
  par$exeLocPath <- file.path(par$dirRelPath, 'mustang', paste0(par$prgmTag,'.R'))
  
  opt$Rscript  <- 'Rscript'
  par$topDir   <- par$macDir
  opt$datDir   <- file.path(par$topDir, 'dat')
  opt$samDir   <- file.path(par$topDir, 'sampleSheets/report-20191112/individual')
  
  opt$expDirName <- 'EXP6'
  opt$expDirName <- 'EXP5'
  opt$outDir   <- file.path(par$topDir, 'workspace', par$prgmTag, opt$expDirName)
  opt$buildDir <- file.path(par$topDir, 'workspace', 'lighthoof_main/builds', opt$expDirName)

  opt$expDirName <- 'BadDELTA'
  opt$outDir   <- file.path(par$topDir, 'workspace', par$prgmTag, opt$expDirName)
  opt$buildDir <- file.path(par$topDir, 'fromCluster/builds', opt$expDirName)

  opt$expDirName <- 'DeltaBetaCore_sesame_v1.7'
  opt$outDir   <- file.path(par$topDir, 'workspace', par$prgmTag, opt$expDirName)
  opt$buildDir <- file.path(par$topDir, 'workspace/docker', opt$expDirName)
  
  opt$expDirName <- 'DeltaBetaCore'
  opt$outDir   <- file.path(par$topDir, 'workspace', par$prgmTag, opt$expDirName)
  opt$buildDir <- file.path(par$topDir, 'workspace/lighthoof_main/builds', opt$expDirName)

  opt$writeSigs  <- FALSE
  opt$writeCall  <- FALSE
  
  opt$fresh      <- FALSE
  opt$autoDetect <- FALSE
  opt$cluster    <- FALSE
  opt$parallel   <- FALSE
  opt$single     <- FALSE
  
} else {
  par$runMode    <- 'CommandLine'
  par$exeLocPath <- base::substring(args.dat[grep("--file=", args.dat)], 8)
  par$dirLocPath <- base::dirname(par$exeLocPath)
  par$dirRelPath <- base::dirname(base::normalizePath(par$dirLocPath) )
  par$prgmTag    <- base::sub('\\.R$', '', base::basename(par$exeLocPath))
  
  args.dat <- commandArgs(trailingOnly = TRUE)
  option_list = list(
    # Executables::
    make_option(c("--Rscript"), type="character", default=opt$Rscript, 
                help="Rscript path [default= %default]", metavar="character"),
    
    # Directories::
    make_option(c("-o", "--outDir"), type="character", default=opt$outDir, 
                help="Output directory [default= %default]", metavar="character"),
    make_option(c("-d", "--datDir"), type="character", default=opt$datDir, 
                help="Dat Directory [default= %default]", metavar="character"),
    make_option(c("-i", "--samDir"), type="character", default=opt$samDir, 
                help="Human currated sample sheet directory [default= %default]", metavar="character"),
    
    # Optional Files::
    make_option(c("--manifestPath"), type="character", default=opt$manifestPath,
                help="Path to manfifest (CSV) otherwise use dat [default= %default]", metavar="character"),
    make_option(c("--addressPath"), type="character", default=opt$addressPath,
                help="Path to address (RDS) otherwise use dat [default= %default]", metavar="character"),
    
    # Platform/Method Parameters::
    make_option(c("--platform"), type="character", default=opt$platform, 
                help="Forced platform [EPIC, 450k, 27k, NZT] otherwise auto-detect [default= %default]", metavar="character"),
    make_option(c("--manifest"), type="character", default=opt$manifest, 
                help="Forced manifest [B1, B2, B4] otherwise auto-detect [default= %default]", metavar="character"),
    make_option(c("--genome"), type="character", default=opt$genome, 
                help="Genome Build [default= %default]", metavar="character"),
    
    # Run Parameters::
    make_option(c("--writeSigs"), action="store_true", default=opt$writeSigs,
                help="Boolean variable to write Signal file [default= %default]", metavar="boolean"),
    make_option(c("--writeCall"), action="store_true", default=opt$writeCall,
                help="Boolean variable to write Calls (Pval/Beta) file [default= %default]", metavar="boolean"),
    
    make_option(c("--fresh"), action="store_true", default=opt$fresh,
                help="Boolean variable force fresh build, e.g. ignore all other boolean parameters [default= %default]", metavar="boolean"),
    make_option(c("--autoDetect"), action="store_true", default=opt$autoDetect,
                help="Boolean variable to attempt auto sample detect [default= %default]", metavar="boolean"),
    make_option(c("--minPval"), type="double", default=opt$minPval, 
                help="Minimum passing detection p-value. Used in AutoSampleSheet cacluclations [default= %default]", metavar="double"),
    make_option(c("--minDelta"), type="double", default=opt$minDelta, 
                help="Minimum passing delta-beta value. [default= %default]", metavar="double"),
    
    make_option(c("--auto_beta_field"), type="character", default=opt$auto_beta_field, 
                help="Auto-Sample Detection beta field [default= %default]", metavar="character"),
    make_option(c("--auto_pval_field"), type="character", default=opt$auto_pval_field, 
                help="Auto-Sample Detection p-val field [default= %default]", metavar="character"),
    
    # Parallel/Cluster Parameters::
    make_option(c("--single"), action="store_true", default=opt$single, 
                help="Boolean variable to run a single sample on a single-core [default= %default]", metavar="boolean"),
    make_option(c("--parallel"), action="store_true", default=opt$parallel, 
                help="Boolean variable to run parallel on multi-core [default= %default]", metavar="boolean"),
    make_option(c("--cluster"), action="store_true", default=opt$cluster,
                help="Boolean variable to run jobs on cluster by chip [default= %default]", metavar="boolean"),
    make_option(c("--clean"), action="store_true", default=opt$clean,
                help="Clean Output Directory before plotting [default= %default]", metavar="boolean"),
    
    # Experiment Options::
    make_option(c("--expNums"), type="character", default=opt$expNums, 
                help="Experiment Numbers in a single comma seperated string [default= %default]", metavar="character"),
    make_option(c("--expVars"), type="character", default=opt$expVars, 
                help="Experimental Variables in a single comma seperated string [default= %default]", metavar="character"),
    
    make_option(c("--pvalPassKey"), type="character", default=opt$pvalPassKey,
                help="Detection P-value Type for Sample filtering (i.e. Negative Controls of Out-Of-Band) [default= %default]", metavar="character"),
    make_option(c("--pvalPassMin"), type="double", default=opt$pvalPassMin, 
                help="Detection P-value Minimum for Sample filtering [default= %default]", metavar="double"),
    
    make_option(c("--poobMinPval"), type="double", default=opt$poobMinPval, 
                help="Detection p-value filter for poob controls [default= %default]", metavar="double"),
    make_option(c("--negsMinPval"), type="double", default=opt$negsMinPval, 
                help="Detection p-value filter for negative controls [default= %default]", metavar="double"),
    
    # Plotting Parameters 
    make_option(c("--plotType"), type="character", default=opt$plotType, 
                help="Plot output type [auto, density, points] [default= %default]", metavar="character"),
    make_option(c("--plotSpread"), type="character", default=opt$plotSpread, 
                help="Plotting Sample Distribution Method [default= %default]", metavar="character"),
    make_option(c("--plotFormat"), type="character", default=opt$plotFormat, 
                help="Plotting output format [default= %default]", metavar="character"),
    make_option(c("--plotMax"), type="double", default=opt$plotMax, 
                help="Max Sample Display Count for Plotting [default= %default]", metavar="double"),
    make_option(c("--plotSub"), type="double", default=opt$plotSub, 
                help="Sub Sample Display Count for Plotting [default= %default]", metavar="double"),
    make_option(c("--dpi"), type="double", default=opt$dpi, 
                help="DPI for plot images Plotting [default= %default]", metavar="double"),
    make_option(c("--plotControls"), action="store_true", default=opt$plotControls,
                help="Boolean flag to plot all controls (slower) [default= %default]", metavar="boolean"),
    
    
    # Verbosity::
    make_option(c("-v", "--verbosity"), type="double", default=opt$verbosity, 
                help="0-5 (5 is very verbosity) [default= %default]", metavar="double")
  )
  opt_parser = OptionParser(option_list=option_list)
  opt = parse_args(opt_parser)
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                            Validate Options::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
if (is.null(par$runMode) || is.null(par$prgmTag) || is.null(par$dirRelPath)) {
  if (par$runMode=='CommandLine') print_help(opt_parser)
  dplyr::bind_rows(par) %>% tidyr::gather("Params", "Value") %>% as.data.frame() %>% print()
  base::stop("Null Parameters!\n\n")
}

if (is.null(opt$outDir)   || is.null(opt$datDir) || is.null(opt$samDir) ||
    is.null(opt$platform) || is.null(opt$manifest) || is.null(opt$minPval) ||
    is.null(opt$verbosity)) {
  if (par$runMode=='CommandLine') print_help(opt_parser)
  dplyr::bind_rows(opt) %>% tidyr::gather("Option", "Value") %>% as.data.frame() %>% print()
  base::stop("Null Options!\n\n")
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                              Source Files::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
par$sourceA <- file.path(par$dirRelPath, 'R/timeTracker.R')
par$sourceB <- file.path(par$dirRelPath, 'R/lighthoof_functions.R')
par$sourceC <- file.path(par$dirRelPath, 'R/lighthoof_sesame_functions.R')
par$sourceD <- file.path(par$dirRelPath, 'R/lighthoof_plotting_functions.R')
par$sourceE <- file.path(par$dirRelPath, 'R/Rcpp/CpgVariation_Rcpp.cpp')

if (!file.exists(par$sourceA)) stop(glue::glue("[{par$prgmTag}]: Source={par$sourceA} does not exist!{RET}"))
if (!file.exists(par$sourceB)) stop(glue::glue("[{par$prgmTag}]: Source={par$sourceB} does not exist!{RET}"))
if (!file.exists(par$sourceC)) stop(glue::glue("[{par$prgmTag}]: Source={par$sourceC} does not exist!{RET}"))
if (!file.exists(par$sourceD)) stop(glue::glue("[{par$prgmTag}]: Source={par$sourceD} does not exist!{RET}"))
if (!file.exists(par$sourceD)) stop(glue::glue("[{par$prgmTag}]: Source={par$sourceE} does not exist!{RET}"))

base::source(par$sourceA)
base::source(par$sourceB)
base::source(par$sourceC)
base::source(par$sourceD)
Rcpp::sourceCpp(par$sourceE)

if (is.null(opt$manifestPath))
  opt$manifestPath <- file.path(opt$datDir, 'manifest', paste0(opt$platform,'-',opt$manifest,'.manifest.sesame.add-sorted.csv.gz') )
if (is.null(opt$addressPath))
  opt$addressPath <- file.path(opt$datDir, 'manifest', paste0(opt$platform,'-',opt$manifest,'.manifest.address.cpg-sorted.csv.gz') )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                              Preprocessing::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
dplyr::bind_rows(opt) %>% tidyr::gather("Option", "Value") %>% print()
dplyr::bind_rows(par) %>% tidyr::gather("Params", "Value") %>% print()

pTracker <- timeTracker$new(verbose=opt$verbosity)

# opt$outDir <- file.path(opt$outDir, opt$platform, opt$manifest)
if (!dir.exists(opt$outDir)) dir.create(opt$outDir, recursive=TRUE)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                             Load Manifest::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
man_tib <- loadManifestSource(opt$manifestPath, verbose=opt$verbosity,tt=pTracker)
add_tib <- loadAddressSource(opt$addressPath, man_tib, fresh=FALSE, split=FALSE, 
                             verbose=opt$verbosity,tt=pTracker)
# TBD:: Concern about overwriting add_tib if fresh and in cluster mode...

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                 Load Auto Sample Sheets From Data::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

opt$verbosity <- 5
opt$pvalDetectFlag   <- TRUE
opt$pvalDetectMinKey <- 'CG_NDI_negs_pval_PassPerc'
opt$pvalDetectMinVal <- 96

opt$flagSampleDetect <- TRUE
opt$addressPath  <- TRUE
opt$filterRef    <- FALSE

auto_ss_tib <- loadAutoSampleSheets(dir=opt$buildDir, platform=opt$platform, manifest=opt$manifest,
                                    addSampleName=TRUE, 
                                    pvalDetectFlag=opt$pvalDetectFlag, pvalDetectMinKey=opt$pvalDetectMinKey, pvalDetectMinVal=opt$pvalDetectMinVal,
                                    flagSampleDetect=opt$flagSampleDetect, filterRef=opt$filterRef,
                                    addPaths=opt$addressPath, 
                                    verbose=opt$verbosity,tt=pTracker)

# file_list <- list.files(opt$samDir, pattern='SampleSheet.csv', full.names=TRUE)
# mman_ss_tib <- suppressMessages(suppressWarnings(lapply(file_list, readr::read_csv) )) %>% 
#   dplyr::bind_rows()

# Remove Unescessary Fields::
# mman_ss_tib <- mman_ss_tib %>% dplyr::select(-c("Catalog_Number", "Sample_Kit_Full_Name", "Sample_Kit_Type") )

mman_cnames <- NULL
# mman_cnames <- mman_ss_tib %>% dplyr::select(-Sentrix_Barcode) %>% names()
mand_cnames <- c('Bead_Pool', 'Chip_Format')

# ann_ss_tib <- mman_ss_tib %>% dplyr::inner_join(auto_ss_tib, by="Sentrix_Barcode")
# ann_ss_tib %>% dplyr::select(mand_cnames, mman_cnames)

ann_ss_tib  <- auto_ss_tib
uniq_cnames <- getUniqueFields(ann_ss_tib, c(mand_cnames, mman_cnames), verbose=opt$verbosity)
# ann_ss_tib %>% dplyr::select(uniq_cnames) %>% distinct()

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                    Group Samples:: i.e. Add Experiment Keys
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
if (length(uniq_cnames)!=0) {
  exp_ss_tibs <- ann_ss_tib %>% 
    tidyr::unite(Experiment_Key, uniq_cnames, sep='_', remove=FALSE) %>%
    tidyr::unite(Split_Key, Sample_Name, Experiment_Key, sep='_', remove=FALSE) %>%
    dplyr::select(Experiment_Key, Split_Key, dplyr::everything())
} else {
  exp_ss_tibs <- ann_ss_tib %>% dplyr::mutate(Split_Key=Sample_Name, Experiment_Key=Sample_Name) %>%
    dplyr::select(Experiment_Key, Split_Key, dplyr::everything())
  
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                            Plot Beta Matrix::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# Process by experiment first::
ss_exp_tib <- exp_ss_tibs %>% dplyr::arrange(-CG_NDI_negs_pval_PassPerc) %>% split(.$Experiment_Key)

for (expKey in names(ss_exp_tib)) {
  expKey <- '8x1_EPIC'
  ss_sam_tib <- ss_exp_tib[[expKey]] %>% split(.$Sample_Name)

  masked_tib <- NULL
  counts_vec <- NULL
  sample_vec <- names(ss_sam_tib)
  
  for (sample in names(ss_sam_tib)) {
    cur_tib <- ss_sam_tib[[sample]]
    
    maxSamples <- NULL
    # maxSamples <- 3
    betas <- loadCallFiles(files=cur_tib$Calls_Path, selKey='Probe_ID', 
                           datKey=opt$auto_beta_field, minKey=opt$auto_pval_field, minVal=opt$minPval, 
                           prefix=sample, max=maxSamples,
                           verbose=opt$verbosity+10,vt=0,tc=0,tt=pTracker)
    if (is.null(masked_tib)) { masked_tib <- betas
    } else { masked_tib <- masked_tib %>% dplyr::full_join(betas, by='Probe_ID') }
    counts_vec <- c(counts_vec, base::ncol(betas) - 1)
  }

  # Rcpp::sourceCpp(par$sourceE)
  loci_sum_tib <- getLociSampleSummary(masked_tib, verbose=opt$verbosity+1,vt=0,tt=pTracker)
  
  # Calculate Delta-Beta Distrbution::sd
  dbd_std <- loci_sum_tib %>% dplyr::select(ends_with('db_sd')) %>% 
    tidyr::gather('metric', 'deltaBeta')
  ggplot2::ggplot(dbd_std, aes(x=deltaBeta, color=metric, fill=metric)) +
    ggplot2::geom_density(alpha=0.3) +
    ggplot2::xlim(-0.01, 0.25)
  
  # Calculate Delta-Beta Distrbution::mean
  dbd_avg <- loci_sum_tib %>% dplyr::select(ends_with('db_avg')) %>% 
    tidyr::gather('metric', 'deltaBeta')
  ggplot2::ggplot(dbd_avg, aes(x=deltaBeta, color=metric, fill=metric)) +
    ggplot2::geom_density(alpha=0.3) +
    ggplot2::xlim(-0.01, 0.25)
  
  # Calculate CSS(Cluster Seperation Scores)::mean
  css_avg <- loci_sum_tib %>% 
    dplyr::mutate(HU_avg=H_avg-U_avg, MH_avg=M_avg-H_avg, MU_avg=M_avg-U_avg) %>%
    dplyr::select(HU_avg, MH_avg, MU_avg) %>% 
    tidyr::gather('metric', 'CSS')
  ggplot2::ggplot(css_avg, aes(x=CSS, color=metric, fill=metric)) +
    ggplot2::geom_density(alpha=0.5)

  # Calculate CSS(Cluster Seperation Scores)::median
  css_med <- loci_sum_tib %>% 
    dplyr::mutate(HU_med=H_med-U_med, MH_med=M_med-H_med, MU_med=M_med-U_med) %>%
    dplyr::select(HU_med, MH_med, MU_med) %>% 
    tidyr::gather('metric', 'CSS')
  
  ggplot2::ggplot(css_med, aes(x=CSS, color=metric, fill=metric)) +
    ggplot2::geom_density(alpha=0.5)
  
  # Cell Line Only Tib::
  cellLine_mat <- masked_tib %>% 
    dplyr::select(-starts_with('U_')) %>% 
    dplyr::select(-starts_with('H_')) %>% 
    dplyr::select(-starts_with('M_')) %>%
    dplyr::select_if(is_numeric) %>% as.matrix()
  
  # dplyr::bind_cols(
  cellLine_sum <- base::cbind(
    matrixStats::rowMeans2(cellLine_mat),
    matrixStats::rowMedians(cellLine_mat),
    matrixStats::rowSds(cellLine_mat),
    matrixStats::rowMads(cellLine_mat)
  ) %>% tibble::as_tibble() %>%
    purrr::set_names('beta_avg','beta_med','beta_sd','beta_mad') %>%
    dplyr::mutate(Probe_ID=dplyr::pull(masked_tib, 1)) %>%
    dplyr::select(Probe_ID, everything())
  
  break
}

if (FALSE) {
  # For now use standard...
  # ss_sam_tib <- exp_ss_tibs %>% dplyr::arrange(-CG_inf_negs_pval_PassPerc) %>% split(.$Sample_Name)
  # ss_sam_tib <- exp_ss_tibs %>% dplyr::arrange(-CG_RNDI_poob_pval_PassPerc) %>% split(.$Sample_Name)
  ss_sam_tib <- exp_ss_tibs %>% dplyr::arrange(-CG_NDI_negs_pval_PassPerc) %>% split(.$Sample_Name)
  cat(glue::glue("[{par$prgmTag}]: Starting Cross Sample Experiment Comparison(linear).{RET}"))
  
  for (sample in names(ss_sam_tib)) {
    sample <- 'HELA'
    
    exp_split_tibs <- ss_sam_tib[[sample]] %>% split(.$Experiment_Key)
    
    betas <- NULL
    for (expKey in names(exp_split_tibs)) { # break }
      cur_tib <- exp_split_tibs[[expKey]]
      
      maxSamples <- NULL
      betas <- loadCallFiles(files=cur_tib$Calls_Path, selKey='Probe_ID', 
                             datKey=opt$auto_beta_field, minKey=opt$auto_pval_field, minVal=opt$minPval, max=maxSamples,
                             verbose=opt$verbosity+10,vt=0,tc=0,tt=pTracker)
      
      btmp <- betas %>% dplyr::select(-Probe_ID) %>% head(n=100) %>% as.matrix()
      # C_lociVariation(btmp, 0.2)
      # Rcpp::sourceCpp(par$sourceE)
      
      
      break
    }
    break
    
    # retTibs <- plotBetaMatrix_bySample(tib=ss_sam_tib, sample=sample, field='Beta', minPval=opt$minPval, 
    #                                    pvalKey=opt$auto_beta_field, betaKey=opt$auto_pval_field,
    #                                    manifest=man_tib, outDir=opt$outDir, 
    #                                    spread=opt$plotSpread, outType=opt$plotType, dpi=opt$dpi, format=opt$plotFormat,
    #                                    max=2, maxCnt=opt$plotSub, retTibs=TRUE,
    #                                    verbose=opt$verbosity+1,vt=0)
    
    ret_dat <- plotBetaMatrix_bySample(tib=ss_sam_tib, sample=sample, field='Beta', minPval=opt$minPval, 
                                       pvalKey=opt$auto_pval_field, betaKey=opt$auto_beta_field,
                                       manifest=man_tib, outDir=opt$outDir, 
                                       spread=opt$plotSpread, outType=opt$plotType, dpi=opt$dpi, format=opt$plotFormat,
                                       max=2, maxCnt=opt$plotSub,
                                       verbose=opt$verbosity+1,vt=0)
    if (opt$single) break
    break
  }
}


# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                                Finished::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
sysTime <- Sys.time()
cat(glue::glue("[{par$prgmTag}]: Finished(time={sysTime}){RET}{RET}"))


# End of file
