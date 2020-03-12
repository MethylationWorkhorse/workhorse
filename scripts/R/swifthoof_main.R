
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                              Source Packages::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

rm(list=ls(all=TRUE))

suppressWarnings(suppressPackageStartupMessages( base::require("optparse",quietly=TRUE) ))

suppressWarnings(suppressPackageStartupMessages( base::require("tidyverse") ))
suppressWarnings(suppressPackageStartupMessages( base::require("stringr") ))
suppressWarnings(suppressPackageStartupMessages( base::require("glue") ))

suppressWarnings(suppressPackageStartupMessages( base::require("matrixStats") ))
suppressWarnings(suppressPackageStartupMessages( base::require("scales") ))
suppressWarnings(suppressPackageStartupMessages( base::require("grid") ))

# Parallel Computing Packages
suppressWarnings(suppressPackageStartupMessages( base::require("doParallel") ))

# Load sesame::
suppressWarnings(suppressPackageStartupMessages( base::require("sesame") ))

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

par$retData     <- FALSE

# Illumina based directories::
par$macDir <- '/Users/bbarnes/Documents/CustomerFacing'
par$lixDir <- '/illumina/scratch/darkmatter/data'

# Executables::
opt$Rscript <- NULL

# Directories::
opt$outDir     <- NULL
opt$datDir     <- NULL
opt$idatsDir   <- NULL

# Optional Files::
opt$subManifest  <- FALSE
opt$manifestPath <- NULL
opt$addressPath  <- NULL
opt$auto_sam_csv <- NULL

# Platform/Method Options::
opt$platform  <- 'EPIC'
opt$manifest  <- 'B4'

# Run Options::
opt$fresh       <- FALSE
opt$buildSubDir <- FALSE
opt$autoDetect  <- FALSE
opt$workflows   <- 'ind'

# Output Options::
opt$loadIdat    <- FALSE
opt$saveIdat    <- FALSE

opt$loadSsets   <- FALSE
opt$saveSsets   <- FALSE
opt$saveRawSset <- FALSE
opt$lightFootPrint <- FALSE

opt$writeSset   <- FALSE
opt$writeSsum   <- FALSE
opt$writeCalls  <- FALSE
opt$writeSsheet <- FALSE
opt$writeAuto   <- FALSE

# Reporting Options::
opt$sigs_sum_field <- NULL

# Threshold Options::
opt$minNegPval   <- 0.02
opt$minOobPval   <- 0.1
opt$minNegPerc   <- 98
opt$minOobPerc   <- 90
opt$minDeltaBeta <- 0.2

opt$percisionSigs <- 1
opt$percisionBeta <- 4
opt$percisionPval <- 6

# Parallel/Cluster Options::
opt$single   <- FALSE
opt$parallel <- FALSE
opt$cluster  <- FALSE

# Plotting Options::
opt$plotSset  <- FALSE
opt$plotCalls <- FALSE
opt$plotAuto  <- FALSE

opt$plotFormat <- 'pdf'
opt$plotFormat <- 'png'

opt$dpi <- 72
opt$dpi <- 120

opt$plotMax <- 10000
opt$plotSub <- 5000

# Verbosity Options::
opt$verbosity <- 3

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                               Parse Options::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
args.dat <- commandArgs(trailingOnly = FALSE)
if (args.dat[1]=='RStudio') {
  
  # Default Parameters for local Mac::
  par$runMode    <- args.dat[1]
  par$dirRelPath <- file.path(par$macDir, 'workhorse/scripts')
  par$prgmTag    <- 'swifthoof_main'
  par$exeLocPath <- file.path(par$dirRelPath, 'R', paste0(par$prgmTag,'.R'))
  par$topDir     <- par$macDir
  
  # Default Options for local Mac::
  opt$Rscript  <- 'Rscript'
  opt$datDir   <- file.path(par$topDir, 'workhorse/dat')

  opt$subManifest  <- FALSE
  
  # Run Options::
  par$retData     <- TRUE
  opt$fresh       <- FALSE
  opt$buildSubDir <- FALSE
  opt$autoDetect  <- TRUE
  
  opt$workflows  <- 'raw,ind,ndi,din'
  opt$workflows  <- 'ind'

  # Output Options::
  opt$loadIdat    <- TRUE
  opt$saveIdat    <- TRUE
  
  opt$loadSsets   <- TRUE
  opt$saveSsets   <- TRUE
  opt$saveRawSset <- TRUE
  
  opt$saveIdat    <- FALSE
  opt$saveSsets   <- FALSE
  opt$saveRawSset <- FALSE

  opt$writeSset   <- FALSE
  opt$writeSsum   <- FALSE
  opt$writeCalls  <- TRUE
  opt$writeSsheet <- TRUE
  opt$writeAuto   <- FALSE
  
  # Output all::
  output_all <- TRUE
  output_all <- FALSE
  if (output_all) {
    opt$loadIdat    <- TRUE
    opt$saveIdat    <- TRUE
    
    opt$loadSsets   <- TRUE
    opt$saveSsets   <- TRUE
    opt$saveRawSset <- TRUE
    
    opt$writeSset   <- TRUE
    opt$writeSsum   <- TRUE
    opt$writeCalls  <- TRUE
    opt$writeSsheet <- TRUE
    opt$writeAuto   <- FALSE
  }
  
  opt$plotAuto    <- FALSE
  
  opt$percisionSigs <- 1
  opt$percisionBeta <- 4
  opt$percisionPval <- 6

  # Parallel/Cluster Options::
  opt$single   <- TRUE
  opt$parallel <- FALSE
  opt$cluster  <- FALSE
  
  # Different Experiments to run locally::
  #
  opt$expRunStr  <- 'EXP5'
  opt$expChipNum <- '203452220020'
  opt$idatsDir <- file.path(par$topDir, 'report-20191112/idats', opt$expRunStr)
  
  opt$expRunStr  <- 'BadDELTA'
  opt$expChipNum <- '203319730003'
  opt$idatsDir   <- file.path(par$topDir, paste('idats', opt$expRunStr, sep='_') )

  opt$expRunStr  <- 'DeltaBetaCore'
  opt$expChipNum <- '202761400007'
  opt$expChipNum <- '204275020006'
  opt$idatsDir   <- '/Users/bbarnes/Documents/Projects/workhorse/idats_DeltaBetaCore'

  opt$expRunStr  <- 'NA12878'
  opt$expChipNum <- '200348350023'
  opt$expChipNum <- '200348350002' # R05C01
  opt$idatsDir   <- file.path('/Users/bbarnes/Documents/Projects/workhorse', paste0('idats/idats_',opt$expRunStr), opt$expChipNum)

  opt$expRunStr  <- 'EX-CAM1'
  opt$expChipNum <- '202915460006'
  opt$idatsDir   <- file.path('/Users/bbarnes/Documents/Projects/workhorse', paste0('idats/idats_',opt$expRunStr), opt$expChipNum)
  
  opt$expRunStr  <- '24x1-EPIC'
  opt$expChipNum <- '203631770004'
  opt$expChipNum <- '204275020006' # _R12C02
  opt$expChipNum <- '203631770004'
  opt$idatsDir   <- file.path('/Users/bbarnes/Documents/Projects/workhorse', paste0('idats/idats_',opt$expRunStr), opt$expChipNum)
  
  opt$expRunStr  <- 'BETA-8x1-EPIC-Core'
  opt$expChipNum <- '202761400007'
  opt$idatsDir   <- file.path('/Users/bbarnes/Documents/Projects/workhorse', paste0('idats/idats_',opt$expRunStr), opt$expChipNum)

  opt$expRunStr  <- 'BETA-EX-SET1'
  opt$expChipNum <- '202915460006'
  opt$idatsDir   <- file.path('/Users/bbarnes/Documents/Projects/workhorse', paste0('idats/idats_',opt$expRunStr), opt$expChipNum)
  
  opt$outDir     <- file.path(par$topDir, 'workspace', par$prgmTag, par$runMode, opt$expRunStr)
  
  # Excalibur Only:: Original [obsolete]
  #  opt$auto_sam_csv <- file.path(opt$datDir, 'ref', 'AutoSampleDetection_EPIC-B4_BETA_EX_median_beta.csv.gz')
  
  # Original Core Betas (3T+6C)
  # opt$auto_sam_csv <- file.path(opt$datDir, 'ref', paste('AutoSampleDetection',pool,'Pval-Betas.csv.gz', sep='_'))
  
  # Try1 Excalibur:: No Detection P-value [obsolete]
  #  opt$auto_sam_csv <- '/Users/bbarnes/Documents/CustomerFacing/rds/refs/AutoSampleDetection_EPIC-B4_BETA_EX_median_beta_noPval.csv.gz'
  
  # Subsets::
  opt$auto_sam_csv <- file.path(opt$datDir, 'ref/AutoSampleDetection_EPIC-B4_8x1_pneg98_Median_beta_noPval_BETA-Dx.csv.gz')
  #  opt$auto_sam_csv <- file.path(opt$datDir, 'ref/AutoSampleDetection_EPIC-B4_8x1_pneg98_Median_beta_noPval_BETA-Zymo.csv.gz')
  #  opt$auto_sam_csv <- file.path(opt$datDir, 'ref/AutoSampleDetection_EPIC-B4_8x1_pneg98_Median_beta_noPval_DELTA-Zymo.csv.gz')

  # Try2 Excalibur:: No Detection P-value
  # opt$auto_sam_csv <- file.path(opt$datDir, '/ref/AutoSampleDetection_EPIC-B4_8x1_pneg98_Median_beta_noPval.csv.gz')
  
  opt$verbosity <- 3
  
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
    make_option(c("-i", "--idatsDir"), type="character", default=opt$idatsDir, 
                help="idats directory [default= %default]", metavar="character"),
    
    # Optional Files::
    make_option(c("--manifestPath"), type="character", default=opt$manifestPath,
                help="Path to manfifest (CSV) otherwise use dat [default= %default]", metavar="character"),
    make_option(c("--addressPath"), type="character", default=opt$addressPath,
                help="Path to address (RDS) otherwise use dat [default= %default]", metavar="character"),
    make_option(c("--subManifest"), action="store_true", default=opt$subManifest,
                help="Boolean variable to use subset manifest instead of subset. [default= %default]", metavar="boolean"),
    make_option(c("--auto_sam_csv"), type="character", default=opt$auto_sam_csv,
                help="Path to auto detect beta values (CSV) [default= %default]", metavar="character"),
    
    # Platform/Method Parameters::
    make_option(c("--platform"), type="character", default=opt$platform, 
                help="Forced platform [EPIC, 450k, 27k, NZT] otherwise auto-detect [default= %default]", metavar="character"),
    make_option(c("--manifest"), type="character", default=opt$manifest, 
                help="Forced manifest [B1, B2, B4] otherwise auto-detect [default= %default]", metavar="character"),
    
    # Run Parameters::
    make_option(c("--fresh"), action="store_true", default=opt$fresh,
                help="Boolean variable to build fresh version of database files [default= %default]", metavar="boolean"),
    make_option(c("--buildSubDir"), action="store_true", default=opt$buildSubDir,
                help="Boolean variable to build subdirectories based on Chip/BeadPool (for R&D purposes) [default= %default]", metavar="boolean"),
    make_option(c("--autoDetect"), action="store_true", default=opt$autoDetect,
                help="Boolean variable to auto detect reference samples. Must provide reference samples. [default= %default]", metavar="boolean"),
    make_option(c("--workflows"), type="character", default=opt$workflows,
                help="Order of operations comma seperated [ raw,ind,ndi,din ] [default= %default]", metavar="character"),
    # make_option(c("--sampleSheet"), type="character", default=opt$sampleSheet, 
    #             help="Target Sample Sheet containing samples/chips to ONLY analyze [default= %default]", metavar="character"),
    
    # Output Options::
    make_option(c("--loadIdat"), action="store_true", default=opt$loadIdat,
                help="Boolean variable to load existing IDAT from RDS file [default= %default]", metavar="boolean"),
    make_option(c("--saveIdat"), action="store_true", default=opt$saveIdat,
                help="Boolean variable to write IDAT RDS file [default= %default]", metavar="boolean"),
    
    make_option(c("--loadSsets"), action="store_true", default=opt$loadSsets,
                help="Boolean variable to load existing Signal Set from RDS file [default= %default]", metavar="boolean"),
    make_option(c("--saveSsets"), action="store_true", default=opt$saveSsets,
                help="Boolean variable to write Signal Set RDS file [default= %default]", metavar="boolean"),
    make_option(c("--saveRawSset"), action="store_true", default=opt$saveRawSset,
                help="Boolean variable to write Raw Signal Set RDS file [default= %default]", metavar="boolean"),
    make_option(c("--lightFootPrint"), action="store_true", default=opt$lightFootPrint,
                help="Boolean variable to NOT save any RDS files [default= %default]", metavar="boolean"),
    
    make_option(c("--writeSset"), action="store_true", default=opt$writeSset,
                help="Boolean variable to write Signal Set file [default= %default]", metavar="boolean"),
    make_option(c("--writeSsum"), action="store_true", default=opt$writeSsum,
                help="Boolean variable to write Signal Set Summary file [default= %default]", metavar="boolean"),
    make_option(c("--writeCalls"), action="store_true", default=opt$writeCalls,
                help="Boolean variable to write Calls (Pval/Beta) file [default= %default]", metavar="boolean"),
    make_option(c("--writeSsheet"), action="store_true", default=opt$writeSsheet,
                help="Boolean variable to Sample Sheet file [default= %default]", metavar="boolean"),
    make_option(c("--writeAuto"), action="store_true", default=opt$writeAuto,
                help="Boolean variable to write Auto-Detection Matricies (Pval/Beta) file [default= %default]", metavar="boolean"),

    # Reporting Options::
    make_option(c("--sigs_sum_field"), type="character", default=opt$sigs_sum_field, 
                help="Signal summary field in AutoSampleSheet [default= %default]", metavar="character"),
    
    # Threshold Options::
    make_option(c("--minNegPval"), type="double", default=opt$minNegPval, 
                help="Minimum passing detection p-value using Negative Controls [default= %default]", metavar="double"),
    make_option(c("--minOobPval"), type="double", default=opt$minOobPval,
                help="Minimum passing detection p-value using Negative Out-Of-Band [default= %default]", metavar="double"),
    
    make_option(c("--minNegPerc"), type="double", default=opt$minNegPerc, 
                help="Minimum percentage of loci passing detection p-value using Negative Controls to flag Requeue of sample. [default= %default]", metavar="double"),
    make_option(c("--minOobPerc"), type="double", default=opt$minOobPerc, 
                help="Minimum percentage of loci passing detection p-value using Out-Of-Band to flag Requeue of sample. [default= %default]", metavar="double"),
    
    make_option(c("--minDeltaBeta"), type="double", default=opt$minDeltaBeta,
                help="Minimum passing delta-beta. Used in AutoSampleSheet cacluclations [default= %default]", metavar="double"),
    
    make_option(c("--percisionSigs"), type="integer", default=opt$percisionSigs,
                help="Rounding percision for signal values in calls output files [default= %default]", metavar="double"),
    make_option(c("--percisionBeta"), type="integer", default=opt$percisionBeta,
                help="Rounding percision for beta values in calls output files [default= %default]", metavar="double"),
    make_option(c("--percisionPval"), type="integer", default=opt$percisionPval,
                help="Rounding percision for detection p-values in calls output files [default= %default]", metavar="double"),
    
    # Parallel/Cluster Parameters::
    make_option(c("--single"), action="store_true", default=opt$single, 
                help="Boolean variable to run a single sample on a single-core [default= %default]", metavar="boolean"),
    make_option(c("--parallel"), action="store_true", default=opt$parallel, 
                help="Boolean variable to run parallel on multi-core [default= %default]", metavar="boolean"),
    make_option(c("--cluster"), action="store_true", default=opt$cluster,
                help="Boolean variable to run jobs on cluster by chip [default= %default]", metavar="boolean"),
    
    # Plotting Options::
    make_option(c("--plotSset"), action="store_true", default=opt$plotSset,
                help="Boolean variable to plot intensity distributions for sset [default= %default]", metavar="boolean"),
    make_option(c("--plotCalls"), action="store_true", default=opt$plotCalls,
                help="Boolean variable to plot detection p-values and beta distributions [default= %default]", metavar="boolean"),
    make_option(c("--plotAuto"), action="store_true", default=opt$plotAuto,
                help="Boolean variable to plot Auto-Detection Matricies (Pval/Beta) file [default= %default]", metavar="boolean"),
    
    make_option(c("--plotFormat"), type="character", default=opt$plotFormat, 
                help="Plotting output format [default= %default]", metavar="character"),
    make_option(c("--dpi"), type="double", default=opt$dpi, 
                help="DPI for plot images Plotting [default= %default]", metavar="double"),
    
    make_option(c("--plotMax"), type="double", default=opt$plotMax, 
                help="Max Sample Display Count for Plotting [default= %default]", metavar="double"),
    make_option(c("--plotSub"), type="double", default=opt$plotSub, 
                help="Sub Sample Display Count for Plotting [default= %default]", metavar="double"),
    
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
if (opt$lightFootPrint) {
  opt$saveRawSset <- FALSE
  opt$saveSsets   <- FALSE
}

if (is.null(par$runMode) || is.null(par$prgmTag) || is.null(par$dirRelPath)) {
  if (par$runMode=='CommandLine') print_help(opt_parser)
  dplyr::bind_rows(par) %>% tidyr::gather("Params", "Value") %>% as.data.frame() %>% print()
  base::stop("Null Parameters!\n\n")
}

if (is.null(opt$outDir)   || is.null(opt$idatsDir) || is.null(opt$datDir) ||
    is.null(opt$platform) || is.null(opt$manifest) ||
    is.null(opt$verbosity)) {
  if (par$runMode=='CommandLine') print_help(opt_parser)
  dplyr::bind_rows(opt) %>% tidyr::gather("Option", "Value") %>% as.data.frame() %>% print()
  base::stop("Null Options!\n\n")
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                              Source Files::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
par$source_dir <- file.path(par$dirRelPath, 'R/workflows')
source_files   <- list.files(path=par$source_dir, pattern='.R$', full.names=TRUE, recursive=TRUE)
for (sfile in source_files) base::source(sfile)

par$sourceA <- file.path(par$dirRelPath, 'R/swifthoof_functions.R')
if (!file.exists(par$sourceA)) stop(glue::glue("[{par$prgmTag}]: Source={par$sourceA} does not exist!{RET}"))
base::source(par$sourceA)

if (is.null(opt$manifestPath))
  opt$manifestPath <- file.path(opt$datDir, 'manifest', paste0(opt$platform,'-',opt$manifest,'.manifest.sesame.cpg-sorted.csv.gz') )
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

workflows_vec <- NULL
if (!is.null(opt$workflows)) workflows_vec <- opt$workflows %>% str_split(pattern=',', simplify=TRUE) %>% as.vector()

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#             Select Chips from idats and/or Target Manifest::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
chipPrefixes <- sesame::searchIDATprefixes(opt$idatsDir)

# if (!is.null(opt$sampleSheet)) { tar_ss_tib <- readr::read_csv(opt$sampleSheet) }

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                                  Main::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
if (opt$cluster) {
  cat(glue::glue("[{par$prgmTag}]: Launching Chips in Cluster Mode!"),"\n", sep='')
  
  par$lan_exe <- ''
  par$isLinux <- FALSE
  if (dir.exists(par$lixDir)) {
    par$isLinux <- TRUE
    par$lan_exe <- 'qsub -cwd -pe threaded 16 -l excl=true -N'
    if (dir.exists(par$macDir)) stop(glue::glue("[{par$prgmTag}]: Linux/Mac directories exist???{RET}{RET}"))
  }
  prefixTib <- prefixesToChipTib(chipPrefixes) %>% split(.$barcode)
  
  par$shellDir <- file.path(opt$outDir, 'shells')
  if (!dir.exists(par$shellDir)) dir.create(par$shellDir, recursive=TRUE)
  
  for (chipName in names(prefixTib)) { # break }
    runShell <- file.path(par$shellDir, paste0('run_',par$prgmTag,'_',chipName,'.sh'))
    lanShell <- file.path(par$shellDir, paste0('lan_',par$prgmTag,'_',chipName,'.sh'))
    
    # Remove Cluster Option
    cmd_bool <- opt %>% bind_rows() %>% gather("Options", "Value") %>% 
      dplyr::filter(Value=='TRUE') %>% 
      dplyr::filter(Options!='cluster') %>% 
      dplyr::select(Options) %>% 
      dplyr::mutate(Options=paste0('--',Options)) %>% 
      dplyr::pull() %>% paste(collapse=" ")
    
    # Add ChipName to idat Directory
    cmd_strs <- opt %>% bind_rows() %>% gather("Options", "Value") %>%
      dplyr::filter(Value!='TRUE' & Value!='FALSE' & Options!='prgmPath') %>%
      dplyr::mutate(Value=case_when(Options=='idatsDir' ~ paste0(Value,'/',chipName), TRUE ~ Value),
                    Value=case_when(Options=='outDir' ~ paste0(Value,'/',chipName), TRUE ~ Value) ) %>%
      dplyr::mutate(Options=paste0('--',Options,'=',Value)) %>%
      dplyr::select(Options) %>% 
      dplyr::pull() %>% paste(collapse=" ")
    
    cmd_full <- paste(opt$Rscript, par$exeLocPath, cmd_bool, cmd_strs,"\n", sep=' ')
    readr::write_file(cmd_full, path=runShell)
    Sys.chmod(runShell, mode="0777")
    
    # Add cluster execute if avialbel (i.e. linux)
    cmd_lan <- paste0(runShell,RET, sep='')
    if (par$isLinux)
      cmd_lan <- paste(par$lan_exe, paste('imWH',chipName, sep='_'), runShell,RET, sep=' ')
    readr::write_file(cmd_lan, path=lanShell)
    Sys.chmod(lanShell, mode="0777")
    
    # Launch Script
    cat(glue::glue("[{par$prgmTag}]: Launching chip={chipName}, shell={lanShell}!"),"\n", sep='')
    system(lanShell)
    
    if (opt$single) break
  }

} else {
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                             Load Manifest::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  man_tib <- loadManifestSource(opt$manifestPath, verbose=opt$verbosity,tt=pTracker)
  add_tib <- loadAddressSource(opt$addressPath, man_tib, fresh=opt$fresh, split=FALSE, 
                               verbose=opt$verbosity,tt=pTracker)
  # TBD:: Concern about overwriting add_tib if fresh and in cluster mode...
  
  auto_opt_tib <- NULL
  auto_ref_tib <- NULL
  auto_can_tib <- NULL
  
  if (opt$autoDetect) {
    
    stime <- system.time({
      auto_sam_tib <- suppressMessages(suppressWarnings(readr::read_csv(opt$auto_sam_csv) ))
    })
    pTracker$addTime(stime,'loadAutoSamples_New')
  }

  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                             Process Chip::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  chipTimes <- NULL
  
  if (opt$paralle) {
    funcTag <- 'sesamizeSingleSample'
    par$retData <- FALSE
    chipTimes <- foreach (prefix=names(chipPrefixes), .inorder=T, .final = function(x) setNames(x, names(chipPrefixes))) %dopar% {
      val <- NULL
      try_str <- ''
      val = tryCatch({
        try_str <- 'Pass'
        sesamizeSingleSample(prefix=chipPrefixes[[prefix]], man=man_tib, add=add_tib, ref=auto_sam_tib, opt=opt, 
                             retData=par$retData, workflows=workflows_vec)
      }, warning = function(w) {
        try_str <- paste('warning',funcTag, sep='-')
        NA
      }, error = function(e) {
        try_str <- paste('error',funcTag, sep='-')
        NA
      }, finally = {
        try_str <- paste('cleanup',funcTag, sep='-')
        NA
      })
      cat(glue::glue("[{par$prgmTag}] parallelFunc={funcTag}: try_str={try_str}. Done.{RET}{RET}"))

      val
    }
  } else {
    
    for (prefix in names(chipPrefixes)) {
      # chipTimes[[prefix]] <-
      #   sesamizeSingleSample(prefix=chipPrefixes[[prefix]], man=man_tib, add=add_tib, autoRef=auto_can_tib, opt=opt, retData=par$retData)
      
      # AnnotationHub( hub=getAnnotationHubOption("URL"),
      #                cache=getAnnotationHubOption("CACHE"),
      #                proxy=getAnnotationHubOption("PROXY"),
      #                localHub=FALSE)
      
      opt$verbosity <- 6
      rdat <- sesamizeSingleSample(prefix=chipPrefixes[[prefix]], man=man_tib, add=add_tib, ref=auto_sam_tib, opt=opt, 
                                   retData=par$retData, workflows=workflows_vec)

      #  GCT Score Example::
      # rdat$raw_sset %>% sesame::bisConversionControl()
      #  'BETA-8x1-EPIC-Core', prefix="202761400007_R01C01", GCT=1.173306 (raw)/ 1.162016 (cur==IND)
      #  'BETA-EX-SET1', prefix="202915460006_R01C01", GCT=1.453281 (raw)/ 1.354524 (cur=IND)
      
      # Output::
      #   1. SampleSheet
      #   2. Calls
      #   3. SignalSet
      #   4. Variants
      
      # Excalibur::
      #   1. Build all calls (no auto-detect)
      #   2. Build Mean Reference Samples
      #   3. Re-build all calls (auto-detect)
      
      
      # mutateSSET
      # sset2mtc -> mtc2sum
      # sset2sig -> sig2sum
      # sset2sss

      # 1. Generate Calls files
      
      # Build Average Reference Sets (pval-mask + row-averages)
      
      # Swift Workflow:: One Pipeline
      #  1. Raw steps -> pval/betas
      #  2. Pipeline -> pvals/betas
      #  3. Build new manifest -> Pipeline -> pvals/betas
      #  4. Delta Pre/Post position (Rcpp) [Requires sorted array]
      #  5. 
      
      # Signature of next base 0001, 0010, 1010, 01010, etc.
      # Signature of Query SNP CpG, is there difference between 
      #   Inf1: Block extention and miss-incorporation
      #   Inf2: Block extension and miss-incorporation
      
      if (opt$single) break
    }
  }
  
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                                Finished::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
sysTime <- Sys.time()
cat(glue::glue("[{par$prgmTag}]: Finished(time={sysTime}){RET}{RET}"))

# End of file
