
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                              Source Packages::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

rm(list=ls(all=TRUE))

suppressWarnings(suppressPackageStartupMessages(require("optparse",quietly=TRUE)))

suppressWarnings(suppressPackageStartupMessages(require("tidyverse")) )
suppressWarnings(suppressPackageStartupMessages(require("stringr")) )
suppressWarnings(suppressPackageStartupMessages(require("glue")) )
# suppressWarnings(suppressPackageStartupMessages(require("scales")) )
suppressWarnings(suppressPackageStartupMessages(require("matrixStats")) )
# suppressWarnings(suppressPackageStartupMessages(require("grid")) )

# Parallel Computing Packages
suppressWarnings(suppressPackageStartupMessages(require("doParallel")) )

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

par$retData <- FALSE

# Executables::
opt$Rscript <- NULL

# Directories::
opt$outDir     <- NULL
opt$datDir     <- NULL
opt$idatsDir   <- NULL

# Optional Files::
opt$manifestPath <- NULL
opt$addressPath  <- NULL

# Platform/Method Parameters::
# opt$platform  <- 'Auto'
# opt$manifest  <- 'Auto'
opt$platform  <- 'EPIC'
opt$manifest  <- 'B4'

# Run Parameters::
opt$sampleSheet <- NULL

opt$writeSigs  <- FALSE
opt$writeCall  <- FALSE
opt$writeAuto  <- FALSE
opt$plotAuto   <- FALSE

opt$fresh      <- FALSE
opt$freshAdd   <- FALSE
opt$autoDetect <- FALSE
opt$minPval    <- 0.02
opt$minDelta   <- 0.2

opt$negsMinCutoff   <- 96
opt$sigs_sum_field  <- 'avg'

# Parallel/Cluster Parameters::
opt$single   <- FALSE
opt$parallel <- FALSE
opt$cluster  <- FALSE

opt$plotMax <- 3
opt$plotSub <- 10000
opt$plotSub <- 5000

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
  par$dirRelPath <- file.path(par$macDir, 'git/workhorse/scripts')
  par$prgmTag    <- 'lighthoof_main'
  par$exeLocPath <- file.path(par$dirRelPath, 'R', paste0(par$prgmTag,'.R'))
  
  opt$Rscript  <- 'Rscript'
  par$topDir   <- par$macDir
  opt$datDir   <- file.path(par$topDir, 'dat')

  opt$writeSigs  <- FALSE
  opt$writeCall  <- TRUE
  opt$writeAuto  <- TRUE
  opt$plotAuto   <- TRUE
  
  # Fast Option::
  opt$writeCall  <- FALSE
  opt$writeAuto  <- FALSE
  opt$plotAuto   <- FALSE
  
  opt$fresh      <- FALSE
  opt$autoDetect <- TRUE
  opt$cluster    <- FALSE
  opt$parallel   <- FALSE
  opt$single     <- TRUE

  opt$expRunStr  <- 'EXP5'
  opt$expChipNum <- '203452220020'
  opt$idatsDir <- file.path(par$topDir, 'report-20191112/idats', opt$expRunStr)
  
  opt$expRunStr  <- 'BadDELTA'
  opt$expChipNum <- '203319730003'
  opt$idatsDir <- file.path(par$topDir, paste('idats', opt$expRunStr, sep='_') )
  
  opt$outDir   <- file.path(par$topDir, 'workspace', par$prgmTag, opt$expRunStr)

  if (!opt$cluster) {
    opt$outDir   <- file.path(opt$outDir, opt$expChipNum)
    opt$idatsDir <- file.path(opt$idatsDir, opt$expChipNum)
  }
  
  # Auto-Detect idats::
  # opt$autoDetect <- TRUE
  # opt$outDir   <- file.path(par$topDir, 'workspace', par$prgmTag, 'autoReference')
  # opt$idatsDir <- file.path(par$topDir, 'idats_refs')
  
  par$retData <- TRUE
  
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
    
    # Platform/Method Parameters::
    make_option(c("--platform"), type="character", default=opt$platform, 
                help="Forced platform [EPIC, 450k, 27k, NZT] otherwise auto-detect [default= %default]", metavar="character"),
    make_option(c("--manifest"), type="character", default=opt$manifest, 
                help="Forced manifest [B1, B2, B4] otherwise auto-detect [default= %default]", metavar="character"),
    
    # Run Parameters::
    make_option(c("--sampleSheet"), type="character", default=opt$sampleSheet, 
                help="Target Sample Sheet containing samples/chips to ONLY analyze [default= %default]", metavar="character"),
    
    make_option(c("--writeSigs"), action="store_true", default=opt$writeSigs,
                help="Boolean variable to write Signal file [default= %default]", metavar="boolean"),
    make_option(c("--writeCall"), action="store_true", default=opt$writeCall,
                help="Boolean variable to write Calls (Pval/Beta) file [default= %default]", metavar="boolean"),
    make_option(c("--writeAuto"), action="store_true", default=opt$writeAuto,
                help="Boolean variable to write Auto-Detection Matricies (Pval/Beta) file [default= %default]", metavar="boolean"),
    make_option(c("--plotAuto"), action="store_true", default=opt$plotAuto,
                help="Boolean variable to plot Auto-Detection Matricies (Pval/Beta) file [default= %default]", metavar="boolean"),

    make_option(c("--freshAdd"), action="store_true", default=opt$freshAdd,
                help="Boolean variable force fresh Address built from Manifest [default= %default]", metavar="boolean"),
    make_option(c("--fresh"), action="store_true", default=opt$fresh,
                help="Boolean variable force fresh build, e.g. ignore all other boolean parameters [default= %default]", metavar="boolean"),
    make_option(c("--autoDetect"), action="store_true", default=opt$autoDetect,
                help="Boolean variable to attempt auto sample detect [default= %default]", metavar="boolean"),
    make_option(c("--minPval"), type="double", default=opt$minPval, 
                help="Minimum passing detection p-value. Used in AutoSampleSheet cacluclations [default= %default]", metavar="double"),
    make_option(c("--minDelta"), type="double", default=opt$minDelta, 
                help="Minimum passing delta-beta. Used in AutoSampleSheet cacluclations [default= %default]", metavar="double"),
    
    make_option(c("--negsMinCutoff"), type="double", default=opt$negsMinCutoff, 
                help="Minimum passing detection p-value percentage of probes used to call requeue of sample from provider [default= %default]", metavar="double"),
    make_option(c("--sigs_sum_field"), type="character", default=opt$sigs_sum_field, 
                help="Signal summary field in AutoSampleSheet [default= %default]", metavar="character"),
    
    make_option(c("--plotFormat"), type="character", default=opt$plotFormat, 
                help="Plotting output format [default= %default]", metavar="character"),
    make_option(c("--plotMax"), type="double", default=opt$plotMax, 
                help="Max Sample Display Count for Plotting [default= %default]", metavar="double"),
    make_option(c("--plotSub"), type="double", default=opt$plotSub, 
                help="Sub Sample Display Count for Plotting [default= %default]", metavar="double"),
    make_option(c("--dpi"), type="double", default=opt$dpi, 
                help="DPI for plot images Plotting [default= %default]", metavar="double"),
    

    # Parallel/Cluster Parameters::
    make_option(c("--single"), action="store_true", default=opt$single, 
                help="Boolean variable to run a single sample on a single-core [default= %default]", metavar="boolean"),
    make_option(c("--parallel"), action="store_true", default=opt$parallel, 
                help="Boolean variable to run parallel on multi-core [default= %default]", metavar="boolean"),
    make_option(c("--cluster"), action="store_true", default=opt$cluster,
                help="Boolean variable to run jobs on cluster by chip [default= %default]", metavar="boolean"),
    
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

if (is.null(opt$outDir)   || is.null(opt$idatsDir) || is.null(opt$datDir) ||
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

if (!file.exists(par$sourceA)) stop(glue::glue("[{par$prgmTag}]: Source={par$sourceA} does not exist!{RET}"))
if (!file.exists(par$sourceB)) stop(glue::glue("[{par$prgmTag}]: Source={par$sourceB} does not exist!{RET}"))
if (!file.exists(par$sourceC)) stop(glue::glue("[{par$prgmTag}]: Source={par$sourceC} does not exist!{RET}"))
if (!file.exists(par$sourceD)) stop(glue::glue("[{par$prgmTag}]: Source={par$sourceD} does not exist!{RET}"))

base::source(par$sourceA)
base::source(par$sourceB)
base::source(par$sourceC)
base::source(par$sourceD)

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
  if (opt$fresh) opt$freshAdd <- TRUE
  
  man_tib <- loadManifestSource(opt$manifestPath, verbose=opt$verbosity,tt=pTracker)
  add_tib <- loadAddressSource(opt$addressPath, man_tib, fresh=opt$freshAdd, split=FALSE, 
                               verbose=opt$verbosity,tt=pTracker)
  # TBD:: Concern about overwriting add_tib if fresh and in cluster mode...
  
  auto_opt_tib <- NULL
  auto_ref_tib <- NULL
  auto_can_tib <- NULL
  
  if (opt$autoDetect) {
    pool <- 'BETA'

    auto_sam_csv <- file.path(opt$datDir, 'ref', paste('AutoSampleDetection',pool,'Pval-Betas.csv.gz', sep='_'))
    stime <- system.time({
      auto_sam_tib <- suppressMessages(suppressWarnings(readr::read_csv(auto_sam_csv) ))
    })
    pTracker$addTime(stime,'loadAutoSamples_New')
  }
  pTracker %>% print()
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                             Process Chip::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  chipTimes <- NULL
  
  if (opt$paralle) {
    chipTimes <- foreach (prefix=names(chipPrefixes), .inorder=T, .final = function(x) setNames(x, names(chipPrefixes))) %dopar% {
      sesamizeSingleSample(prefix=chipPrefixes[[prefix]], man=man_tib, add=add_tib, autoRef=auto_sam_tib, opt=opt, retData=par$retData)
    }
  } else {
    
    for (prefix in names(chipPrefixes)) {
      # chipTimes[[prefix]] <-
      #   sesamizeSingleSample(prefix=chipPrefixes[[prefix]], man=man_tib, add=add_tib, autoRef=auto_can_tib, opt=opt, retData=par$retData)
      
      rdat <- sesamizeSingleSample(prefix=chipPrefixes[[prefix]], man=man_tib, add=add_tib, autoRef=auto_sam_tib, opt=opt, retData=par$retData)

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
