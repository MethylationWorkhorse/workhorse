
rm(list=ls(all=TRUE))

suppressWarnings(suppressPackageStartupMessages(require("optparse",quietly=TRUE)))

suppressWarnings(suppressPackageStartupMessages(require("tidyverse")) )
suppressWarnings(suppressPackageStartupMessages(require("stringr")) )
suppressWarnings(suppressPackageStartupMessages(require("glue")) )
suppressWarnings(suppressPackageStartupMessages(require("scales")) )
suppressWarnings(suppressPackageStartupMessages(require("matrixStats")) )

# Parallel Computing Packages
suppressWarnings(suppressPackageStartupMessages(require("doParallel")) )

doParallel::registerDoParallel()
num_cores   <- detectCores()
num_workers <- getDoParWorkers()

# suppressWarnings(suppressPackageStartupMessages(require("sesame")) )
# suppressWarnings(suppressPackageStartupMessages(require("sesameData")) )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                             Default Options::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

opt <- NULL
opt$prgmTag <- 'workhorse'
opt$runMode  <- NULL

opt$autoDetect <- FALSE
opt$single     <- FALSE
opt$parallel   <- FALSE
opt$cluster    <- FALSE

opt$fresh       <- FALSE

opt$loadIDATS   <- FALSE
opt$loadSSETS   <- FALSE

opt$writeIDATS  <- FALSE
opt$writeSSETS  <- FALSE

opt$writePrbCSV <- TRUE
opt$writePrbRDS <- TRUE
opt$writeSSheet <- TRUE

# opt$writeMIDMAN <- FALSE
# opt$writeSESMAN <- FALSE
# opt$readfMIDMAN <- FALSE

# These write
opt$writeRDS    <- FALSE
opt$writeCSV    <- FALSE




opt$loadMAN     <- FALSE
opt$loadGRP     <- FALSE
opt$loadRDS     <- FALSE

# Save intermediate manifest csv files
opt$saveIDATS  <- FALSE
opt$saveMAN    <- FALSE
opt$saveGRP    <- FALSE

# Save intermediat RDS SSET files
opt$saveRDS     <- FALSE
opt$overRDS     <- FALSE

opt$useGRP      <- FALSE

opt$retSSET     <- FALSE
opt$retPRBS     <- FALSE

opt$Rscript  <- NULL
opt$topDir   <- NULL
opt$srcDir   <- NULL
opt$outDir   <- NULL
opt$idatsDir <- NULL
opt$funcScript <- 'workhorse_functions.R'
opt$sesaScript <- 'sesame_functions.R'
opt$idatScript <- 'idat_functions.R'

opt$sampleSheet <- NULL

opt$platform  <- 'EPIC-B4'
opt$build     <- 'hg19'
opt$method    <- 'both'

opt$poobMinPval <- 0.2
opt$negsMinPval <- 0.02

opt$poobMinCutoff_Stringent <- 85
opt$negsMinCutoff_Stringent <- 98

opt$poobMinCutoff_Relaxed <- 80
opt$negsMinCutoff_Relaxed <- 96


opt$verbosity <- 4

# QC Variable
tarSample <- NULL
prgmPath  <- NULL

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                                 Starting::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
ctime <- Sys.time()
cat(glue::glue("[{opt$prgmTag}]: Starting(time={ctime})"),"\n\n",sep='')

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                             Parse Options::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
args.dat <- commandArgs(trailingOnly = FALSE)
if (args.dat[1]=='RStudio') {
  opt$Rscript <- '/usr/local/bin/Rscript'
  opt$topDir  <- '/Users/bbarnes/Documents/Projects/workhorse'
  opt$srcDir  <- '/Users/bbarnes/Documents/Projects/workhorse/git/workhorse/scripts'
  # opt$srcDir  <- '/Users/bbarnes/Documents/git-test2/workhorse/scripts'
  prgmPath <- file.path(opt$srcDir, 'R', paste0(opt$prgmTag,'.R'))
  stopifnot(file.exists(prgmPath))
  
  opt$fresh <- TRUE
  
  opt$runMode  <- args.dat[1]
  opt$outDir   <- file.path(opt$topDir, 'workspace', opt$prgmTag)
  opt$idatsDir <- file.path(opt$topDir, 'idats')
  
  # opt$platform  <- 'NZT'
  opt$platform  <- 'EPIC-B4'
  if (opt$platform=='EPIC-B4') {
    opt$idatsDir  <- file.path(opt$topDir, 'idats', 'Open-ReferenceBETA')
    opt$idatsDir  <- file.path(opt$topDir, 'idats', 'ReferenceBETA')
    opt$idatsDir  <- file.path(opt$topDir, 'idats', 'ReferenceBETA', '201502830033')

    # tarSample <- '202761400009_R07C01'  # M
    # tarSample <- '202761400009_R06C01'  # H
    # tarSample <- '202761400009_R05C01'  # U
    # tarSample <- '203452220015_R06C01'  # H miscalled U, but with 842126 CpGs instead of full 862927
    opt$outDir <- file.path(opt$outDir, 'testing')
    # opt$verbosity <- 10
  } else if (opt$platform=='NZT') {
    opt$idatsDir  <- file.path(opt$topDir, 'idats', 'NZT')
  } else {
    stop("\nUNSUPPORTED PLATFORM!!!\n\n")
  }
  
  opt$autoDetect <- TRUE
  opt$single     <- TRUE
  opt$parallel   <- FALSE
  opt$cluster    <- FALSE

  opt$writeSSheet <- TRUE
  opt$writeMIDMAN <- FALSE
  
  opt$saveMAN    <- FALSE
  opt$saveGRP    <- FALSE
  opt$writeRDS    <- TRUE
  opt$writeCSV    <- FALSE
  
  opt$loadMAN     <- TRUE
  opt$loadGRP     <- FALSE
  opt$loadRDS     <- TRUE
  opt$saveRDS     <- TRUE
  opt$overRDS     <- FALSE
  
  opt$retSSET     <- FALSE
  opt$retPRBS     <- FALSE
} else {
  prgmPath <- substring(args.dat[grep("--file=", args.dat)], 8)
  
  opt$runMode <- 'CommandLine'
  opt$prgmTag <- sub('\\.R$', '', basename(substring(args.dat[grep("--file=", args.dat)], 8)))
  opt$srcDir  <- dirname(normalizePath(dirname(substring(args.dat[grep("--file=", args.dat)], 8)) ))
  # opt$topDir  <- dirname(opt$srcDir)
  
  args.dat <- commandArgs(trailingOnly = TRUE)
  option_list = list(
    make_option(c("--Rscript"), type="character", default=opt$Rscript, 
                help="Rscript path [default= %default]", metavar="character"),
    make_option(c("--prgmTag"), type="character", default=opt$prgmTag, 
                help="program name [default= %default]", metavar="character"),
    make_option(c("--runMode"), type="character", default=opt$runMode, 
                help="program run mode [default= %default]", metavar="character"),
    make_option(c("--sampleSheet"), type="character", default=opt$sampleSheet, 
                help="Sample Sheet provided by the user [default= %default]", metavar="character"),
    
    # Directories/Files::
    make_option(c("--topDir"), type="character", default=opt$topDir, 
                help="top data directory location [default= %default]", metavar="character"),
    make_option(c("--srcDir"), type="character", default=opt$srcDir, 
                help="source data directory location [default= %default]", metavar="character"),
    make_option(c("-o", "--outDir"), type="character", default=opt$outDir, 
                help="output directory [default= %default]", metavar="character"),
    make_option(c("-i", "--idatsDir"), type="character", default=opt$idatsDir, 
                help="idats directory [default= %default]", metavar="character"),
    
    # Standard Variables::
    make_option(c("--funcScript"), type="character", default=opt$funcScript, 
                help="Name of general functions script to soruce [default= %default]", metavar="character"),
    make_option(c("--sesaScript"), type="character", default=opt$sesaScript, 
                help="Name of Sesame functions script to soruce [default= %default]", metavar="character"),
    make_option(c("--idatScript"), type="character", default=opt$idatScript, 
                help="Name of Idat functions script to soruce [default= %default]", metavar="character"),
    
    # Platform/Method Parameters::
    make_option(c("--platform"), type="character", default=opt$platform, 
                help="Platform to use for analysis [default= %default]", metavar="character"),
    make_option(c("--build"), type="character", default=opt$build, 
                help="Genome build to use for analysis [default= %default]", metavar="character"),
    make_option(c("--method"), type="character", default=opt$method, 
                help="Method to use for analysis [default= %default]", metavar="character"),
    
    # Detection Parameters::
    make_option(c("--poobMinPval"), type="double", default=opt$poobMinPval, 
                help="Detection p-value filter for poob controls [default= %default]", metavar="double"),
    make_option(c("--negsMinPval"), type="double", default=opt$negsMinPval, 
                help="Detection p-value filter for negative controls [default= %default]", metavar="double"),
    
    make_option(c("--poobMinCutoff_Stringent"), type="double", default=opt$poobMinCutoff_Stringent, 
                help="Stringent detection p-value threshold percent for poob controls [default= %default]", metavar="double"),
    make_option(c("--negsMinCutoff_Stringent"), type="double", default=opt$negsMinCutoff_Stringent, 
                help="Stringent detection p-value threshold percent for negative controls [default= %default]", metavar="double"),
    
    make_option(c("--poobMinCutoff_Relaxed"), type="double", default=opt$poobMinCutoff_Relaxed, 
                help="Relaxed detection p-value threshold percent for poob controls [default= %default]", metavar="double"),
    make_option(c("--negsMinCutoff_Relaxed"), type="double", default=opt$negsMinCutoff_Relaxed, 
                help="Relaxed detection p-value threshold percent for negative controls [default= %default]", metavar="double"),
    
    # Boolean Variables::
    make_option(c("--fresh"), action="store_true", default=opt$fresh,
                help="Boolean variable force fresh build, e.g. ignore all other boolean parameters [default= %default]", metavar="boolean"),
    
    make_option(c("--autoDetect"), action="store_true", default=opt$autoDetect,
                help="Boolean variable to attempt auto sample detect [default= %default]", metavar="boolean"),
    
    #
    # Functional Start::
    #
    make_option(c("--saveMAN"), action="store_true", default=opt$saveMAN,
                help="Boolean to write Master Manifest [default= %default]", metavar="boolean"),
    
    make_option(c("--loadMAN"), action="store_true", default=opt$loadMAN,
                help="Boolean to load previous Manifest files [default= %default]", metavar="boolean"),
    make_option(c("--loadIDATS"), action="store_true", default=opt$loadIDATS,
                help="Boolean to load previously processed IDATS files [default= %default]", metavar="boolean"),
    make_option(c("--loadSSETS"), action="store_true", default=opt$loadSSETS,
                help="Boolean to load previously processed SSETS files [default= %default]", metavar="boolean"),
    
    make_option(c("--writeIDATS"), action="store_true", default=opt$writeIDATS,
                help="Boolean to write intermediate IDATS [default= %default]", metavar="boolean"),
    make_option(c("--writeSSETS"), action="store_true", default=opt$writeSSETS,
                help="Boolean to write intermediate SSETS [default= %default]", metavar="boolean"),
    make_option(c("--writePrbCSV"), action="store_true", default=opt$writePrbCSV,
                help="Boolean to write Probe Data (CSV) [default= %default]", metavar="boolean"),
    make_option(c("--writePrbRDS"), action="store_true", default=opt$writePrbRDS,
                help="Boolean to write Probe Data (RDS) [default= %default]", metavar="boolean"),
    make_option(c("--writeSSheet"), action="store_true", default=opt$writeSSheet,
                help="Boolean to write Auto-SampleSheet [default= %default]", metavar="boolean"),
    
    # Processing Parameters::
    make_option(c("--single"), action="store_true", default=opt$single, 
                help="Boolean variable to run a single sample on a single-core [default= %default]", metavar="boolean"),
    make_option(c("--parallel"), action="store_true", default=opt$parallel, 
                help="Boolean variable to run parallel on multi-core [default= %default]", metavar="boolean"),
    make_option(c("--cluster"), action="store_true", default=opt$cluster,
                help="Boolean variable to run jobs on cluster by chip [default= %default]", metavar="boolean"),

    #    
    # Function End::
    #
    
    make_option(c("--writeMIDMAN"), action="store_true", default=opt$writeMIDMAN,
                help="Boolean to write Sample Mid-Manifest [default= %default]", metavar="boolean"),
    make_option(c("--writeRDS"), action="store_true", default=opt$writeRDS,
                help="Boolean to write full-RDS data files [default= %default]", metavar="boolean"),
    make_option(c("--writeCSV"), action="store_true", default=opt$writeCSV,
                help="Boolean to write full-CSV data files [default= %default]", metavar="boolean"),
    
    make_option(c("--loadGRP"), action="store_true", default=opt$loadGRP,
                help="Boolean to load previous Master Grouped Manifest [default= %default]", metavar="boolean"),
    make_option(c("--loadRDS"), action="store_true", default=opt$loadRDS,
                help="Boolean to load previous RDS files [default= %default]", metavar="boolean"),
    
    # Save intermediate Manifest Files::
    make_option(c("--saveGRP"), action="store_true", default=opt$saveGRP,
                help="Boolean to write Master Grouped Manifest [default= %default]", metavar="boolean"),
  
    # Save intermediate SSET RDS files::  
    make_option(c("--saveRDS"), action="store_true", default=opt$saveRDS,
                help="Boolean to save current RDS files [default= %default]", metavar="boolean"),
    make_option(c("--overRDS"), action="store_true", default=opt$overRDS,
                help="Boolean to overwrite previous RDS files [default= %default]", metavar="boolean"),

    make_option(c("--useGRP"), action="store_true", default=opt$useGRP,
                help="Boolean to load CpG Grouping Variables R&D Purposes Only [default= %default]", metavar="boolean"),
    make_option(c("--retSSET"), action="store_true", default=opt$retSSET,
                help="Boolean to return SSET intermediate for debuging [default= %default]", metavar="boolean"),
    make_option(c("--retPRBS"), action="store_true", default=opt$retPRBS,
                help="Boolean to return PRBS intermediate for debuging [default= %default]", metavar="boolean"),
    
    # Verbosity::
    make_option(c("-v", "--verbosity"), type="double", default=opt$verbosity, 
                help="0-5 (5 is very verbosity) [default= %default]", metavar="double")
  )
  opt_parser = OptionParser(option_list=option_list)
  opt = parse_args(opt_parser)
}
opt$prgmPath <- prgmPath
if (!is.null(opt$funcScript) && !file.exists(opt$funcScript)) opt$funcScript <- file.path(opt$srcDir, 'R', opt$funcScript)
if (!is.null(opt$sesaScript) && !file.exists(opt$sesaScript)) opt$sesaScript <- file.path(opt$srcDir, 'R', opt$sesaScript)
if (!is.null(opt$idatScript) && !file.exists(opt$idatScript)) opt$idatScript <- file.path(opt$srcDir, 'R', opt$idatScript)

if (is.null(opt$Rscript)  || 
    is.null(opt$prgmTag)  || 
    is.null(opt$runMode)  ||
    is.null(opt$topDir)   ||
    is.null(opt$srcDir)   ||
    is.null(opt$idatsDir) ||
    is.null(opt$prgmPath) ||
    is.null(opt$platform) || 
    is.null(opt$build)    || 
    is.null(opt$method)   ||
    is.null(opt$poobMinPval) || is.null(opt$negsMinPval) ||
    is.null(opt$verbosity)) {
  print_help(opt_parser)
  optTib <- bind_rows(opt) %>% gather("Options", "Value")
  print(optTib %>% as.data.frame())
  stop("Null arguments!\n\n")
} else if (!dir.exists(opt$topDir)   ||
           !dir.exists(opt$srcDir)   ||
           !dir.exists(opt$idatsDir) ||
           !file.exists(opt$prgmPath) ||
           !file.exists(opt$Rscript)  ||
           !file.exists(opt$funcScript) ||
           !file.exists(opt$sesaScript) ||
           !file.exists(opt$idatScript)) {
  print_help(opt_parser)
  optTib <- bind_rows(opt) %>% gather("Options", "Value")
  print(optTib %>% as.data.frame())
  stop("Directory arguments do not exist!\n\n")
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                                  Begin::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

sysTime <- Sys.time()
cat(glue::glue("[{opt$prgmTag}]: Begin(time={sysTime})"),"\n\n",sep='')

if (opt$verbosity>2) {
  cat(glue::glue("[{opt$prgmTag}]: Options::"),"\n", sep='')
  opt.tib <- opt %>% bind_rows() %>% gather("Options", "Value")
  print(opt.tib)
}
source(opt$funcScript)
source(opt$sesaScript)
source(opt$idatScript)
opt$datDir <- file.path(opt$topDir, 'dat')
if (!dir.exists(opt$datDir)) {
  stop(glue::glue("[{opt$prgmTag}]: Critical ERROR: datDir={datDir} does not exist!"),"\n", sep='')
  q()
}

if (opt$fresh) {
  opt$writeMIDMAN <- FALSE
  opt$writeRDS    <- FALSE
  opt$writeCSV    <- FALSE
  
  opt$loadGRP     <- FALSE
  opt$loadRDS     <- FALSE
  
  # Save intermediate manifest csv files
  opt$saveMAN    <- FALSE
  opt$saveGRP    <- FALSE
  
  # Save intermediat RDS SSET files
  opt$saveRDS     <- FALSE
  opt$overRDS     <- FALSE
  
  ## CURRENT VARIABLES::
  opt$loadMAN     <- TRUE
  opt$loadIDATS   <- FALSE
  opt$loadSSETS   <- FALSE
  
  opt$writeIDATS  <- TRUE
  opt$writeSSETS  <- TRUE
  
  opt$writePrbCSV <- FALSE
  opt$writePrbRDS <- TRUE
  opt$writeSSheet <- TRUE

  rm.dir <- '/Users/bbarnes/Documents/Projects/workhorse/workspace/workhorse/testing'  
  if (dir.exists(rm.dir)) unlink(rm.dir, recursive=TRUE)
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                              Sample Sheet::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
ss.hum.tib <- NULL
if (!is.null(opt$sampleSheet)) {
  if (file.exists(opt$sampleSheet)) {
    ss.hum.tib <- readr::read_csv(ss.hum.tib)
  } else {
    cat("\n",glue::glue("[{opt$prgmTag}]: Warning: Human SampleSheet={opt$sampleSheet} does not exist! Skipping..."),"\n\n", sep='')
    stop("\n",glue::glue("[{opt$prgmTag}]: ERROR: Human SampleSheet={opt$sampleSheet} does not exist! Skipping..."),"\n\n")
  }
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                              Chip Prefixes::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
chipPrefixes <- sesame::searchIDATprefixes(opt$idatsDir, recursive=TRUE)
chipPrefixes.len <- length(chipPrefixes)
cat(glue::glue("[{opt$prgmTag}]: searchIDATprefixes={chipPrefixes.len}"),"\n",sep='')

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                                  Main::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

if (opt$cluster) {
  cat(glue::glue("[{opt$prgmTag}]: Launching Chips in Cluster Mode!"),"\n", sep='')
  
  isLinux <- FALSE
  lan.exe <- ''
  if (length(grep("^/illumina", opt$topDir)==1)) {
    lan.exe <- 'qsub -cwd -pe threaded 16 -l excl=true -N'
    isLinux <- TRUE
    # stop("\n",glue::glue("[{opt$prgmTag}]: ERROR: Cluster Mode is only allowed on Linux Cluster!!!"),"\n\n", sep='')
  }
  
  prefixTib <- tibble::tibble()
  for (prefixKey in names(chipPrefixes)) {
    prefixStr <- chipPrefixes[[prefixKey]]
    prefixDir <- basename(dirname(prefixStr))
    # cat("Dir=",prefixDir, ", Key=", prefixKey, ", PrefixStr=", prefixStr,"\n", sep='')
    loc.tib <- tibble::tibble( Barcode=prefixDir, Sentrix_Name=prefixKey, Prefix=prefixStr )
    prefixTib <- prefixTib %>% dplyr::bind_rows(loc.tib)
  }
  prefixTib <- prefixTib %>% split(.$Barcode)
  opt$shellDir <- file.path(opt$outDir, 'shells')
  if (!dir.exists(opt$shellDir)) dir.create(opt$shellDir, recursive=TRUE)
  
  for (chipName in names(prefixTib)) { # break }
    runShell <- file.path(opt$shellDir, paste0('run_',opt$prgmTag,'_',chipName,'.sh'))
    lanShell <- file.path(opt$shellDir, paste0('lan_',opt$prgmTag,'_',chipName,'.sh'))

    # Remove Cluster Option
    cmd.bool <- opt.tib %>% dplyr::filter(Value=='TRUE') %>% 
      dplyr::filter(Options!='cluster') %>%
      dplyr::select(Options) %>% 
      dplyr::mutate(Options=paste0('--',Options)) %>% .$Options %>% paste(collapse=" ")

    # Add ChipName to idat Directory
    cmd.strs <- opt.tib %>% dplyr::filter(Value!='TRUE' & Value!='FALSE' & Options!='prgmPath') %>%
      dplyr::mutate(Value=case_when(Options=='idatsDir' ~ paste0(Value,'/',chipName), TRUE ~ Value)) %>%
      dplyr::mutate(Options=paste0('--',Options,'=',Value)) %>%
      dplyr::select(Options) %>% .$Options %>% paste(collapse=" ")

    cmd.full <- paste(opt$Rscript, opt$prgmPath, cmd.bool, cmd.strs,"\n", sep=' ')
    readr::write_file(cmd.full, path=runShell)
    Sys.chmod(runShell, mode="0777")
    
    # Add cluster execute if avialbel (i.e. linux)
    cmd.lan <- runShell
    if (isLinux) 
      cmd.lan <- paste(lan.exe, paste('imWH',chipName, sep='_'), runShell,"\n", sep=' ')
    readr::write_file(cmd.lan, path=lanShell)
    Sys.chmod(lanShell, mode="0777")

    # Launch Script
    cat(glue::glue("[{opt$prgmTag}]: Launching chip={chipName}, shell={lanShell}!"),"\n", sep='')
    system(lanShell)
    
    if (opt$single) break
  }
  
} else {
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                              Preprocessing::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

  man.grp.rds <- file.path(opt$datDir, 'manifest', paste0(opt$platform,'.manifest.sesame.grouped.cpg-sorted.rds'))
  man.grp.csv <- file.path(opt$datDir, 'manifest', paste0(opt$platform,'.manifest.sesame.grouped.cpg-sorted.csv.gz'))

  man.cpg.rds <- file.path(opt$datDir, 'manifest', paste0(opt$platform,'.manifest.sesame.cpg-sorted.rds'))
  man.add.rds <- file.path(opt$datDir, 'manifest', paste0(opt$platform,'.manifest.sesame.add-sorted.rds'))
  
  man.cpg.csv <- file.path(opt$datDir, 'manifest', paste0(opt$platform,'.manifest.sesame.cpg-sorted.csv.gz'))
  man.add.csv <- file.path(opt$datDir, 'manifest', paste0(opt$platform,'.manifest.sesame.add-sorted.csv.gz'))
  
  man.cpg.tib <- NULL
  man.add.tib <- NULL
  man.grp.tib <- NULL
  
  # Load previous versions if both are defined
  if (opt$loadMAN &&
      file.exists(man.cpg.rds) && 
      file.exists(man.cpg.csv) &&
      file.exists(man.add.rds) &&
      file.exists(man.add.csv)) {
    
    cat(glue::glue("[{opt$prgmTag}]: Loading Sesame Manifest({opt$platform}) RDS(CPG)={man.cpg.rds}"),"\n", sep='')
    man.cpg.tib <- readr::read_rds(man.cpg.rds)
    cat(glue::glue("[{opt$prgmTag}]: Loading Sesame Manifest({opt$platform}) RDS(ADD)={man.add.rds}"),"\n", sep='')
    man.add.tib <- readr::read_rds(man.add.rds)
    cat(glue::glue("[{opt$prgmTag}]: Loaded CSV/RDS Manifests({opt$platform})."),"\n\n", sep='')
    
    opt$saveMAN <- FALSE
  } else {
    if (opt$platform=='EPIC' ||
        opt$platform=='HM450' ||
        opt$platform=='HM27') {
      stop("\n[Warning]: Not supported for standard manifests yet, platform =",opt$platform,"!!!\n\n")
    } else if (opt$platform=='EPIC-B4') {
      man.cpg.tib <- sesameManifest(platform='EPIC', build=opt$build, verbose=opt$verbosity)
      man.add.tib <- man.cpg.tib %>% dplyr::arrange(U)
    } else if (opt$platform=='NZT') {
      # epic.ctl.csv <- file.path(opt$topDir, 'dat/ctl/EPIC-B0.controls.sesame.csv.gz')
      # epic.ctl_tib <- suppressMessages(suppressWarnings(readr::read_csv(epic.ctl.csv) ))
      
      man.cpg.csv <- file.path(opt$topDir, 'dat/manifest/NZT.manifest.ann.csv.gz')
      man.cpg.tib <- suppressMessages(suppressWarnings(readr::read_csv(man.cpg.csv) )) %>%
        dplyr::arrange(Probe_ID)
      man.add.tib <- man.cpg.tib %>% dplyr::arrange(U)
    } else {
      stop("\n\nERROR UNSUPPORTED PLATFORM 3!!!\n\n")
    }
  }

  # Seperate Case for Loading Grouped Manifest::
  if (opt$loadGRP &&
      file.exists(man.grp.rds)) {
    cat(glue::glue("[{opt$prgmTag}]: Loading Sesame Grouped Manifest({opt$platform}) RDS(CPG)={man.grp.rds}"),"\n", sep='')
    man.grp.tib <- readr::read_rds(man.grp.rds)
  } else {
    if (opt$useGRP) man.grp.tib <- addManifestGroups(platform='EPIC', build=opt$build, manCPG=man.cpg.tib, verbose=opt$verbosity)
  }

  if (opt$saveMAN) {
    cat(glue::glue("[{opt$prgmTag}]: Writing CSV/RDS Manifests datDir={opt$datDir}"),"\n", sep='')
    readr::write_csv(man.cpg.tib, man.cpg.csv)
    readr::write_rds(man.cpg.tib, man.cpg.rds, compress='gz')
    readr::write_csv(man.add.tib, man.add.csv)
    readr::write_rds(man.add.tib, man.add.rds, compress='gz')
    cat(glue::glue("[{opt$prgmTag}]: Done writing CSV/RDS Manifests"),"\n\n", sep='')
  }
  if (opt$saveGRP) {
    if (!is.null(man.grp.tib)) {
      cat(glue::glue("[{opt$prgmTag}]: Writing CSV/RDS Grouped Manifests datDir={opt$datDir}"),"\n", sep='')
      readr::write_csv(man.grp.tib, man.grp.csv)
      readr::write_rds(man.grp.tib, man.grp.rds, compress='gz')
      cat(glue::glue("[{opt$prgmTag}]: Done writing CSV/RDS Grouped Manifests"),"\n\n", sep='')
    }
  }
  
  # Load data for auto detect
  #
  # TBD: Replace ref.sam.tib with::
  #  - can.tib[[decoder]] => Canonical Reference
  #  - ref.tib[[decoder]] => Sample Reference
  can.tib <- NULL
  ref.tib <- NULL
  # ref.sam.tib <- NULL
  if (opt$autoDetect) {
    Pools <- c('BETA', 'DELTA')
    for (pool in Pools) {
      can.rds <- file.path(opt$topDir, 'dat/ref', paste0('Canonical_',pool,'_Betas-Mvals.rds'))
      ref.rds <- file.path(opt$topDir, 'dat/ref', paste0('Reference_',pool,'_Betas-Mvals.rds'))
      
      if (opt$verbosity>1) cat(glue::glue("[{opt$prgmTag}]: Loading Canonical({pool}) RDS={can.rds}"),"\n", sep='')
      can.tib[[pool]] <- readr::read_rds(can.rds)
      
      if (FALSE) {
        if (opt$verbosity>1) cat(glue::glue("[{opt$prgmTag}]: Loading Reference({pool}) RDS={ref.rds}"),"\n", sep='')
        ref.tib[[pool]] <- readr::read_rds(ref.rds)
      }
      if (opt$verbosity>1) cat(glue::glue("[{opt$prgmTag}]: Loaded both Canonical/Reference({pool}) Samples for Auto-Detect!"),"\n\n", sep='')
    }
  }
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                                Main::Chips
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

  ssheets <- NULL
  if (opt$paralle) {
    cat(glue::glue("[{opt$prgmName}]: Will use multi-core parallel processing..."),"\n", sep='')
    ssheets <- foreach (prefix=names(chipPrefixes), .inorder=T, .final = function(x) setNames(x, names(chipPrefixes))) %dopar% {
      ss.ret <- singleSampleWorkflow(prefix=chipPrefixes[[prefix]],
                                     opt=opt,
                                     man.cpg = man.cpg.tib,
                                     man.add = man.add.tib,
                                     man.ses = NULL,
                                     can.tib = can.tib,
                                     ref.tib = ref.tib,
                                     ctls_tib = NULL,
                                     vt=1)
      ss.ret
    }
  } else {
    cat(glue::glue("[{opt$prgmName}]: Will use single-core linear processing..."),"\n", sep='')
    for (prefix in names(chipPrefixes)) {
      if (!is.null(tarSample) && prefix!=tarSample) next
      
      # opt$retSSET <- TRUE
      # opt$fullData <- TRUE
      # opt$writeMIDMAN <- FALSE
      # opt$saveRDS  <- FALSE
      ss.ret <- singleSampleWorkflow(prefix=chipPrefixes[[prefix]],
                                     opt=opt,
                                     man.cpg = man.cpg.tib,
                                     man.add = man.add.tib,
                                     man.ses = NULL,
                                     can.tib = can.tib,
                                     ref.tib = ref.tib,
                                     ctls_tib = NULL,
                                     vt=1)

      # if (opt$fullData) {
      #   ss.print <- ss.ret$sam_ss_tib %>% gather() %>% as.data.frame() %>% head(n=50)
      #   print(ss.print)
      # }
      
      ssheets[[prefix]] <- ss.ret
      if (opt$single) break
    }
  }
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                                Finished::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
sysTime <- Sys.time()
cat(glue::glue("[{opt$prgmTag}]: Finished(time={sysTime})"),"\n\n",sep='')

# End of file
