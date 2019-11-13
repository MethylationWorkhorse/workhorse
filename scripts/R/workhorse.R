
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

opt$writeSSheet <- TRUE
opt$writeMIDMAN <- TRUE

opt$writeMAN    <- TRUE
opt$writeGRP    <- TRUE
opt$writeRDS    <- TRUE
opt$writeCSV    <- TRUE

opt$loadMAN     <- FALSE
opt$loadGRP     <- FALSE
opt$loadRDS     <- FALSE
opt$saveRDS     <- TRUE
opt$overRDS     <- TRUE

opt$retSSET     <- FALSE
opt$retPRBS     <- FALSE

opt$Rscript <- '/usr/local/bin/Rscript'
opt$topDir  <- '/Users/bbarnes/Documents/Projects/workhorse'
opt$srcDir  <- '/Users/bbarnes/Documents/Projects/workhorse/scripts'
opt$outDir  <- NULL
opt$idatsDir <- NULL
opt$funcScript <- 'workhorse_functions.R'
opt$sampleSheet <- NULL

opt$platform  <- 'EPIC-B4'
opt$build     <- 'hg19'
opt$method    <- 'both'

opt$poobMinPval <- 0.2
opt$negsMinPval <- 0.02

opt$verbosity <- 3

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
  prgmPath <- file.path(opt$srcDir, 'R', paste0(opt$prgmTag,'.R'))
  stopifnot(file.exists(prgmPath))
  
  opt$runMode  <- args.dat[1]
  opt$topDir   <- dirname(opt$srcDir)
  opt$outDir   <- file.path(opt$topDir, 'workspace', opt$prgmTag)
  opt$idatsDir <- file.path(opt$topDir, 'idats')
  
  # opt$platform  <- 'NZT'
  opt$platform  <- 'EPIC-B4'
  if (opt$platform=='EPIC-B4') {
    opt$idatsDir  <- file.path(opt$topDir, 'idats', 'DeltaBetaCore')
    
    # opt$idatsDir  <- file.path(opt$topDir, 'idats', 'DeltaBetaCore','202761400009')
    tarSample <- '202761400009_R07C01'  # M
    tarSample <- '202761400009_R06C01'  # H
    tarSample <- '202761400009_R05C01'  # U
  } else if (opt$platform=='NZT') {
    opt$idatsDir  <- file.path(opt$topDir, 'idats', 'NZT')
  } else {
    stop("\nUNSUPPORTED PLATFORM!!!\n\n")
  }
  
  opt$autoDetect <- TRUE
  opt$single     <- FALSE
  opt$parallel   <- FALSE
  opt$cluster    <- FALSE

  opt$writeSSheet <- TRUE
  opt$writeMIDMAN <- FALSE
  
  opt$writeMAN    <- FALSE
  opt$writeRDS    <- TRUE
  opt$writeCSV    <- FALSE
  
  opt$loadMAN     <- TRUE
  opt$loadRDS     <- TRUE
  opt$saveRDS     <- TRUE
  opt$overRDS     <- TRUE
  
  opt$retSSET     <- FALSE
  opt$retPRBS     <- FALSE
  
} else {
  prgmPath <- substring(args.dat[grep("--file=", args.dat)], 8)
  
  opt$runMode <- 'CommandLine'
  opt$prgmTag <- sub('\\.R$', '', basename(substring(args.dat[grep("--file=", args.dat)], 8)))
  opt$srcDir  <- dirname(normalizePath(dirname(substring(args.dat[grep("--file=", args.dat)], 8)) ))
  opt$topDir  <- dirname(opt$srcDir)
  
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
                help="name of functions script to soruce [default= %default]", metavar="character"),
    
    # Platform/Method Parameters::
    make_option(c("--platform"), type="character", default=opt$platform, 
                help="platform to use for analysis [default= %default]", metavar="character"),
    make_option(c("--build"), type="character", default=opt$build, 
                help="genome build to use for analysis [default= %default]", metavar="character"),
    make_option(c("--method"), type="character", default=opt$method, 
                help="method to use for analysis [default= %default]", metavar="character"),
    
    # Detection Parameters::
    make_option(c("--poobMinPval"), type="double", default=opt$poobMinPval, 
                help="detection p-value filter for poob controls [default= %default]", metavar="double"),
    make_option(c("--negsMinPval"), type="double", default=opt$negsMinPval, 
                help="detection p-value filter for negative controls [default= %default]", metavar="double"),
    
    # Boolean Variables::
    make_option(c("--autoDetect"), action="store_true", default=opt$autoDetect,
                help="Boolean variable to attempt auto sample detect [default= %default]", metavar="boolean"),
    
    make_option(c("--writeSSheet"), action="store_true", default=opt$writeSSheet,
                help="Boolean to write Auto-SampleSheet [default= %default]", metavar="boolean"),
    make_option(c("--writeMIDMAN"), action="store_true", default=opt$writeMIDMAN,
                help="Boolean to write Sample Mid-Manifest [default= %default]", metavar="boolean"),
    make_option(c("--writeMAN"), action="store_true", default=opt$writeMAN,
                help="Boolean to write Master Manifest [default= %default]", metavar="boolean"),
    make_option(c("--writeGRP"), action="store_true", default=opt$writeGRP,
                help="Boolean to write Master Grouped Manifest [default= %default]", metavar="boolean"),
    make_option(c("--writeRDS"), action="store_true", default=opt$writeRDS,
                help="Boolean to write full-RDS data files [default= %default]", metavar="boolean"),
    make_option(c("--writeCSV"), action="store_true", default=opt$writeCSV,
                help="Boolean to write full-CSV data files [default= %default]", metavar="boolean"),
    
    make_option(c("--loadMAN"), action="store_true", default=opt$loadMAN,
                help="Boolean to load previous Manifest files [default= %default]", metavar="boolean"),
    make_option(c("--loadGRP"), action="store_true", default=opt$loadGRP,
                help="Boolean to load previous Master Grouped Manifest [default= %default]", metavar="boolean"),
    make_option(c("--loadRDS"), action="store_true", default=opt$loadRDS,
                help="Boolean to load previous RDS files [default= %default]", metavar="boolean"),
    make_option(c("--saveRDS"), action="store_true", default=opt$saveRDS,
                help="Boolean to save current RDS files [default= %default]", metavar="boolean"),
    make_option(c("--overRDS"), action="store_true", default=opt$overRDS,
                help="Boolean to overwrite previous RDS files [default= %default]", metavar="boolean"),
    make_option(c("--retSSET"), action="store_true", default=opt$retSSET,
                help="Boolean to return SSET intermediate for debuging [default= %default]", metavar="boolean"),
    make_option(c("--retPRBS"), action="store_true", default=opt$retPRBS,
                help="Boolean to return PRBS intermediate for debuging [default= %default]", metavar="boolean"),
    
    # Processing Parameters::
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
opt$prgmPath <- prgmPath
if (!is.null(opt$funcScript) && !file.exists(opt$funcScript))
  opt$funcScript <- file.path(opt$srcDir, 'R', opt$funcScript)

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
           !file.exists(opt$funcScript) ) {
  print_help(opt_parser)
  optTib <- bind_rows(opt) %>% gather("Options", "Value")
  print(optTib %>% as.data.frame())
  stop("Directory arguments do not exist!\n\n")
}

if (opt$verbosity>2) {
  cat(glue::glue("[{opt$prgmTag}]: Options::"),"\n", sep='')
  opt.tib <- opt %>% bind_rows() %>% gather("Options", "Value")
  print(opt.tib)
}
source(opt$funcScript)
opt$topDir <- dirname(opt$srcDir)
opt$datDir <- file.path(opt$topDir, 'dat')
if (!dir.exists(opt$datDir)) dir.create(opt$topDir, recursive=TRUE)

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
  if (opt$loadMAN==TRUE &&
      file.exists(man.cpg.rds) && 
      file.exists(man.cpg.csv) &&
      file.exists(man.add.rds) &&
      file.exists(man.add.csv)) {
    
    cat(glue::glue("[{opt$prgmTag}]: Found all four CSV/RDS Manifests(opt$platform)."),"\n", sep='')
    
    if (file.exists(man.grp.rds)) {
      cat(glue::glue("[{opt$prgmTag}]: Loading Sesame Grouped Manifest({opt$platform}) RDS(CPG)={man.grp.rds}"),"\n", sep='')
      man.grp.tib <- readr::read_rds(man.grp.rds)
    }
    
    cat(glue::glue("[{opt$prgmTag}]: Loading Sesame Manifest({opt$platform}) RDS(CPG)={man.cpg.rds}"),"\n", sep='')
    man.cpg.tib <- readr::read_rds(man.cpg.rds)
    cat(glue::glue("[{opt$prgmTag}]: Loading Sesame Manifest({opt$platform}) RDS(ADD)={man.add.rds}"),"\n", sep='')
    man.add.tib <- readr::read_rds(man.add.rds)
    cat(glue::glue("[{opt$prgmTag}]: Loaded CSV/RDS Manifests({opt$platform})."),"\n\n", sep='')
    
    opt$writeMAN <- FALSE
  } else {
    if (opt$platform=='EPIC' ||
        opt$platform=='HM450' ||
        opt$platform=='HM27') {
      stop("\n[Warning]: Not supported for standard manifests yet, platform =",opt$platform,"!!!\n\n")
    } else if (opt$platform=='EPIC-B4') {
      man.cpg.tib <- sesameManifest(platform='EPIC', build=opt$build, verbose=opt$verbosity)
      man.grp.tib <- addManifestGroups(platform='EPIC', build=opt$build, manCPG=man.cpg.tib, verbose=opt$verbosity)
      man.add.tib <- man.cpg.tib %>% dplyr::arrange(U)
    } else if (opt$platform=='NZT') {
      # epic.ctl.csv <- file.path(opt$topDir, 'dat/ctl/EPIC-B0.controls.sesame.csv.gz')
      # epic.ctl.tib <- suppressMessages(suppressWarnings(readr::read_csv(epic.ctl.csv) ))
      
      man.cpg.csv <- file.path(opt$topDir, 'dat/manifest/NZT.manifest.ann.csv.gz')
      man.cpg.tib <- suppressMessages(suppressWarnings(readr::read_csv(man.cpg.csv) )) %>%
        dplyr::arrange(Probe_ID)
      man.add.tib <- man.cpg.tib %>% dplyr::arrange(U)
    } else {
      stop("\n\nERROR UNSUPPORTED PLATFORM 3!!!\n\n")
    }
  }
  
  if (opt$writeMAN) {
    cat(glue::glue("[{opt$prgmTag}]: Writing CSV/RDS Manifests datDir={opt$datDir}"),"\n", sep='')
    readr::write_csv(man.cpg.tib, man.cpg.csv)
    readr::write_rds(man.cpg.tib, man.cpg.rds, compress='gz')
    readr::write_csv(man.add.tib, man.add.csv)
    readr::write_rds(man.add.tib, man.add.rds, compress='gz')
    if (!is.null(man.grp.tib)) {
      cat(glue::glue("[{opt$prgmTag}]: Writing CSV/RDS Grouped Manifests datDir={opt$datDir}"),"\n", sep='')
      readr::write_csv(man.grp.tib, man.grp.csv)
      readr::write_rds(man.grp.tib, man.grp.rds, compress='gz')
    }
    cat(glue::glue("[{opt$prgmTag}]: Done writing CSV/RDS Manifests"),"\n\n", sep='')
  }
  
  # Load data for auto detect
  ref.sam.tib <- NULL
  if (opt$autoDetect) {
    ref.rds <- file.path(opt$topDir, 'dat/ref/Reference_betasMvals_Pval0.02.rds')
    cat(glue::glue("[{opt$prgmTag}]: Loading Reference RDS={ref.rds}"),"\n", sep='')
    ref.sam.tib <- readr::read_rds(ref.rds)
    cat(glue::glue("[{opt$prgmTag}]: Loaded Reference Samples"),"\n\n", sep='')
  }
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                                  Main::
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
                                     ref.tib = ref.sam.tib,
                                     ctl.tib = NULL,
                                     vt=1)
      ss.ret
    }
  } else {
    cat(glue::glue("[{opt$prgmName}]: Will use single-core linear processing..."),"\n", sep='')
    for (prefix in names(chipPrefixes)) {
      # if (!is.null(tarSample) && prefix!=tarSample) next
      
      ss.ret <- singleSampleWorkflow(prefix=chipPrefixes[[prefix]],
                                     opt=opt,
                                     man.cpg = man.cpg.tib,
                                     man.add = man.add.tib,
                                     man.ses = NULL,
                                     ref.tib = ref.sam.tib,
                                     ctl.tib = NULL,
                                     vt=1)
      ssheets[[prefix]] <- ss.ret
      if (opt$single) break
    }
  }
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                                Finished::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
ctime <- Sys.time()
cat(glue::glue("[{opt$prgmTag}]: Finished(time={ctime})"),"\n\n",sep='')

# End of file
