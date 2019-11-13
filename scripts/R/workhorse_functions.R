
suppressWarnings(suppressPackageStartupMessages(require("tidyverse")) )
suppressWarnings(suppressPackageStartupMessages(require("stringr")) )
suppressWarnings(suppressPackageStartupMessages(require("glue")) )
suppressWarnings(suppressPackageStartupMessages(require("scales")) )
suppressWarnings(suppressPackageStartupMessages(require("matrixStats")) )

COM <- ","
TAB <- "\t"
RET <- "\n"

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                            Workflow Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

singleSampleWorkflow = function(prefix, 
                                opt,
                                man.cpg, 
                                man.add,
                                man.ses=NULL, 
                                ref.tib=NULL, 
                                ctl.tib=NULL, 
                                vt=1) {
  funcTag <- 'singleSampleWorkflow'
  stopifnot(length(opt$platform)>0)
  verbose <- opt$verbosity
  
  if (is.null(prefix) || length(prefix)==0) 
    stop("\n", glue::glue("[{funcTag}]: ERROR prefix is NULL"),"\n", sep='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{TAB} Starting prefix={prefix}."),"\n\n", sep='')
  
  # Sanity Checking::
  stopifnot(!is.null(opt$platform))
  
  # Get Raw Idat Data::
  idat.list <- prefixToTib(prefix, verbose=verbose)
  # return(idat.list)
  
  # Build Local Output Directory::
  loc.outDir <- file.path(opt$outDir, idat.list$Barcode, opt$platform)
  if (!dir.exists(loc.outDir)) dir.create(loc.outDir, recursive=TRUE)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{TAB} Local OutDir={loc.outDir}"),"\n", sep='')
  
  # Build Output Files::
  opt$outName <- paste(idat.list$Sentrix_Name,opt$platform,opt$method, sep='_')
  opt$prbs1_CSV  <- file.path(loc.outDir, paste0(opt$outName, '.ProbesI.tib.csv.gz'))
  opt$prbs2_CSV  <- file.path(loc.outDir, paste0(opt$outName, '.ProbesII.tib.csv.gz'))
  opt$prbs_RDS   <- file.path(loc.outDir, paste0(opt$outName, '.Probes.tib.rds'))
  opt$ss_CSV     <- file.path(loc.outDir, paste0(opt$outName, '.SampleSheet.csv.gz'))
  opt$man_CSV    <- file.path(loc.outDir, paste0(opt$outName, '.SampleSesameManifest.csv.gz'))
  opt$sset_RDS   <- file.path(loc.outDir, paste0(opt$outName, '.sset.rds'))
  
  # Get Sesame Sample Manifest::
  # if (is.null(man.ses) || !opt$cleanManifest)
  #   man.ses <- getManifest(man.add, idat.list$tib, verbose=verbose)
  #   man.ses <- getManifest(man.cpg, idat.list$tib, verbose=verbose)
  man.ses <- getManifest(man.add, idat.list$tib, verbose=verbose)
  stopifnot(length(man.ses)>0)
  if (opt$writeMIDMAN) {
    cat(glue::glue("[{funcTag}]:{TAB} Writing Sample-Mid-Manifest={opt$man_CSV}"),"\n", sep='')
    readr::write_csv(man.ses, opt$man_CSV)
  }
  
  # Build SSET::
  # ssets <- fetchSSet(prefix, platform='custom', manifest=man.cpg, method=opt$method, desplat=opt$platform,
  # ssets <- fetchSSet(prefix, platform='custom', manifest=man.add, method=opt$method, desplat=opt$platform,
  #                    ssetRDS=opt$sset_RDS, loadRDS=opt$loadRDS, saveRDS=opt$saveRDS, overRDS=opt$overRDS,
  #                    verbose=verbose)
  ssets <- fetchSSet(prefix, platform='custom', manifest=man.add, method=opt$method, desplat=opt$platform,
                     ssetRDS=opt$sset_RDS, loadRDS=opt$loadRDS, saveRDS=opt$saveRDS, overRDS=opt$overRDS,
                     verbose=verbose)
  
  if (opt$retSSET) return(ssets)
  ssets.tib   <- NULL
  ssets.merge <- NULL
  
  prbs <- ssetsToTibs(ssets, verbose=verbose)
  if (verbose>=vt+5) print(prbs)
  if (opt$retPRBS) return(prbs)
  
  if (is.null(ctl.tib))
    ctl.tib <- man.add %>% dplyr::filter(stringr::str_starts(Probe_ID, 'ctl')) %>% dplyr::select(Probe_ID, Probe_Type)
  if (verbose>=vt) cat("\t", glue::glue("[{funcTag}]: Extracted Control Probes"),"\n", sep='')
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                          Auto Sample Detection::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  auto.tib <- tibble(sam.ref=5,
                     Auto_Sample_DB_Key='Unknown',
                     Auto_Sample_DB_Val=NA,
                     Auto_Sample_R2_Key='Unknown',
                     Auto_Sample_R2_Val=NA)
  
  if (opt$autoDetect && !is.null(ref.tib)) {
    if (verbose>=vt) cat("\t", glue::glue("[{funcTag}]: Starting Auto Sample Detection..."),"\n", sep='')
    # New method that selects common columns before binding row::
    prbsI  <- prbs$I  %>% dplyr::select(Probe_ID, Probe_Type, Beta, Mval, NegsDetP, PoobDetP)
    prbsII <- prbs$II %>% dplyr::select(Probe_ID, Probe_Type, Beta, Mval, NegsDetP, PoobDetP)
    sam <- dplyr::bind_rows(prbsI, prbsII) %>% dplyr::arrange(Probe_ID)
    
    # Specail Case for NZT Pool::
    #
    if (opt$platform=='NZT')
      sam <- sam %>% dplyr::filter(str_detect(Probe_ID, '_C_I')) %>% 
      dplyr::mutate(Probe_ID=str_remove(Probe_ID, '_[FR].*$'))
    
    auto.tib <- predictAutoSample(sam=sam, ref=ref.tib, minPval=opt$negsMinPval, minDelta=0.2, verbose=verbose)
    if (verbose>=vt) cat("\t", glue::glue("[{funcTag}]: Detected Auto Samples(DB)={auto.tib$Auto_Sample_DB_Key}..."),"\n", sep='')
    
    if (verbose>=vt+4) print(auto.tib)
    # return(auto.tib)
  }
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                        Output Analytical Tibble::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  # Write merged sample::
  if (opt$writeRDS) {
    if (verbose>=vt) cat("\t", glue::glue("[{funcTag}]: Writing Probe RDS={opt$prbs_RDS}"),"\n", sep='')
    readr::write_rds(prbs, opt$prbs_RDS, compress="gz")
  }
  if (opt$writeCSV) {
    if (verbose>=vt) cat("\t", glue::glue("[{funcTag}]: Writing Probe(I) CSV={opt$prbs1_CSV}"),"\n", sep='')
    readr::write_csv(prbs$I,  opt$prbs1_CSV)
    if (verbose>=vt) cat("\t", glue::glue("[{funcTag}]: Writing Probe(II) CSV={opt$prbs2_CSV}"),"\n", sep='')
    readr::write_csv(prbs$II, opt$prbs2_CSV)
  }
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                   Auto Sample Sheet/Summary Generation::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  if (verbose>=vt) cat("\t", glue::glue("[{funcTag}]: Building Sample Sheet/Summary Stats..."),"\n", sep='')
  prbs.sum <- summarizePrbs(prbs=prbs, ctl=ctl.tib, 
                            negsMin=opt$negsMinPval, poobMin=opt$poobMinPval, verbose=verbose)
  if (verbose>=vt+4) print(prbs.sum)
  
  if (verbose>=vt) cat("\t", glue::glue("[{funcTag}]: {TAB}Building core.sum..."),"\n", sep='')
  core.lst <- NULL
  core.lst$analysisPlatform <- opt$platform
  core.lst$Sentrix_Name <- idat.list$Sentrix_Name
  core.lst$Barcode      <- idat.list$Barcode
  core.lst$Poscode      <- idat.list$Poscode
  core.lst$Row          <- idat.list$row
  core.lst$Column       <- idat.list$col
  core.lst$ChipType     <- idat.list$ChipType
  core.lst$ChipFormat   <- idat.list$ChipFormat
  core.sum <- dplyr::bind_rows(core.lst)
  if (verbose>=vt+4) print(core.sum)
  
  if (verbose>=vt) cat("\t", glue::glue("[{funcTag}]: {TAB}Merging Summaries..."),"\n", sep='')
  sam.ss <- core.sum %>%
    dplyr::bind_cols(auto.tib) %>%
    # dplyr::bind_cols(phen) %>%
    # dplyr::bind_cols(ages) %>%
    dplyr::bind_cols(prbs.sum) %>%
    dplyr::mutate_if(is.numeric, list(round), 6)
  if (verbose>=vt+4) print(sam.ss)
  
  if (opt$writeSSheet) {
    if (verbose>=vt)
      cat("\t", glue::glue("[{funcTag}]: Writing Auto Sample Sheet CSV={opt$outCsv}"),"\n", sep='')
    readr::write_csv(sam.ss, opt$ss_CSV)
  }
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{TAB} Finished."),"\n\n", sep='')
  
  sam.ss
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                         Auto Sample Functions ::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

firstVsRestDelta = function(mat, minDelta=0.2, verbose=0, vt=5) {
  funcTag <- 'firstVsRestDelta'
  if (verbose>=vt) cat("\t",glue::glue("[{funcTag}]: {TAB}Starting..."),"\n", sep='')
  
  ncols  <- base::ncol(mat)
  nrows  <- base::nrow(mat)
  cnames <- colnames(mat)
  
  stopifnot(nrows>0)
  
  max_cnt <- 0
  max_idx <- 0
  max_per <- 0
  max_key <- 'Unknown'
  for (i in c(2:ncols)) {
    cur_cnt <- length(which(rowDiffs(mat[,c(1,i)]) <= minDelta))
    if (verbose>=vt) cat("\t",glue::glue("[{funcTag}]: i={i}, matchCnt={cur_cnt}"),"\n", sep='')
    if (cur_cnt > max_cnt) {
      max_cnt <- cur_cnt
      max_idx <- i
      if (verbose>=vt) cat("\t",glue::glue("[{funcTag}]: {TAB}{TAB} maxCnt={max_cnt}, maxIdx={max_idx}"),"\n", sep='')
    }
  }
  if (max_idx != 0) max_key <- cnames[max_idx]
  max_per <- round(100*max_cnt/nrows, 3)
  
  mtib <- tibble("Auto_Sample_DB_Key" = max_key,
                 "Auto_Sample_DB_Val" = max_per)
  if (verbose>=vt) cat("\t",glue::glue("[{funcTag}]: {TAB}Finished."),"\n\n", sep='')
  
  mtib
}

predictAutoSample = function(sam, ref, minPval=0.02, minDelta=0.2, verbose=0, vt=5) {
  funcTag <- 'predictAutoSample'
  if (verbose>=vt) cat("\t",glue::glue("[{funcTag}]: {TAB}Starting..."),"\n", sep='')
  
  # cmpTypes <- c('Beta_base', 'Beta_swap', 'Mval_base', 'Mval_swap')
  # probeTypes <- c('cg', 'rs', 'ch')
  
  # cmpTypes <- c('Beta_base')
  cmpTypes <- c('Beta')
  probeTypes <- c('cg')
  
  sam.split <- sam %>% split(.$Probe_Type)
  
  # Delta Mval(swap)
  # R2 Mval(swap)
  cur.merge.dat <- NULL
  for (pt in probeTypes) {
    if (length(sam.split[[pt]])==0 || length(ref[[pt]])==0) {
      if (opt$verbosity>=vt) cat("\t", glue::glue("[{funcTag}]: Skipping Auto Detect for probeType={pt}..."),"\n", sep='')
      next
    }
    
    # Remove poor performing detection p-values using the min(NEGs between base/swap)
    #  sam.split[[pt]] %>% dplyr::select(Probe_ID,starts_with("Beta"),starts_with("Mval"))
    #  sam.split[[pt]] %>% dplyr::select(contains("NegsDetP_base")) 
    
    # Remove Poorly Detected Probes based on Negs DetP for now::
    #  Below is the old method with all options::
    #  - baseFailNegs_idx <- which(sam.split[[pt]]$NegsDetP_base > minPval)
    baseFailNegs_idx <- which(sam.split[[pt]]$NegsDetP > minPval)
    
    atib <- NULL
    dtib <- NULL
    rtib <- NULL
    for (cmpType in cmpTypes) {
      sdat <- sam.split[[pt]] %>% dplyr::select(Probe_ID,starts_with(cmpType))
      sdat[[cmpType]][baseFailNegs_idx] <- NA
      
      # Get Betas For Reference::
      rdat <- ref[[pt]] %>% select(Probe_ID, ends_with(cmpType))
      
      sr_mat <- sdat %>% dplyr::left_join(rdat, by="Probe_ID") %>% 
        dplyr::select(-Probe_ID) %>% as.matrix()
      
      # Get Delta Beta::
      dtib <- firstVsRestDelta(sr_mat, verbose=verbose, vt=vt)
      
      # Get R-Squared::
      sr_cor <- NULL
      sr_cor = tryCatch({
        cor(sr_mat, method="pearson", use="pairwise.complete.obs")
      }, warning = function(w) {
        'warning-Cor'
      }, error = function(e) {
        'error-Cor'
      }, finally = {
        'cleanup-Cor'
      })
      # sr_cor <- cor(sr_mat, method="pearson", use="pairwise.complete.obs")
      
      mkey <- NA
      mval <- NA
      
      if (!is.null(sr_cor)) {
        if (verbose>=vt) print(sr_cor[1,])
        ncols <- length(sr_cor[1,])
        mval <- max(sr_cor[1,2:ncols])
        midx <- which(sr_cor[1,2:ncols]==mval)
        mkey <- names(sr_cor[1,2:ncols])[midx]
        mval <- round(max(sr_cor[1,2:ncols]),6)
        
        if (verbose>=vt) {
          cat("\tmval=",mval,"\n", sep='')
          cat("\tmidx=",midx,"\n", sep='')
          cat("\tmkey=",mkey,"\n", sep='')
        }
      }
      rtib <- tibble("Auto_Sample_R2_Key" = mkey,
                     "Auto_Sample_R2_Val" = mval)
      
      atib <- dplyr::bind_cols(dtib, rtib)
    }
  }
  if (verbose>=vt) cat("\t",glue::glue("[{funcTag}]: {TAB}Finished."),"\n\n", sep='')
  
  atib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                         Sample Summary Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

summarizePrbs = function(prbs, ctl, negsMin=0.02, poobMin=0.2, verbose=0, vt=3) {
  funcTag <- 'summarizePrbs'
  
  if (verbose>=vt) cat("\t", glue::glue("[{funcTag}]: Starting..."),"\n", sep='')
  
  # prbsI   <- prbs$I  %>% dplyr::filter(Probe_Type!='ctl') %>% dplyr::group_by(Probe_Type,Design_Type)
  prbsI   <- prbs$I  %>% dplyr::filter(Probe_Type!='ctl') %>% dplyr::group_by(Probe_Type)
  if (verbose>=vt+2) print(prbsI)
  prbsII  <- prbs$II %>% dplyr::filter(Probe_Type!='ctl') %>% dplyr::group_by(Probe_Type)
  if (verbose>=vt+2) print(prbsII)
  prbsCtl <- prbs$II %>% dplyr::filter(Probe_Type=='ctl') %>% dplyr::select(-Probe_Type) %>%
    dplyr::left_join(ctl, by="Probe_ID") %>%
    dplyr::select(Probe_ID, Probe_Type, everything()) %>% dplyr::group_by(Probe_Type)
  if (verbose>=vt+2) print(prbsCtl)
  
  # Get Overall Detection P-value for CGs only::
  cgs1 <- prbsI  %>% dplyr::ungroup() %>% dplyr::filter(Probe_Type=='cg') %>% 
    dplyr::select(NegsDetP, PoobDetP, Swapped)
  cgs2 <- prbsII %>% dplyr::ungroup() %>% dplyr::filter(Probe_Type=='cg') %>% 
    dplyr::select(NegsDetP, PoobDetP) %>% dplyr::mutate(Swapped=NA)
  sumCG <- cgs1 %>% dplyr::bind_rows(cgs2) %>%
    dplyr::summarise(Total_Count_CG=n(), 
                     PassDetpNegs_Count_CG=count(NegsDetP<negsMin), 
                     PassDetpNegs_Percent_CG=round(100*PassDetpNegs_Count_CG/Total_Count_CG,3),
                     PassDetpPoob_Count_CG=count(PoobDetP<poobMin), 
                     PassDetpPoob_Percent_CG=round(100*PassDetpPoob_Count_CG/Total_Count_CG,3),
                     SwappedCount_CG=count(!is.na(Swapped)) )
  # return(sumCG)
  
  # Infinium I::
  sumA <- prbsI %>%
    dplyr::summarise(Total_Count=n(), 
                     PassDetpNegs_Count=count(NegsDetP<negsMin), 
                     PassDetpNegs_Percent=round(100*PassDetpNegs_Count/Total_Count,3),
                     PassDetpPoob_Count=count(PoobDetP<poobMin), 
                     PassDetpPoob_Percent=round(100*PassDetpPoob_Count/Total_Count,3),
                     SwappedCount=count(!is.na(Swapped)))
  
  sumB <- prbsI %>%
    dplyr::summarise_at(vars(Beta, Mval, GM, GU, RM, RU), list(min=min, 
                                                               median=median, 
                                                               max=max,
                                                               mean=mean,
                                                               sd=sd), na.rm=TRUE)
  
  sumI <- sumA %>% dplyr::full_join(sumB, by="Probe_Type") %>%
    tidyr::gather(Stat, Value, -Probe_Type) %>% 
    tidyr::unite(Probe_Type, Probe_Type, Stat, sep='_') %>% 
    tidyr::spread(Probe_Type, Value) 
  
  # Infinium II::
  sumA <- prbsII %>%
    dplyr::summarise(Total_Count=n(), 
                     PassDetpNegs_Count=count(NegsDetP<negsMin), 
                     PassDetpNegs_Percent=round(100*PassDetpNegs_Count/Total_Count,3),
                     PassDetpPoob_Count=count(PoobDetP<poobMin), 
                     PassDetpPoob_Percent=round(100*PassDetpPoob_Count/Total_Count,3) )
  
  sumB <- prbsII %>%
    dplyr::summarise_at(vars(Beta, Mval, G, R), list(min=min, 
                                                     median=median, 
                                                     max=max,
                                                     mean=mean,
                                                     sd=sd), na.rm=TRUE)
  
  sumII <- sumA %>% dplyr::full_join(sumB, by="Probe_Type") %>%
    tidyr::gather(Stat, Value, -Probe_Type) %>% 
    tidyr::unite(Probe_Type, Probe_Type, Stat, sep='_') %>% 
    tidyr::spread(Probe_Type, Value) 
  
  # Controls::
  sumA <- prbsCtl %>%
    dplyr::summarise(Total_Count=n(), 
                     PassDetpNegs_Count=count(NegsDetP<0.02), PassDetpNegs_Percent=round(100*PassDetpNegs_Count/Total_Count,3),
                     PassDetpPoob_Count=count(PoobDetP<0.2), PassDetpPoob_Percent=round(100*PassDetpPoob_Count/Total_Count,3) )
  if (verbose>=vt+2) print(sumA)
  sumB <- prbsCtl %>%
    dplyr::summarise_at(vars(Beta, Mval, G, R), list(min=min, 
                                                     median=median, 
                                                     max=max,
                                                     mean=mean,
                                                     sd=sd), na.rm=TRUE)
  if (verbose>=vt+2) print(sumB)
  
  sumCt <- sumA %>% dplyr::full_join(sumB, by="Probe_Type") %>% 
    tidyr::gather(Stat, Value, -Probe_Type) %>% 
    tidyr::unite(Probe_Type, Probe_Type, Stat, sep='_') %>% 
    tidyr::spread(Probe_Type, Value) 
  
  sumAll <- dplyr::bind_cols(sumCG, sumI, sumII, sumCt) %>%
    dplyr::select(Total_Count_CG, PassDetpNegs_Count_CG, PassDetpNegs_Percent_CG, 
                  PassDetpPoob_Count_CG, PassDetpPoob_Percent_CG, SwappedCount_CG,
                  everything())
  
  # Missing Stats::
  #  - Number of SWAPS
  if (verbose>=vt) cat("\t", glue::glue("[{funcTag}]: Finished."),"\n\n", sep='')
  
  sumAll
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                         SSET Conversion Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

ssetToSigTib = function(sset, sig, verbose=0, vt=4) {
  funcTag='ssetToTib'
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{TAB} Starting..."),"\n", sep='')
  
  dat <- NULL
  if (sig=='IG') {
    dat <- sset@IG %>% tibble::as_tibble(rownames='Probe_ID',.name_repair = "unique") %>%
      dplyr::mutate(U=as.integer(U), M=as.integer(M)) %>%
      dplyr::rename(GU=U, GM=M) %>%
      dplyr::mutate(Design_Type=sig)
  } else if (sig=='IR') {
    dat <- sset@IR %>% tibble::as_tibble(rownames='Probe_ID',.name_repair = "unique") %>%
      dplyr::mutate(U=as.integer(U), M=as.integer(M)) %>%
      dplyr::rename(RU=U, RM=M) %>%
      dplyr::mutate(Design_Type=sig)
  } else if (sig=='oobG') {
    dat <- sset@oobG %>% tibble::as_tibble(rownames='Probe_ID',.name_repair = "unique") %>%
      dplyr::mutate(U=as.integer(U), M=as.integer(M)) %>%
      dplyr::rename(GU=U, GM=M) %>%
      dplyr::mutate(Design_Type=sig)
  } else if (sig=='oobR') {
    dat <- sset@oobR %>% tibble::as_tibble(rownames='Probe_ID',.name_repair = "unique") %>%
      dplyr::mutate(U=as.integer(U), M=as.integer(M)) %>%
      dplyr::rename(RU=U, RM=M) %>%
      dplyr::mutate(Design_Type=sig)
  } else if (sig=='II') {
    dat <- sset@II %>% tibble::as_tibble(rownames='Probe_ID',.name_repair = "unique") %>%
      dplyr::mutate(U=as.integer(U), M=as.integer(M)) %>%
      dplyr::rename(RII=U, GII=M) %>%
      dplyr::mutate(Design_Type=sig)
  } else {
    stop("\n", glue::glue("[{funcTag}]: ERROR Unsupported signal request ({sig}) from sset!"),"\n\n", sep='')
  }
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{TAB} Finished."),"\n\n", sep='')
  
  dat
}

ssetToPvals = function(sset, verbose=0, vt=4) {
  funcTag <- 'ssetToPvals'
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{TAB}{TAB} Starting..."),"\n", sep='')
  
  poob <- sesame::pval(sesame::pOOBAH(sset=sset) )
  negs <- sesame::pval(sesame::detectionPnegEcdf(sset=sset) )
  detp <- cbind(poob, negs) %>%
    tibble::as_tibble(rownames='Probe_ID') %>%
    dplyr::rename(PoobDetP=poob, NegsDetP=negs)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{TAB}{TAB} Finished."),"\n\n", sep='')
  
  detp
}

bindBands = function(ssets, verbose=0, vt=4) {
  funcTag <- 'bindBands'
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{TAB}{TAB} Starting..."),"\n", sep='')
  baseI.R <- ssetToSigTib(ssets, 'oobG') %>% 
    dplyr::select(-Design_Type) %>%
    dplyr::full_join(ssetToSigTib(ssets, 'IR'), by="Probe_ID")
  
  baseI.G <- ssetToSigTib(ssets, 'IG') %>% 
    dplyr::full_join(dplyr::select(ssetToSigTib(ssets, 'oobR'), -Design_Type), by="Probe_ID")
  
  baseI <- dplyr::bind_rows(baseI.R, baseI.G) %>%
    dplyr::arrange(Probe_ID)
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{TAB}{TAB} Finished."),"\n\n", sep='')
  
  baseI
}

ssetsToTibs = function(ssets, verbose=0, vt=3) {
  funcTag <- 'ssetsToTibs'
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{TAB}{TAB} Starting..."),"\n", sep='')
  if (verbose>=vt+5) print(ssets)
  
  base <- 'base'
  swap <- 'swap'
  stopifnot(length(ssets[[base]])==1)
  stopifnot(length(ssets[[base]])==1)
  
  # Calculate Signals::
  sigsII <- ssetToSigTib(ssets[[base]], 'II') %>% dplyr::rename(G=GII, R=RII) %>%
    dplyr::arrange(Probe_ID)
  if (verbose>=vt+3) cat(glue::glue("[{funcTag}]:{TAB}{TAB} Calculated SigsII::"),"\n", sep='')
  if (verbose>=vt+3) print(sigsII)

  base_sigs <- bindBands(ssets[[base]])
  swap_sigs <- bindBands(ssets[[swap]])
  
  # Calculate Swapped Probes::
  join_sigs <- dplyr::full_join(base_sigs, swap_sigs, by="Probe_ID", suffix=c("_base", "_swap"))
  swap_cpgs <- join_sigs %>% dplyr::filter(Design_Type_base != Design_Type_swap) %>% 
    dplyr::select(Probe_ID) %>%
    dplyr::mutate(Swapped=TRUE)
  if (verbose>=vt+3) cat(glue::glue("[{funcTag}]:{TAB}{TAB} Calculated swap_cpgs::"),"\n", sep='')
  if (verbose>=vt+3) print(swap_cpgs)
  
  sigsI <- swap_cpgs %>% dplyr::right_join(swap_sigs, by="Probe_ID")
  if (verbose>=vt+3) cat(glue::glue("[{funcTag}]:{TAB}{TAB} Calculated SigsI::"),"\n", sep='')
  if (verbose>=vt+3) print(sigsI)
  
  # Calculate Swap DetP Diff::
  base_pval <- ssetToPvals(ssets[[base]])
  swap_pval <- ssetToPvals(ssets[[swap]])
  
  join_pval <- dplyr::full_join(base_pval, swap_pval, by="Probe_ID", suffix=c("_base", "_swap")) %>%
    dplyr::mutate(PoobDetP_Delta=PoobDetP_swap-PoobDetP_base,
                  NegsDetP_Delta=NegsDetP_swap-NegsDetP_base) %>%
    dplyr::select(Probe_ID, NegsDetP_swap, PoobDetP_swap, NegsDetP_Delta, PoobDetP_Delta) %>%
    dplyr::rename(NegsDetP=NegsDetP_swap, PoobDetP=PoobDetP_swap)
  
  if (verbose>=vt+2) cat(glue::glue("[{funcTag}]:{TAB}{TAB}{TAB} join_pval::"),"\n", sep='')
  if (verbose>=vt+2) print(join_pval)
  if (verbose>=vt+2) cat(glue::glue("[{funcTag}]:{TAB}{TAB}{TAB} END join_pval::"),"\n\n", sep='')
  
  # Calculate Beta Values::
  beta <- sesame::getBetas(ssets[[swap]], quality.mask=FALSE, nondetection.mask=FALSE, mask.use.tcga=FALSE, sum.TypeI=FALSE)
  mval <- sesame::BetaValueToMValue(beta)
  
  if (verbose>=vt+2) cat(glue::glue("[{funcTag}]:{TAB}{TAB} Calculated beta/mvals::"),"\n", sep='')
  # if (verbose>=vt+3) print(beta)
  # if (verbose>=vt+3) cat(glue::glue("[{funcTag}]:{TAB}{TAB} mval::"),"\n", sep='')
  # if (verbose>=vt+3) print(mval)
  # if (verbose>=vt+3) cat(glue::glue("[{funcTag}]:{TAB}{TAB} END beta/mval::"),"\n\n", sep='')
  
  call <- cbind(beta, mval) %>%
    tibble::as_tibble(rownames='Probe_ID') %>%
    dplyr::rename(Beta=beta, Mval=mval) %>%
    dplyr::full_join(join_pval, by="Probe_ID") %>%
    dplyr::mutate_if(purrr::is_double, list(round), 6) %>%
    dplyr::mutate(Probe_Type=case_when(
      stringr::str_starts(Probe_ID, 'cg')  ~ 'cg',
      stringr::str_starts(Probe_ID, 'ch')  ~ 'ch',
      stringr::str_starts(Probe_ID, 'rs')  ~ 'rs',
      stringr::str_starts(Probe_ID, 'ctl') ~ 'ctl',
      TRUE ~ 'unk'
    ))
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{TAB}{TAB} Calls::"),"\n", sep='')
  if (verbose>=vt) print(call)

  # Join Data by Design_Type::
  data <- NULL
  data$II <- call %>% dplyr::right_join(sigsII, by="Probe_ID") %>% 
    dplyr::select(Probe_ID, Probe_Type, Beta, Mval, NegsDetP, PoobDetP, NegsDetP_Delta, PoobDetP_Delta, Probe_Type, G, R)
  data$I  <- call %>% dplyr::right_join(sigsI, by="Probe_ID") %>% 
    dplyr::select(Probe_ID, Probe_Type, Design_Type, Swapped, Beta, Mval, NegsDetP, PoobDetP, NegsDetP_Delta, PoobDetP_Delta, Probe_Type, GM, GU, RM, RU)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{TAB} Finished."),"\n\n", sep='')
  
  data
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                                fetchSSet::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

fetchSSet = function(prefix, platform, manifest, method='both', desplat=NULL,
                     ssetRDS=NULL, loadRDS=TRUE, saveRDS=TRUE, overRDS=TRUE,
                     verbose=0, vt=4) {
  funcTag='fetchSSet'
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{TAB} Prefix={prefix}"),"\n", sep='')
  
  if (is.null(ssetRDS)) ssetRDS <- '/Users/bbarnes/Documents/tmp.sset.rds'
  # outDir <- dirname(ssetRDS)
  # if (!dir.exists(outDir)) dir.create(outDir, recursive=TRUE)
  # TBD: Allow Platforms that are not just EPIC::
  #  Currently we will make all internal platforms EPIC for Sesame
  #   This is to address internal Sesame Calculations
  
  loadedRDS <- FALSE
  overRDS   <- FALSE
  
  if (loadRDS && file.exists(ssetRDS)) {
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{TAB} Loading={ssetRDS}"),"\n", sep='')
    sset <- readr::read_rds(ssetRDS)
    loadedRDS <- TRUE
  } else {
    if (method=='both') {
      sset <- NULL
      
      if (verbose>=vt+5) cat("\n\nMANIFEST::\n")
      if (verbose>=vt+5) print(manifest)
      if (verbose>=vt+5) cat("MANIFEST\n\n")
      
      sset$base <- sesame::readIDATpair(prefix, platform=platform, manifest=manifest) %>%
        sesame::dyeBiasCorrTypeINorm()
      
      # sset$base <- sesame::readIDATpair(prefix, platform=platform, manifest=manifest) %>%
      #   sesame::dyeBiasCorrTypeINorm() %>%
      sset$base <- sset$base %>% sesame::noob()
      sset$base@platform <- 'EPIC'
      
      if (verbose>=vt+3) cat("\t\t", glue::glue("[{funcTag}]: BASE platform=..."),"\n", sep='')
      if (verbose>=vt+3) cat("\t\t", glue::glue("[{funcTag}]: BASE platform={sset$base@platform}"),"\n", sep='')
      if (verbose>=vt+3) print(sset$base)
      
      sset$swap <- sesame::readIDATpair(prefix, platform=platform, manifest=manifest) %>%
        sesame::dyeBiasCorrTypeINorm() %>%
        sesame::noob() %>%
        sesame::inferTypeIChannel(switch_failed=FALSE, verbose=FALSE, summary=FALSE)
      sset$swap@platform <- 'EPIC'
    } else {
      stop(glue::glue("[{funcTag}]:{TAB} Unsupported method={method}."),"\n\n", sep='')
    }
    
    if (overRDS || (saveRDS && !file.exists(ssetRDS) && !loadedRDS) ) {
      if (FALSE) {
        if (verbose>=vt) cat(glue::glue("[{funcTag}]:{TAB} Writing={ssetRDS}"),"\n", sep='')
        readr::write_rds(sset, ssetRDS, compress="gz")
      }
    }
  }
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{TAB} Loaded SSET"),"\n\n", sep='')
  
  sset
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                             idat Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

parseIdat = function(file, channel, verbose=0, vt=4) {
  funcTag <- 'parseIdat'
  stopifnot(file.exists(file))
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{TAB} Starting Parse, channel={channel}, file={file}"),"\n", sep='')
  
  dat <- NULL
  idat <- illuminaio::readIDAT(file)
  meanName  <- paste0(channel,'_Mean')
  sdName    <- paste0(channel,'_SD')
  nbeadName <- paste0(channel,'_NBeads')
  tib <- idat$Quants %>% tibble::as_tibble(rownames="Address") %>%
    dplyr::mutate(Address=as.integer(Address))
  colnames(tib) <- c('Address',meanName,sdName,nbeadName)
  if (verbose>=vt+3) print(tib)
  
  tib.nrow <- base::nrow(tib)
  dat$tib <- tib
  dat$ChipType <- idat$ChipType
  
  if (idat$ChipType=='BeadChip 8x5') dat$ChipFormat <- '8x1'
  else if (idat$ChipType=='BeadChip 12x8') dat$ChipFormat <- '12x1'
  else if (idat$ChipType=='BeadChip 24x1x4') dat$ChipFormat <- '24x1'
  else if (idat$ChipType=='BeadChip 24x1x2') dat$ChipFormat <- '24x1'
  else stop("\n",glue::glue("[{funcTag}]: ERROR: Unrecognized ChipType={idat$ChipType}"),"\n\n", sep='')
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{TAB} nrows={tib.nrow}, ChipFormat={dat$ChipFormat}, file={file}"),"\n", sep='')
  
  # TBD::Extract full code and row/column
  dat$Barcode <- idat$Barcode
  dat$Poscode <- idat$Unknowns$MostlyA
  dat$Sentrix_Name <- paste(dat$Barcode,dat$Poscode, sep='_')
  s <- stringr::str_split(gsub("^R","", dat$Poscode),'C',simplify=T,n=2)
  dat$row <- as.integer(s[1])
  dat$col <- as.integer(s[2])
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{TAB} Done Parsing row={dat$row}, col={dat$col}"),"\n", sep='')
  
  dat
}

prefixToTib = function(prefix, verbose=0, vt=4) {
  funcTag <- 'prefixToTib'
  
  dat <- NULL
  for (channel in c('Grn','Red')) {
    # Ensure Idat is Gzipped::
    channel.file <- paste0(prefix,"_",channel,".idat")
    if (file.exists(channel.file)) 
      system(paste0("gzip ",channel.file))
    channel.file <- paste0(prefix,"_",channel,".idat.gz")
    stopifnot(file.exists(channel.file))
    
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{TAB} Channel={channel}, file={channel.file}"),"\n", sep='')
    cur <- parseIdat(channel.file, channel, verbose=verbose)
    
    if (is.null(dat)) {
      dat <- cur
    } else {
      if (dat$ChipType != cur$ChipType) 
        stop("\n",glue::glue("[{funcTag}]: ERROR ChipType's don't match {dat$ChipType} != {cur$ChipType}"),"\n", sep='')
      if (dat$ChipFormat != cur$ChipFormat) 
        stop("\n",glue::glue("[{funcTag}]: ERROR ChipFormat's don't match {dat$ChipFormat} != {cur$ChipFormat}"),"\n", sep='')
      if (dat$Barcode != cur$Barcode)
        stop("\n",glue::glue("[{funcTag}]: ERROR Barcode's don't match {dat$Barcode} != {cur$Barcode}"),"\n", sep='')
      
      dat$tib <- dat$tib %>% dplyr::full_join(cur$tib, by="Address")
    }
  }
  if (verbose>=vt) cat("\t",glue::glue("[{funcTag}]: Done"),"\n\n", sep='')
  
  dat
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                          Sesame Manifest Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

getManifest = function(man, dat, verbose=0, vt=4) {
  funcTag <- 'getManifest'
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{TAB} Starting"),"\n", sep='')
  if (verbose>=vt+5) print(man)
  
  tib <- man %>% dplyr::right_join(dat, by=c("U"="Address")) %>% 
    dplyr::filter(!is.na(Probe_ID)) %>%
    dplyr::select(Probe_ID, M, U, DESIGN, COLOR_CHANNEL, col, Probe_Type, Probe_Source, Next_Base) %>%
    dplyr::arrange(Probe_ID)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{TAB} Done"),"\n\n", sep='')
  
  tib
}

addManifestGroups = function(platform='EPIC',build='hg19', manCPG=NULL, maxCnt=NULL, verbose=0, vt=3) {
  funcTag <- 'addManifestGroups'
  
  if (verbose>=vt)
    cat("\t",glue::glue("[{funcTag}]:{TAB} Starting platform={platform}, build={build}..."),"\n", sep='')
  
  stime <- system.time({
    minLen <- 100
    grpIdx <- 0
    if (is.null(manCPG))
      manCPG <- sesameManifest(platform='EPIC', build=opt$build, verbose=opt$verbosity)
    
    man.chr.tib <- manCPG %>% dplyr::arrange(seqnames, start) %>%
      dplyr::filter(!str_starts(Probe_ID, 'ctl')) %>%
      dplyr::filter(str_starts(seqnames, 'chr')) %>%
      dplyr::mutate(Local_Group_Key=0, Local_Group_Idx=0) %>% 
      dplyr::mutate(Local_Group_Key=as.integer(Local_Group_Key), 
                    Local_Group_Idx=as.integer(Local_Group_Idx)) %>%
      split(.$seqnames)
    
    if (length(man.chr.tib[['*']])>0) man.chr.tib[['*']] <- NULL
    if (verbose>=vt+2) print(names(man.chr.tib))
    
    man.ret.tib <- NULL
    # for (chr in names(man.chr.tib)) {
    man.ret.tib <- foreach (chr=names(man.chr.tib), .inorder=T, .final = function(x) setNames(x, names(man.chr.tib))) %dopar% {
      man.tib <- man.chr.tib[[chr]]
      if (!is.null(maxCnt)) man.tib <- man.chr.tib[[chr]] %>% head(n=maxCnt)
      
      nrows <- base::nrow(man.tib)
      if (verbose>=vt)
        cat("\t",glue::glue("[{funcTag}]:{TAB} chr={chr}, grpIdx={grpIdx}, nrows={nrows}"),"\n", sep='')
      
      prePos <- 0
      grpCnt <- 0
      for (idx in seq(1:nrows)) {
        curPos <- man.tib$start[idx]
        
        if (curPos-prePos < minLen) {
          grpCnt <- grpCnt + 1
        } else {
          grpCnt <- 1
          grpIdx <- grpIdx + 1
          prePos <- curPos
        }
        man.tib$Local_Group_Key[idx] <- grpIdx
        man.tib$Local_Group_Idx[idx] <- grpCnt
        
        # if (grpIdx>10) break
      }
      man.tib <- man.tib %>% dplyr::add_count(Local_Group_Key, name="Local_Group_Cnt") %>%
        dplyr::mutate(Local_Group_Key=as.integer(Local_Group_Key), 
                      Local_Group_Idx=as.integer(Local_Group_Idx))
      
      if (verbose>=vt)
        cat("\t",glue::glue("[{funcTag}]:{TAB} grpIdx={grpIdx}, prePos={prePos}"),"\n", sep='')
      # break
      man.tib
    }
  })
  man.ret.tib <- dplyr::bind_rows(man.ret.tib)
  if (verbose>=vt) {
    print(stime)
    cat("\t",glue::glue("[{funcTag}]:{TAB} Done."),"\n", sep='')
  }
  if (verbose>=vt+2) print(man.ret.tib)
  
  man.ret.tib
}

sesameManifest = function(platform='EPIC', build='hg19', verbose=0, vt=4) {
  funcTag <- 'sesameManifest'
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{TAB} Starting..."),"\n", sep='')
  name <- paste(platform,build,'manifest', sep='.')
  # tib <- suppressMessages(suppressWarnings(sesameData::sesameDataGet(name) )) %>% 
  tib <- sesameData::sesameDataGet(name) %>% 
    as.data.frame() %>% 
    rownames_to_column('Probe_ID') %>% 
    tibble::as_tibble() %>% 
    dplyr::mutate(col=NA, Masked=NA) %>%
    dplyr::mutate(
      col=case_when(
        stringr::str_starts(channel, 'G') ~ 'G',
        stringr::str_starts(channel, 'R') ~ 'R'),
      Masked=case_when(
        MASK_mapping==TRUE ~ 'Mapping',
        MASK_general==TRUE ~ 'General'
      ) ) %>%
    dplyr::mutate(M=address_B,
                  U=address_A,
                  COLOR_CHANNEL=channel,
                  Probe_Type=probeType,
                  Probe_Source='EPIC-B4',
                  Next_Base=nextBaseRef,
                  Probe_CpG_Cnt=probeCpGcnt) %>%
    dplyr::rename(DESIGN=designType) %>%
    dplyr::select(Probe_ID, M, U, DESIGN, COLOR_CHANNEL, col, Probe_Type, Probe_Source, Next_Base, 
                  seqnames, start, end, width, strand, Probe_CpG_Cnt, Masked) %>%
    dplyr::arrange(Probe_ID)
  
  # return(tib)
  
  name <- paste(platform,'address', sep='.')
  ses.all.que <- sesameData::sesameDataGet(name)
  
  # Probe_ID        M        U DESIGN COLOR_CHANNEL  col
  # ses.org.tib <- ses.all.que$ordering %>% tibble::as_tibble()
  ses.ctl.tib <- ses.all.que$controls
  ses.ctl.tib <- ses.ctl.tib %>% 
    dplyr::rename(U=Address) %>%
    dplyr::mutate(U=as.character(U)) %>%
    dplyr::mutate(U=as.integer(U),
                  M=NA,
                  M=as.integer(M),
                  Probe_ID=str_replace_all(Name, "[-\\(\\) ]+",'_'),
                  Probe_ID=str_replace_all(Probe_ID, "_+$",''),
                  Probe_ID=paste('ctl',Probe_ID,  sep='_'),
                  M=NA, 
                  COLOR_CHANNEL=as.character(Color_Channel),
                  col=NA,
                  Probe_Type=str_replace_all(Type, "\\s", '_'),
                  Probe_Source='EPIC-B4',
                  Next_Base=NA) %>% 
    dplyr::mutate(col=as.character(col),
                  Next_Base=as.character(Next_Base)) %>%
    tibble::as_tibble() %>%
    dplyr::mutate(DESIGN='II') %>%
    dplyr::select(Probe_ID, M, U, DESIGN, COLOR_CHANNEL, col, Probe_Type, Probe_Source, Next_Base)
  
  ses.all.man <- ses.ctl.tib %>% dplyr::bind_rows(tib) %>%
    dplyr::arrange(Probe_ID)
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{TAB} Done."),"\n", sep='')
  
  ses.all.man
}

# End of file
