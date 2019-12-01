
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
                                can.tib=NULL,
                                ref.tib=NULL,
                                ctls_tib=NULL,
                                vt=1, tc=0, tabsStr='') {
  funcTag <- 'singleSampleWorkflow'
  for (i in c(1:tc)) tabsStr <- paste0(tabsStr, TAB)
  
  stopifnot(length(opt$platform)>0)
  stopifnot(length(opt$verbosity)>0)
  stopifnot(!is.null(opt$platform))
  stopifnot(!is.null(opt$method))
  
  if (is.null(prefix) || length(prefix)==0) 
    stop("\n", glue::glue("[{funcTag}]: ERROR prefix is NULL"),"\n", sep='')
  
  verbose <- opt$verbosity
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting prefix={prefix}."),"\n", sep='')
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #              Extract Expected Sentrix Variables from Prefix::
  opt$sentrix_name <- BiocGenerics::basename(prefix)
  sentrix_codes <- opt$sentrix_name %>% stringr::str_split(pattern='_', simplify=TRUE, n=2) %>% 
    as.data.frame() %>% purrr::set_names('Barcode','Poscode')
  opt$sentrix_barcode <- sentrix_codes$Barcode
  opt$sentrix_poscode <- sentrix_codes$Poscode
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #               Build Current (Local) Output Directory::
  opt$curDir <- file.path(opt$outDir, opt$sentrix_barcode, opt$platform)
  if (!dir.exists(opt$curDir)) dir.create(opt$curDir, recursive=TRUE)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Current (Local) OutDir={opt$curDir}"),"\n", sep='')
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                     Explicitly Name Output Files::
  opt$outName <- paste(opt$sentrix_name,opt$platform,opt$method, sep='_')
  
  opt$idat_sig_CSV <- file.path(opt$curDir, paste0(opt$outName, '.IdatSignals.csv.gz'))
  opt$idat_sss_CSV <- file.path(opt$curDir, paste0(opt$outName, '.IdatSampleSheetShort.csv.gz'))
  
  opt$man_CSV    <- file.path(opt$curDir, paste0(opt$outName, '.SampleSesameManifest.csv.gz'))
  opt$sset_RDS   <- file.path(opt$curDir, paste0(opt$outName, '.sset.rds'))
  
  opt$prbs1_CSV  <- file.path(opt$curDir, paste0(opt$outName, '.ProbesI.tib.csv.gz'))
  opt$prbs2_CSV  <- file.path(opt$curDir, paste0(opt$outName, '.ProbesII.tib.csv.gz'))
  opt$prbs_RDS   <- file.path(opt$curDir, paste0(opt$outName, '.Probes.tib.rds'))
  
  opt$ss_CSV     <- file.path(opt$curDir, paste0(opt$outName, '.SampleSheet.csv.gz'))
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                             fetch IDAT Data::
  idat_dat <- fetchPrefix(prefix, opt$idat_sig_CSV, opt$idat_sss_CSV, 
                          loadIdat=opt$loadIDATS, saveIdat=opt$writeIDATS,
                          verbose=verbose,vt=vt+3,tc=tc+1)
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                             fetch SSET Data::
  ssets <- fetchSSet(prefix, opt=opt, man_add=man.add, idat=idat_dat, platform='custom', method=opt$method, 
                     verbose=verbose,vt=vt+4,tc=tc+1)
  if (opt$retSSET) return(ssets)
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                             fetch PRBS Data::
  prbs <- ssetsToTibs(ssets, verbose=verbose,vt=vt+4,tc=tc+1)
  if (opt$retPRBS) return(prbs)
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                             fetch CTLS Data::
  if (is.null(ctls_tib))
    ctls_tib <- man.add %>% dplyr::filter(stringr::str_starts(Probe_ID, 'ctl')) %>% dplyr::select(Probe_ID, Probe_Type)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{TAB} Extracted Control Probes"),"\n", sep='')
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                        Output Analytical Tibble::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  if (opt$writePrbRDS) {
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Writing Probe RDS={opt$prbs_RDS}"),"\n", sep='')
    readr::write_rds(prbs, opt$prbs_RDS, compress="gz")
  }
  if (opt$writePrbCSV) {
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Writing Probe(I) CSV={opt$prbs1_CSV}"),"\n", sep='')
    readr::write_csv(prbs$I,  opt$prbs1_CSV)
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Writing Probe(II) CSV={opt$prbs2_CSV}"),"\n", sep='')
    readr::write_csv(prbs$II, opt$prbs2_CSV)
  }
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                        Rejoining of Infinium Types::
  #
  #  TBD:: Optimization in Sesame or previous steps could remove this
  #
  prbs_tib <- dplyr::bind_rows(
    dplyr::select(prbs$I,  Probe_ID, Probe_Type, Beta, Mval, NegsDetP, PoobDetP, Swapped),
    dplyr::select(prbs$II, Probe_ID, Probe_Type, Beta, Mval, NegsDetP, PoobDetP) %>% dplyr::mutate(Swapped=NA)
  )

  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                         Detection P-value Summary::
  detP_sum_tib <- NULL
  detP_sum_tib <- summariseDetP(tib=prbs_tib, opt=opt, onlyCG=TRUE, verbose=verbose,vt=vt+3,tc=tc+1)

  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                      Auto Sample Detection Summary::
  auto_sum_tib <- NULL
  auto_sum_tib <- summariseAuto(tib=prbs_tib, ref=can.tib, opt=opt, onlyCG=TRUE, verbose=verbose,vt=vt+3,tc=tc+1)
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                      Bisulfite Conversion (GCT Score)::
  #
  #  TBD:: Need of optimization here!!!
  #
  GCT_list <- list()
  for (s in names(ssets)) {
    GCT_list[[s]] <- safeGCT(ssets[[s]], verbose=verbose,vt=vt+4,tc=tc+1)
  }
  GCTs_sum_tib <- dplyr::bind_rows(GCT_list) %>% set_names(paste('GCT', names(.), sep='_'))
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                      Phenotype Summary (ethnicity, sex, etc)::TBD
  phen_sum_tib <- list()
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                       Horvath Age Summary (353, 513)::TBD
  ages_sum_tib <- list()
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                   Probe Type Full Summary::TBD JUNK SEE CODE
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Building Sample Sheet/Summary Stats..."),"\n", sep='')
  prbs_sum_tib <- NULL
  prbs_sum_tib <- summarisePrbs(prbs=prbs, ctls=ctls_tib,
                                negsMin=opt$negsMinPval, poobMin=opt$poobMinPval, verbose=verbose)
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                      Combine Full Auto SampleSheet::
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Merging Summaries..."),"\n", sep='')
  sam_ss_tib <- idat_dat$sss %>%
    dplyr::bind_cols(auto_sum_tib) %>%
    dplyr::bind_cols(detP_sum_tib) %>%
    dplyr::bind_cols(GCTs_sum_tib) %>%
    dplyr::bind_cols(phen_sum_tib) %>%
    dplyr::bind_cols(ages_sum_tib) %>%
    dplyr::bind_cols(prbs_sum_tib) %>%
    dplyr::mutate_if(is.numeric, list(round), 6) %>%
    tibble::add_column(AnalysisPlatform=opt$platform) %>%
    dplyr::select(AnalysisPlatform, dplyr::everything())
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                       Write Full Auto SampleSheet::
  if (opt$writeSSheet) {
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Writing Auto Sample Sheet CSV={opt$ss_CSV}"),"\n", sep='')
    readr::write_csv(sam_ss_tib, opt$ss_CSV)
  }
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                      Return Full Data Set:: TESTING
  if (!is.null(opt$fullData) && length(opt$fullData)>0 && opt$fullData) {
    full_data <- list()
    full_data$sam_ss_tib   <- sam_ss_tib
    full_data$core_sum_tib <- idat_dat$sss
    full_data$auto_sum_tib <- auto_sum_tib
    full_data$detP_sum_tib <- detP_sum_tib
    full_data$GCTs_sum_tib <- GCTs_sum_tib
    full_data$prbs_sum_tib <- prbs_sum_tib
    full_data$core_dat_tib <- idat_dat$sig
    full_data$ssets <- ssets
    full_data$prbs  <- prbs
    full_data$ctls  <- ctls_tib
    return(full_data)
  }
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Finished."),"\n\n", sep='')

  return(sam_ss_tib)
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                    Auto Sample Detection Functions::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

firstVsRestDelta = function(mat, method, minDelta=0.2, verbose=0, vt=5) {
  funcTag <- 'firstVsRestDelta'
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{TAB}{TAB} Starting({method})..."),"\n", sep='')
  
  ncols  <- base::ncol(mat)
  nrows  <- base::nrow(mat)
  cnames <- colnames(mat)
  
  stopifnot(nrows>0)
  
  max_cnt <- 0
  max_idx <- 0
  max_per <- 0
  max_key <- 'Unknown'
  for (i in c(2:ncols)) {
    cur_cnt <- length(which(abs(rowDiffs(mat[,c(1,i)])) <= minDelta))
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{TAB}{TAB}{TAB} i={i}, matchCnt={cur_cnt}"),"\n", sep='')
    if (cur_cnt > max_cnt) {
      max_cnt <- cur_cnt
      max_idx <- i
      if (verbose>=vt) cat(glue::glue("[{funcTag}]:{TAB}{TAB}{TAB} maxCnt={max_cnt}, maxIdx={max_idx}"),"\n", sep='')
    }
  }
  if (max_idx != 0) max_key <- cnames[max_idx]
  max_per <- round(100*max_cnt/nrows, 3)
  
  keyName <- paste0('Auto_Sample_',method,'_bDB_Key')
  valName <- paste0('Auto_Sample_',method,'_bDB_Val')
  mtib <- tibble(!!keyName := max_key,
                 !!valName := max_per)
  # mtib <- tibble("Auto_Sample_DB_Key" = max_key,
  #                "Auto_Sample_DB_Val" = max_per)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{TAB}{TAB} Finished."),"\n\n", sep='')
  
  mtib
}

safeCorR2 = function(mat, verbose=0, vt=4, tc=1) {
  funcTag <- 'safeCorR2'
  tabsStr <- ''
  for (i in c(1:tc)) tabsStr <- paste0(tabsStr, TAB)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting..."),"\n", sep='')
  
  cor_val = NA
  cor_val = tryCatch({
    cor(mat, method="pearson", use="pairwise.complete.obs")
  }, warning = function(w) {
    'warning-SBM-Cor'
  }, error = function(e) {
    'error-SBM-Cor'
  }, finally = {
    'cleanup-SBM-Cor'
  })
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done."),"\n\n", sep='')  
  
  cor_val
}

predictAutoSample = function(tib, ref, minPval=0.02, minDelta=0.2, verbose=0,vt=4,tc=1,tabsStr='') {
  funcTag <- 'predictAutoSample'
  for (i in c(1:tc)) tabsStr <- paste0(tabsStr, TAB)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting..."),"\n", sep='')
  
  sbb_tib <- tib %>% dplyr::select(Probe_ID, Beta) %>%
    dplyr::left_join(ref[['BETA']]$Beta, by="Probe_ID")
  sbm_tib <- tib %>% dplyr::select(Probe_ID, Mval) %>%
    dplyr::left_join(ref[['BETA']]$Mval, by="Probe_ID")
  
  sdb_tib <- tib %>% dplyr::select(Probe_ID, Beta) %>%
    dplyr::left_join(ref[['DELTA']]$Beta, by="Probe_ID")
  sdm_tib <- tib %>% dplyr::select(Probe_ID, Mval) %>%
    dplyr::left_join(ref[['DELTA']]$Mval, by="Probe_ID")
  
  sbb_mat <- sbb_tib %>% dplyr::select(-Probe_ID) %>% as.matrix()
  sdb_mat <- sdb_tib %>% dplyr::select(-Probe_ID) %>% as.matrix()
  
  sbm_mat <- sbm_tib %>% dplyr::select(-Probe_ID) %>% as.matrix()
  sdm_mat <- sdm_tib %>% dplyr::select(-Probe_ID) %>% as.matrix()
  
  # Get Delta Beta::BETA
  dbb_tib <- firstVsRestDelta(sbb_mat, method='BETA', verbose=verbose, vt=vt)
  if (verbose>=vt) print(dbb_tib)
  
  # Get Delta Beta::DELTA
  dbd_tib <- firstVsRestDelta(sdb_mat, method='DELTA', verbose=verbose, vt=vt)
  if (verbose>=vt) print(dbd_tib)
  
  # Get R-Squared::BETA-Mval
  sbm_cor <- safeCorR2(sbm_mat, verbose=verbose, vt=vt, tc=tc)
  
  # Get R-Squared::DELTA-Mval
  sdm_cor <- safeCorR2(sdm_mat, verbose=verbose, vt=vt, tc=tc)
  
  # R-Squared with Mvales for BETA
  mbrkey <- NA
  mbrval <- NA
  if (!is.null(sbm_cor)) {
    if (verbose>=vt) print(sbm_cor[1,])
    ncols <- length(sbm_cor[1,])
    mbrval <- max(sbm_cor[1,2:ncols])
    mbridx <- which(sbm_cor[1,2:ncols]==mbrval)
    mbrkey <- names(sbm_cor[1,2:ncols])[mbridx]
    mbrval <- round(max(sbm_cor[1,2:ncols]),6)
  }
  rmb_tib <- tibble("Auto_Sample_BETA_mR2_Key" = mbrkey,
                    "Auto_Sample_BETA_mR2_Val" = mbrval)
  
  # R-Squared with Mvales for DELTA
  mdrkey <- NA
  mdrval <- NA
  if (!is.null(sdm_cor)) {
    if (verbose>=vt) print(sdm_cor[1,])
    ncols <- length(sdm_cor[1,])
    mdrval <- max(sdm_cor[1,2:ncols])
    mdridx <- which(sdm_cor[1,2:ncols]==mdrval)
    mdrkey <- names(sdm_cor[1,2:ncols])[mdridx]
    mdrval <- round(max(sdm_cor[1,2:ncols]),6)
  }
  rmd_tib <- tibble("Auto_Sample_DELTA_mR2_Key" = mdrkey,
                    "Auto_Sample_DELTA_mR2_Val" = mdrval)
  
  atib <- NULL
  atib <- dplyr::bind_cols(dbb_tib, dbd_tib,
                           rmb_tib, rmd_tib)
  
  if (verbose>=vt) print(atib)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done."),"\n\n", sep='')
  
  atib
}

summariseAuto = function(tib, ref=NULL, opt, onlyCG=TRUE,
                         verbose=0,vt=4,tc=1,tabsStr='') {
  funcTag <- 'summariseAuto'
  for (i in c(1:tc)) tabsStr <- paste0(tabsStr, TAB)
  
  auto_sum_tib <- tibble('Auto_Sample_BETA_bDB_Key'='Uknown',
                         'Auto_Sample_BETA_bDB_Val'=NA,
                         'Auto_Sample_DELTA_bDB_Key'='Uknown',
                         'Auto_Sample_DELTA_bDB_Val'=NA,
                         
                         'Auto_Sample_BETA_mR2_Key'='Unknown',
                         'Auto_Sample_BETA_mR2_Val'=NA,
                         'Auto_Sample_DELTA_mR2_Key'='Unknown',
                         'Auto_Sample_DELTA_mR2_Val'=NA
  )
  
  if (length(opt$autoDetect) && !is.null(opt$autoDetect) && opt$autoDetect && !is.null(ref)) {
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting(onlyCG={onlyCG})..."),"\n", sep='')
    
    if (onlyCG) tib <- tib %>% dplyr::filter(Probe_Type=='cg')
    tib <- tib %>% dplyr::arrange(Probe_ID)
    
    # Specail Case for NZT Pool::
    #
    if (opt$platform=='NZT') {
      tib <- tib %>% dplyr::filter(str_detect(Probe_ID, '_C_I')) %>% 
        dplyr::mutate(Probe_ID=str_remove(Probe_ID, '_[FR].*$'))
    }
    
    # 1. Detect Auto-Samples::
    auto_sum_tib <- predictAutoSample(tib=tib, ref=ref, minPval=opt$negsMinPval, minDelta=0.2, 
                                      verbose=verbose,vt=vt+1,tc=tc+1)
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Detected Auto Samples(DB)={auto_sum_tib$Auto_Sample_BETA_bDB_Key}"),"\n", sep='')
    
    # 2. Upload Reference Samples::
    
    
    # 3. Generate R2/DB Reference Comparison Stats::
    
  }
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done."),"\n\n", sep='')
  
  auto_sum_tib
}

summariseDetP = function(tib, opt, onlyCG=TRUE,
                         verbose=0,vt=4,tc=1,tabsStr='') {
  funcTag <- 'summariseDetP'
  for (i in c(1:tc)) tabsStr <- paste0(tabsStr, TAB)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting(onlyCG={onlyCG}..."),"\n", sep='')
  
  if (onlyCG) tib <- tib %>% dplyr::filter(Probe_Type=='cg')

  detP_sum_tib <- tib %>% 
    dplyr::summarise(Total_Count_CG=n(), 
                     PassDetpNegs_Count_CG=count(NegsDetP<opt$negsMinPval), 
                     PassDetpNegs_Percent_CG=round(100*PassDetpNegs_Count_CG/Total_Count_CG,3),
                     PassDetpPoob_Count_CG=count(PoobDetP<opt$poobMinPval), 
                     PassDetpPoob_Percent_CG=round(100*PassDetpPoob_Count_CG/Total_Count_CG,3),
                     SwappedCount_CG=count(!is.na(Swapped)),
                     Sample_Pass_Stringent=case_when(
                       PassDetpNegs_Percent_CG < opt$negsMinCutoff_Stringent ~ FALSE,
                       TRUE ~ TRUE
                     ),
                     Sample_Pass_Relaxed=case_when(
                       PassDetpNegs_Percent_CG < opt$negsMinCutoff_Stringent ~ FALSE,
                       TRUE ~ TRUE
                     ))
  if (verbose>=vt+5) print(detP_sum_tib)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done."),"\n\n", sep='')

  detP_sum_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                         Sample Summary Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

#
# TBD:: This functions needs to be seriously rewritten!!!
#
summarisePrbs = function(prbs, ctls, negsMin=0.02, poobMin=0.2, verbose=0, vt=5) {
  funcTag <- 'summarisePrbs'
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{TAB} Starting..."),"\n", sep='')
  
  prbsI   <- prbs$I  %>% dplyr::filter(Probe_Type!='ctl') %>% dplyr::group_by(Probe_Type)
  if (verbose>=vt+2) print(prbsI)
  
  prbsII  <- prbs$II %>% dplyr::filter(Probe_Type!='ctl') %>% dplyr::group_by(Probe_Type)
  if (verbose>=vt+2) print(prbsII)
  
  prbsCtl <- prbs$II %>% dplyr::filter(Probe_Type=='ctl') %>% dplyr::select(-Probe_Type) %>%
    dplyr::left_join(ctls, by="Probe_ID") %>%
    dplyr::select(Probe_ID, Probe_Type, everything()) %>% dplyr::group_by(Probe_Type)
  if (verbose>=vt+2) print(prbsCtl)
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                               Infinium I::
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
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                               Infinium II::
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
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                            Infinium Controls::
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
  
  sumAll <- dplyr::bind_cols(sumI, sumII, sumCt)
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{TAB} Finished."),"\n\n", sep='')
  
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

ssetsToTibs = function(ssets, verbose=0,vt=5,tc=1,tabsStr='') {
  funcTag <- 'ssetsToTibs'
  for (i in c(1:tc)) tabsStr <- paste0(tabsStr, TAB)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting..."),"\n", sep='')
  
  base <- 'base'
  swap <- 'swap'
  stopifnot(length(ssets[[base]])==1)
  stopifnot(length(ssets[[base]])==1)
  
  # Calculate Signals::
  sigsII <- ssetToSigTib(ssets[[base]], 'II') %>% dplyr::rename(G=GII, R=RII) %>%
    dplyr::arrange(Probe_ID)
  if (verbose>=vt+1) cat(glue::glue("[{funcTag}]:{tabsStr} Calculated SigsII::"),"\n", sep='')
  if (verbose>=vt+2) print(sigsII)
  
  base_sigs <- bindBands(ssets[[base]])
  swap_sigs <- bindBands(ssets[[swap]])
  
  # Calculate Swapped Probes::
  join_sigs <- dplyr::full_join(base_sigs, swap_sigs, by="Probe_ID", suffix=c("_base", "_swap"))
  swap_cpgs <- join_sigs %>% dplyr::filter(Design_Type_base != Design_Type_swap) %>% 
    dplyr::select(Probe_ID) %>%
    dplyr::mutate(Swapped=TRUE)
  if (verbose>=vt+1) cat(glue::glue("[{funcTag}]:{tabsStr} Calculated swap_cpgs::"),"\n", sep='')
  if (verbose>=vt+2) print(swap_cpgs)
  
  sigsI <- swap_cpgs %>% dplyr::right_join(swap_sigs, by="Probe_ID")
  if (verbose>=vt+1) cat(glue::glue("[{funcTag}]:{tabsStr} Calculated SigsI::"),"\n", sep='')
  if (verbose>=vt+2) print(sigsI)
  
  # Calculate Swap DetP Diff::
  base_pval <- ssetToPvals(ssets[[base]])
  swap_pval <- ssetToPvals(ssets[[swap]])
  
  join_pval <- dplyr::full_join(base_pval, swap_pval, by="Probe_ID", suffix=c("_base", "_swap")) %>%
    dplyr::mutate(PoobDetP_Delta=PoobDetP_swap-PoobDetP_base,
                  NegsDetP_Delta=NegsDetP_swap-NegsDetP_base) %>%
    dplyr::select(Probe_ID, NegsDetP_swap, PoobDetP_swap, NegsDetP_Delta, PoobDetP_Delta) %>%
    dplyr::rename(NegsDetP=NegsDetP_swap, PoobDetP=PoobDetP_swap)
  
  if (verbose>=vt+3) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB}{TAB} join_pval::"),"\n", sep='')
  if (verbose>=vt+3) print(join_pval)
  if (verbose>=vt+3) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB}{TAB} END join_pval::"),"\n\n", sep='')
  
  # Calculate Beta Values::
  beta <- sesame::getBetas(ssets[[swap]], quality.mask=FALSE, nondetection.mask=FALSE, mask.use.tcga=FALSE, sum.TypeI=FALSE)
  mval <- sesame::BetaValueToMValue(beta)
  
  if (verbose>=vt+1) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Calculated beta/mvals::"),"\n", sep='')
  
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
  if (verbose>=vt+1) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Calls::"),"\n", sep='')
  if (verbose>=vt+2) print(call)
  
  # Join Data by Design_Type::
  data <- NULL
  data$II <- call %>% dplyr::right_join(sigsII, by="Probe_ID") %>% 
    dplyr::select(Probe_ID, Probe_Type, Beta, Mval, NegsDetP, PoobDetP, NegsDetP_Delta, PoobDetP_Delta, Probe_Type, G, R)
  data$I  <- call %>% dplyr::right_join(sigsI, by="Probe_ID") %>% 
    dplyr::select(Probe_ID, Probe_Type, Design_Type, Swapped, Beta, Mval, NegsDetP, PoobDetP, NegsDetP_Delta, PoobDetP_Delta, Probe_Type, GM, GU, RM, RU)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Finished."),"\n\n", sep='')
  
  data
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                                fetchSSet::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

fetchSSet = function(prefix, opt, man_add, idat=NULL, 
                     method=NULL, platform=NULL, setPlatform='EPIC',
                     verbose=0,vt=4,tc=1,tabsStr='') {
  funcTag='fetchSSet'
  for (i in c(1:tc)) tabsStr <- paste0(tabsStr, TAB)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting prefix={prefix}"),"\n", sep='')
  
  # TBD::Forcing EPIC for Now...
  stopifnot(setPlatform=='EPIC')
  
  if (is.null(platform)) platform <- opt$platform
  if (is.null(method)) method <- opt$method
  
  # First Validate that idat_sss_CSV exist others loadSSET <- FALSE
  if (is.null(idat)) {
    opt$loadSSET <- FALSE
    idat <- fetchIdats(prefix=prefix, sig_csv=opt$idat_sig_CSV, sss_csv=opt$idat_sss_CSV, 
                       loadIdat=opt$loadIDATS, saveIdat=opt$writeIDATS,
                       verbose=verbose,vt=vt+1,tc=tc+1)
  }
  
  if (opt$loadSSET && !is.null(opt$sset_RDS) && file.exists(opt$sset_RDS) && 
      file.mtime(opt$idat_sss_CSV) < file.mtime(opt$sset_RDS)) {
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Loading={opt$sset_RDS}"),"\n", sep='')
    sset <- readr::read_rds(opt$sset_RDS)
  } else {
    # TBD:: Currently Always Calculate Manifest
    manifest <- getManifest(man=man_add, dat=idat$sig, verbose=verbose,vt=vt+1,tc=tc+1)
    stopifnot(!is.null(manifest), length(manifest)>0)
    # if (opt$writeMIDMAN) {
    #   cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Writing Sample-Mid-Manifest={opt$man_CSV}"),"\n", sep='')
    #   readr::write_csv(manifest, opt$man_CSV)
    # }
    
    if (method=='both') {
      sset <- NULL
      
      sset$base <- sesame::readIDATpair(prefix, platform=platform, manifest=manifest) %>%
        sesame::dyeBiasCorrTypeINorm() %>%
        sesame::noob()
      sset$base@platform <- setPlatform
      
      sset$swap <- sesame::readIDATpair(prefix, platform=platform, manifest=manifest) %>%
        sesame::dyeBiasCorrTypeINorm() %>%
        sesame::noob() %>%
        sesame::inferTypeIChannel(switch_failed=FALSE, verbose=FALSE, summary=FALSE)
      sset$swap@platform <- setPlatform
    } else {
      stop(glue::glue("[{funcTag}]: Fatal ERROR Unsupported method={method}."),"\n\n", sep='')
      q()
    }
    
    if (opt$saveRDS) {
      if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Writing={opt$sset_RDS}"),"\n", sep='')
      readr::write_rds(sset, opt$sset_RDS, compress="gz")
    }
  }
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Loaded SSET"),"\n\n", sep='')
  
  sset
}

# End of file
