
suppressWarnings(suppressPackageStartupMessages(require("tidyverse")) )
suppressWarnings(suppressPackageStartupMessages(require("stringr")) )
suppressWarnings(suppressPackageStartupMessages(require("glue")) )

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
#                       Single Sample Workflow::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
sesamizeSingleSample = function(prefix, man, add, autoRef, opt, retData=FALSE, del='_', vt=2,tc=0) {
  funcTag <- 'sesamizeSingleSample'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (opt$verbosity>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} prefix={prefix}.{RET}"))
  
  # TBD:: Validate all options are present in opt!!!
  #
  
  # tTracker <- timeTracker$new(verbose=opt$verbosity)
  tTracker <- timeTracker$new()
  stime <- system.time({
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                              Define Outputs::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    basecode <- basename(prefix)
    out_name <- paste(basecode, opt$platform, opt$manifest, sep=del)
    
    call_csv <- file.path(opt$outDir, paste(out_name, 'call.csv.gz', sep=del) )
    sigs_csv <- file.path(opt$outDir, paste(out_name, 'sigs.csv.gz', sep=del) )
    samp_csv <- file.path(opt$outDir, paste(out_name, 'AutoSampleSheet.csv.gz', sep=del) )
    ssum_csv <- file.path(opt$outDir, paste(out_name, 'SignalSummary.csv.gz', sep=del) )
    time_csv <- file.path(opt$outDir, paste(out_name, 'runTime.csv.gz', sep=del) )
    
    stampBeg_txt <- file.path(opt$outDir, paste(out_name, 'timestamp.beg.txt', sep=del) )
    stampEnd_txt <- file.path(opt$outDir, paste(out_name, 'timestamp.end.txt', sep=del) )
    system(glue::glue('touch {stampBeg_txt}') )
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                           Extract Raw idat::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #
    # TBD:: Validate Total-Beads, Average-Bead-Reps
    idat <- prefixToIdat(prefix, verbose=opt$verbosity,vt=vt+1,tc=tc+1,tt=tTracker)
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                     Get Sesame Manifest/Address Tables::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    ses_man_tib <- getSesameManifest(man, idat$sig, verbose=opt$verbosity,vt=1,tc=tc+1,tt=tTracker)
    ses_add_tib <- ses_man_tib %>% dplyr::select(Probe_ID) %>% dplyr::inner_join(add, by='Probe_ID')

    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                     Initialize Data (Structures) Tables::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    # Raw Data Tables::
    #  - call_tib: table by [ Probe_ID,Design_Type,Probe_Type ] -> [ detp-vals, beta-values] i.e. pval/beta
    #  - sigs_tib: table by [ Probe_ID,Address,Man_Col,Design_Type,Probe_Type ] -> [ Swap, G, R ] i.e. intensities
    call_tib  <- ses_man_tib %>% dplyr::select(Probe_ID, DESIGN, Probe_Type) %>% dplyr::rename(Design_Type=DESIGN)
    sigs_tib  <- ses_add_tib
    
    # Summary and Intermediate Summary Tables::
    #  - samp_sum_tib: complete summary table. i.e. Automatic Sample Sheet
    #  - sigs_sum_tib: complete intensity summary table.
    #  - phen_sum_tib: intermediate summary table [ GCT, HorvathAge, Ethnicity, etc. ] i.e. predictive measuements
    #  - bead_sum_tib: intermeidate summary table [ CG_Bead_Count,CG_Bead_Total,CG_Bead_AvgRep ] i.e. idat summary
    #  - pool_sum_tib: intermeidate summary table [ Bead_Pool,cg_Manifest_Count,ch_Manifest_Count,rs_Manifest_Count ] i.e. manifest summary
    samp_sum_tib <- NULL
    sigs_sum_tib <- NULL
    phen_sum_tib <- NULL
    bead_sum_tib <- ses_add_tib %>% dplyr::filter(Probe_Type=='cg') %>% 
      dplyr::select(Address) %>% dplyr::left_join(idat$sig, by="Address") %>% 
      summarise(CG_Bead_Count=n(), CG_Bead_Total=sum(Bead_Grn,Bead_Red), CG_Bead_AvgRep=CG_Bead_Total/CG_Bead_Count/2)
    pool_sum_tib <- ses_man_tib %>% dplyr::filter(Probe_Type=='cg' | Probe_Type=='ch' | Probe_Type=='rs') %>% 
      dplyr::mutate(Probe_Type=stringr::str_to_upper(Probe_Type)) %>%
      dplyr::group_by(Probe_Type) %>%
      dplyr::summarise(Count=n()) %>% tidyr::spread(Probe_Type,Count) %>% 
      purrr::set_names(paste(names(.),'Manifest_Count',sep='_') ) %>% 
      addBeadPoolToSampleSheet(field='CG_Manifest_Count') %>% dplyr::ungroup()
    
    beadPool <- NULL
    beadPool <- pool_sum_tib %>% select(Bead_Pool) %>% dplyr::pull()
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                             Build Raw SSET::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    raw_sset <- initSesameRaw(prefix=prefix, platform=opt$platform, manifest=ses_man_tib,
                              verbose=opt$verbosity,vt=vt+1,tc=tc+1,tt=tTracker)
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                     'Dye-Swap-Noob' Order of Operations::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    if (FALSE) {
      sCalls <- c('dyeBiasCorrTypeINorm', 'inferTypeIChannel', 'noob')
      nCalls <- c(FALSE, TRUE,  TRUE)
      pCalls <- c(FALSE, TRUE,  TRUE)
      bCalls <- c(FALSE, FALSE, TRUE)
      iCalls <- c(FALSE, FALSE, TRUE)
      fCalls <- c(FALSE, FALSE, FALSE)
      
      ses_data <- sesameWorkflow(sset=raw_sset, add=add, call=call_tib, sigs=sigs_tib, pheno=phen_sum_tib,
                                 stepCalls=sCalls, negsCalls=nCalls,poobCalls=pCalls, betaCalls=bCalls, intsCalls=iCalls, fenoCalls=fCalls,
                                 verbose=opt$verbosity,vt=vt+1,tc=tc+1,tt=tTracker)
      
      call_tib <- ses_data[[1]]
      sigs_tib <- ses_data[[2]]
      phen_sum_tib <- ses_data[[3]]
    }
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                  SwapOpen:: 'Swap-Noob-Dye' Order of Operations::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    if (FALSE) {
      sCalls <- c('inferTypeIChannel','noob','dyeBiasCorrTypeINorm')
      
      nCalls <- c(TRUE,  FALSE, TRUE)
      pCalls <- c(TRUE,  FALSE, TRUE)
      bCalls <- c(FALSE, FALSE, TRUE)
      iCalls <- c(FALSE, FALSE, TRUE)
      fCalls <- c(FALSE, FALSE, FALSE)
      
      ses_data <- sesameWorkflow(sset=raw_sset, add=add, call=call_tib, sigs=sigs_tib, pheno=phen_sum_tib,
                                 stepCalls=sCalls, negsCalls=nCalls,poobCalls=pCalls, betaCalls=bCalls, intsCalls=iCalls, fenoCalls=fCalls,
                                 verbose=opt$verbosity,vt=vt+1,tc=tc+1,tt=tTracker)
      
      call_tib <- ses_data[[1]]
      sigs_tib <- ses_data[[2]]
      phen_sum_tib <- ses_data[[3]]
    }
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #               OpenSesame:: 'Noob-Dye-Swap' Order of Operations::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    if (opt$verbosity>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Building 'OpenSesame'.{RET}"))
    
    if (FALSE) {
      # NOTES: This is what OpenSesame does. The main difference is extracting detection p-values first. 
      #
      core_oo_str  <- 'RNDI'
      
      sCalls <- c('Raw' ,'noob', 'dyeBiasCorrTypeINorm', 'inferTypeIChannel')
      nCalls <- c(TRUE,  FALSE, TRUE,  TRUE)
      pCalls <- c(TRUE,  FALSE, TRUE,  TRUE)
      bCalls <- c(FALSE, FALSE, TRUE,  TRUE)
      iCalls <- c(FALSE, FALSE, FALSE, TRUE)
      fCalls <- c(FALSE, FALSE, TRUE,  FALSE)
    } else {
      # NOTES: This is the same as above, but not calculating detection p-values until after backgorund(noob) and
      #  dye-bias-correction. Not sure there is a big difference. For simplicity we'll stick with this one. 
      #
      core_oo_str  <- 'NDI'
      
      sCalls <- c('noob', 'dyeBiasCorrTypeINorm', 'inferTypeIChannel')
      nCalls <- c(FALSE, TRUE,  TRUE)
      pCalls <- c(FALSE, TRUE,  TRUE)
      bCalls <- c(FALSE, TRUE,  TRUE)
      iCalls <- c(FALSE, FALSE, TRUE)
      fCalls <- c(FALSE, TRUE,  FALSE)
    }
    requeue_key  <- paste('CG',core_oo_str,'negs_pval_PassPerc', sep=del)
    core_oo_beta <- paste(core_oo_str,'beta',sep=del)
    core_oo_negs <- paste(core_oo_str,'negs','pval',sep=del)
    core_oo_poob <- paste(core_oo_str,'poob','pval',sep=del)
    
    ses_data <- sesameWorkflow(sset=raw_sset, add=add, call=call_tib, sigs=sigs_tib, pheno=phen_sum_tib, beadPool=beadPool,
                               stepCalls=sCalls, negsCalls=nCalls,poobCalls=pCalls, betaCalls=bCalls, intsCalls=iCalls, fenoCalls=fCalls,
                               verbose=opt$verbosity,vt=vt+1,tc=tc+1,tt=tTracker)

    call_tib <- ses_data[[1]]
    sigs_tib <- ses_data[[2]]
    phen_sum_tib <- ses_data[[3]]
    # sset_ret <- ses_data[[4]] # For Development Mode...
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                            Build Sample Sheet::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    if (opt$verbosity>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Building Auto-Sample-Sheet: [idat/bead/pheno].{RET}"))
    
    samp_sum_tib <- idat$ann
    samp_sum_tib <- samp_sum_tib %>% dplyr::bind_cols(bead_sum_tib)
    samp_sum_tib <- samp_sum_tib %>% dplyr::bind_cols(pool_sum_tib)
    samp_sum_tib <- samp_sum_tib %>% dplyr::bind_cols(phen_sum_tib)
    samp_sum_tib <- samp_sum_tib %>% add_column(minPvalUsed = opt$minPval)
    samp_sum_tib <- samp_sum_tib %>% add_column(minDeltaUsed = opt$minDelta)
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                           Auto-Detect Sample::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    if (opt$autoDetect) {
      if (opt$verbosity>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Detecting Auto-Sample-Sheet: [SampleName, rsquared,deltaBeta].{RET}"))
      
      auto_data <- sampleDetect(can=call_tib, ref=autoRef, minPval=opt$minPval, minDelta=opt$minDelta,
                                dname='Design_Type', pname='Probe_Type', ptype='cg',
                                join='Probe_ID', field=core_oo_beta, pval=core_oo_negs, suffix='beta', del=del,
                                outDir=opt$outDir, sname=out_name, plotMatrix=opt$plotAuto, writeMatrix=opt$writeAuto,
                                dpi=opt$dpi, format=opt$plotFormat, datIdx=4,
                                verbose=opt$verbosity,vt=vt+1,tc=tc+1,tt=tTracker)
      
      samp_sum_tib <- samp_sum_tib %>% dplyr::bind_cols(auto_data)
    }
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                        Summarize Sample Sheet Stats::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    if (opt$verbosity>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Summarizing Auto-Sample-Sheet: [pval/beta].{RET}"))
    
    samp_sum_tib <- samp_sum_tib %>%
      dplyr::bind_cols(
        pval_sums <- call_tib %>% dplyr::filter(Probe_Type=='cg') %>%
          dplyr::select(ends_with('_pval')) %>%
          dplyr::summarise_all(list(cntPer_lte), min=opt$minPval) %>%
          purrr::set_names(paste0('CG_',names(.),'_PassPerc') ),
        
        beta_sums <- call_tib %>% dplyr::filter(Probe_Type=='cg') %>%
          dplyr::select(ends_with('_beta')) %>%
          # dplyr::summarise_all(list(min=min, avg=mean, med=median, max=max, SD=sd), na.rm=TRUE ) %>%
          dplyr::summarise_all(list(avg=mean, med=median, SD=sd), na.rm=TRUE ) %>%
          purrr::set_names(paste0('CG_',names(.),'_Beta') )
      )
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                        Summarize Sample Intensities::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    if (opt$verbosity>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Summarizing Intensities: [min/avg/med/max/sd].{RET}"))
    
    sigs_sum_tib <- sigsToSummary(sigs_tib, name=core_oo_str, verbose=opt$verbosity,vt=vt+1,tc=tc+1,tt=tTracker)
    
    # NOTES: This is a bit of hacking to get ONLY the fields (defualt=avg) for each probe/design type into the 
    #  Auto-SampleSheet. It puts everything in a human readable order with only the most relevant field (default=avg).
    #  However, all this stuff is fully summarized in the intentisity summary table above. 
    #
    sigField <- paste0(del,opt$sigs_sum_field)
    samp_sum_tib <- samp_sum_tib %>% dplyr::bind_cols(
      dplyr::bind_cols(
        sigs_sum_tib %>% dplyr::filter(Design_Type!='II') %>%
          dplyr::select(Design_Type,Probe_Type,Man_Col, ends_with('_Count'), ends_with(sigField)) %>%
          tidyr::unite(Type, Design_Type,Probe_Type,Man_Col, sep=del) %>% 
          tidyr::gather(Metric, Value, -Type) %>% 
          tidyr::unite(Key, Type,Metric, sep=del) %>% 
          spread(Key, Value) %>%
          purrr::set_names( stringr::str_replace(names(.), '^([MU]I)_([chgrs][chgrs])_', '\\$2_\\$1_') %>% stringr::str_remove_all('\\\\') ),
        
        sigs_sum_tib %>% dplyr::filter(Design_Type=='II') %>% dplyr::select(-contains('_Swap_')) %>%
          dplyr::filter(Probe_Type=='cg' | Probe_Type=='ch' | Probe_Type=='rs') %>%
          dplyr::select(Design_Type,Probe_Type, ends_with('_Count'), ends_with(sigField)) %>%
          tidyr::unite(Type, Design_Type,Probe_Type, sep=del) %>% 
          tidyr::gather(Metric, Value, -Type) %>% 
          tidyr::unite(Key, Type,Metric, sep=del) %>% 
          spread(Key, Value) %>%
          purrr::set_names( stringr::str_replace(names(.), '^(II)_([chgrs][chgrs])_', '\\$2_\\$1_') %>% stringr::str_remove_all('\\\\') ),
        
        sigs_sum_tib %>% dplyr::filter(Design_Type=='II') %>% dplyr::select(-contains('_Swap_')) %>%
          dplyr::filter(Probe_Type!='cg', Probe_Type!='ch', Probe_Type!='rs') %>%
          dplyr::select(Design_Type,Probe_Type, ends_with('_Count'), ends_with(sigField)) %>%
          tidyr::unite(Type, Design_Type,Probe_Type, sep=del) %>% 
          tidyr::gather(Metric, Value, -Type) %>% 
          tidyr::unite(Key, Type,Metric, sep=del) %>% 
          spread(Key, Value) %>%
          purrr::set_names( stringr::str_remove(names(.), '^(II)_') )
      ) %>% dplyr::select(starts_with('cg_'), starts_with('ch_'), starts_with('rs_'), everything())
    )

    pass_str <- 'PASSED'
    fail_str <- paste0('FAILED_',requeue_key,'<',opt$negsMinCutoff)
    requeue_key <- requeue_key %>% rlang::sym()
    samp_sum_tib <- samp_sum_tib %>% dplyr::mutate(Requeue=case_when(
      !!requeue_key < opt$negsMinCutoff ~ fail_str,
      TRUE ~ pass_str
    )) %>% dplyr::select(Requeue, everything())
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                               Write Outputs::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    if (opt$verbosity>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Writing Outputs: [call,sigs,signalSummary,sampleSheet,time].{RET}"))
    
    time_data <- tTracker$time
    
    roundData <- TRUE
    if (roundData) {
      call_tib     <- call_tib %>% dplyr::ungroup() %>% dplyr::mutate_if(is.numeric, round, 6)
      sigs_tib     <- sigs_tib %>% dplyr::ungroup() %>% dplyr::mutate_if(is.numeric, round, 6)
      samp_sum_tib <- samp_sum_tib %>% dplyr::ungroup() %>% dplyr::mutate_if(is.numeric, round, 6)
      sigs_sum_tib <- sigs_sum_tib %>% dplyr::ungroup() %>% dplyr::mutate_if(is.numeric, round, 6)
      time_data    <- time_data %>% dplyr::ungroup() %>% dplyr::mutate_if(is.numeric, round, 6)
    }
    
    if (opt$writeCall) readr::write_csv(call_tib, call_csv)
    if (opt$writeSigs) readr::write_csv(sigs_tib, sigs_csv)
    
    readr::write_csv(samp_sum_tib, samp_csv)
    readr::write_csv(sigs_sum_tib, ssum_csv)
    readr::write_csv(time_data, time_csv)
    
    system(glue::glue('touch {stampEnd_txt}') )
  })
  if (opt$verbosity>=vt) tTracker %>% print()
  if (opt$verbosity>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}{RET}"))
  
  if (retData) {
    ret <- NULL
    ret$idat   <- idat$sig
    # ret$sset   <- sset_ret
    
    ret$call   <- call_tib
    ret$sigs   <- sigs_tib
    
    ret$samSum  <- samp_sum_tib
    ret$sigSum  <- sigs_sum_tib
    ret$fenSum  <- phen_sum_tib
    
    ret$time   <- time_data
    
    return(ret)
  }
  
  # tTracker
  samp_sum_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                     Summarizing Counting Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
cntPer_lte = function(x, min) {
  round(100*length(x[which(x<=min)])/length(x),3)
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                              Pval Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

maskTibs = function(tib, field, pval, minPval, del='_',
                    verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'maskTibs'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  pval_del <- paste0(del,pval)
  snames <- tib %>% dplyr::select(ends_with(pval_del)) %>% names() %>% stringr::str_remove(pval_del)
  
  # mask_tib <- tib
  for (sname in snames) {
    fstr <- paste(sname,field, sep=del)
    pstr <- paste(sname,pval, sep=del)
    tib <- maskTib(tib, field=fstr, pval=pstr, minPval=minPval,
                   verbose=opt$verbosity,vt=vt+1,tc=tc+1,tt=tTracker)
  }
  tib
}

maskTib = function(tib, field, pval, minPval,
                   verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'maskTib'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  field <- field %>% rlang::sym()  
  pval  <- pval %>% rlang::sym()
  # minPval <- minPval %>% rlang::sym()
  
  stime <- system.time({
    tib <- tib %>% dplyr::mutate(!!field := case_when(!!pval >= !!minPval ~ NA_real_, TRUE ~ !!field))
  })
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}{RET}"))
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  
  tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                             Manifest I/O::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
loadManifestSource = function(file, verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'loadManifestSource'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Loading file={file}.{RET}"))
  
  stime <- system.time({
    tib <- NULL
    tib <- suppressMessages(suppressWarnings(readr::read_csv(file)))
  })
  if (verbose>vt+2) print(tib)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}"))
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  
  tib
}

loadAddressSource = function(file, man, fresh=FALSE, save=TRUE, split=FALSE,
                             verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'loadAddressSource'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  stime <- system.time({
    tibs <- NULL
    if (!fresh && file.exists(file)) {
      if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Loading file={file}.{RET}"))
      ## TBD:: CHeck Type!!!
      if (stringr::str_ends(file,'.rds')) {
        tibs <- suppressMessages(suppressWarnings(readr::read_rds(file) ))
      } else {
        tibs <- suppressMessages(suppressWarnings(readr::read_csv(file) ))
      }
    } else {
      if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Building file={file}.{RET}"))
      
      if (split) {
        tibs[['I1']] <- man %>% dplyr::filter(DESIGN=='I') %>% 
          dplyr::select(Probe_ID,Probe_Type,U,M) %>% 
          tidyr::gather(Design_Type,Address, -Probe_Type, -Probe_ID) %>% 
          dplyr::arrange(Address)
        
        tibs[['I2']] <- man %>% dplyr::filter(DESIGN=='II') %>% 
          dplyr::rename(Address=U) %>% dplyr::select(Probe_ID,Probe_Type,Address) %>%
          dplyr::arrange(Address)
      } else {
        tibs <- dplyr::bind_rows(
          man %>% dplyr::filter(DESIGN=='I') %>% 
            dplyr::select(Probe_ID,Probe_Type,U,M, col) %>% 
            tidyr::gather(Design_Type,Address, -Probe_Type, -Probe_ID, -col) %>%
            dplyr::rename(Man_Col=col) %>% 
            dplyr::mutate(Design_Type=paste0(Design_Type,'I')) %>%
            dplyr::select(Probe_ID,Address,Man_Col,Design_Type,Probe_Type),
          
          man %>% dplyr::filter(DESIGN=='II') %>% 
            dplyr::rename(Address=U, Man_Col=col) %>%
            dplyr::mutate(Design_Type='II') %>%
            dplyr::select(Probe_ID,Address,Man_Col,Design_Type,Probe_Type)
        ) %>% dplyr::arrange(Probe_ID)
      }
      
      if (save) {
        if (stringr::str_ends(file,'.rds')) {
          if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Saving file(RDS)={file}.{RET}"))
          readr::write_rds(tibs, file, compress='gz')
        } else {
          if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Saving file(CSV)={file}.{RET}"))
          readr::write_csv(tibs, file)
        }
      }
    }
  })
  if (verbose>=vt+2) print(tibs)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}"))
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  
  tibs
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                             Idat File I/O::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

prefixesToChipTib = function(prefixes) {
  basecodes <- names(prefixes)
  
  tib <- NULL
  for (basecode in basecodes) {
    sentrix_codes <- basecode %>% stringr::str_split(pattern='_', n=2, simplify=TRUE)
    tib <- tib %>% dplyr::bind_rows(
      tibble::tibble(barcode=sentrix_codes[1],
                     poscode=sentrix_codes[2],
                     basecode=basecode,
                     path=prefixes[[basecode]])
    )
  }
  
  tib
}

prefixToIdat = function(prefix, gzip=TRUE, validate=TRUE, verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'prefixToIdat'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Loading prefix={prefix}.{RET}"))
  
  stime <- system.time({
    dat <- NULL
    
    grn_idat <- loadIdat(prefix, 'Grn', gzip=gzip, verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    red_idat <- loadIdat(prefix, 'Red', gzip=gzip, verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                          Extract Signal Data::
    grn_sig <- getIdatSignalTib(grn_idat, channel='Grn', verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    red_sig <- getIdatSignalTib(grn_idat, channel='Red', verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    
    sigs <- dplyr::full_join(grn_sig, red_sig, by="Address")
    
    grn_ann <- dplyr::bind_cols(
      getIdatBarcodeTib(grn_idat, verbose=verbose,vt=vt+1,tc=tc+1,tt=tt),
      getIdatFormatTib(grn_idat, verbose=verbose,vt=vt+1,tc=tc+1,tt=tt),
      getIdatTimeStampTib(grn_idat, method='Decoding', verbose=verbose,vt=vt+1,tc=tc+1,tt=tt),
      getIdatTimeStampTib(grn_idat, method='Extract', verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    )
    
    if (validate) {
      red_ann <- dplyr::bind_cols(
        getIdatBarcodeTib(red_idat, verbose=verbose,vt=vt+1,tc=tc+1,tt=tt),
        getIdatFormatTib(red_idat, verbose=verbose,vt=vt+1,tc=tc+1,tt=tt),
        getIdatTimeStampTib(red_idat, method='Decoding', verbose=verbose,vt=vt+1,tc=tc+1,tt=tt),
        getIdatTimeStampTib(red_idat, method='Extract', verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
      )
      
      unique_cnt <- dplyr::bind_rows(grn_ann, red_ann) %>% dplyr::distinct() %>% base::nrow()
      stopifnot(unique_cnt==1)
      if (verbose>=vt+1) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Passed idat pair validation!{RET}"))
    }
    
    dat$sig <- sigs
    dat$ann <- grn_ann
  })
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}{RET}"))
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  
  dat
}

loadIdat = function(prefix, col, gzip=TRUE, verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'loadIdat'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Loading prefix={prefix}.{RET}"))
  
  stime <- system.time({
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                        Ensure Idats are Gzipped::
    idat_file <- paste0(prefix,"_",col,".idat")
    if (file.exists(idat_file)) {
      if (gzip) {
        if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Compressing file={idat_file}.{RET}"))
        system(paste0("gzip ",idat_file))
        idat_file <- paste0(prefix,"_",col,".idat.gz")
      }
    } else {
      idat_file <- paste0(prefix,"_",col,".idat.gz")
    }
    if (!file.exists(idat_file)) stop(glue::glue("{RET}[{funcTag}]: ERROR: idat_file={idat_file} does NOT exist!!!{RET}{RET}"))
    stopifnot(file.exists(idat_file))
    
    idat <- illuminaio::readIDAT(idat_file)
  })
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}{RET}"))
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  
  idat
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                             Idat Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

getIdatSignalTib = function(idat, channel, del='_', verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'getIdatSignalTib'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  stime <- system.time({
    datTag   <- 'sig'
    meanName <- paste('Raw', channel,datTag, sep=del)
    sdName   <- paste('SD',  channel,datTag, sep=del)
    beadName <- paste('Bead',channel, sep=del)
    
    tib <- idat$Quants %>% tibble::as_tibble(rownames="Address") %>%
      dplyr::mutate(Address=as.integer(Address)) %>%
      purrr::set_names('Address',meanName,sdName,beadName)
    tib_nrow <- base::nrow(tib)
  })
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done. Idat nrows={tib_nrow}.{RET}{RET}"))
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  
  tib
}

getIdatBarcodeTib = function(idat, verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'getIdatBarcodeTib'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  stime <- system.time({
    sentrixName <- paste(idat$Barcode,idat$Unknowns$MostlyA, sep='_')
    rowcol_df <- idat$Unknowns$MostlyA %>% stringr::str_remove("^R") %>% stringr::str_split('C', simplify=TRUE, n=2) %>% as.numeric()
    tib <- tibble::tibble(Sentrix_Name=sentrixName,
                          Sentrix_Barcode=idat$Barcode,
                          Sentrix_Poscode=idat$Unknowns$MostlyA,
                          Sentrix_Row=as.integer(rowcol_df[1]),
                          Sentrix_Col=as.integer(rowcol_df[2]) )
  })
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done. Sentrix_Name(R{tib$Sentrix_Row}:C{tib$Sentrix_Col})={tib$Sentrix_Name}{RET}{RET}"))
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  
  tib # %>% gather()
}

getIdatFormatTib = function(idat, verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'getIdatFormatTib'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  stime <- system.time({
    chipType <- idat$ChipType
    if (chipType=='BeadChip 8x5') chipFormat         <- '8x1'
    else if (chipType=='BeadChip 12x8') chipFormat   <- '12x1'
    else if (chipType=='BeadChip 24x1x4') chipFormat <- '24x1'
    else if (chipType=='BeadChip 24x1x2') chipFormat <- '24x1'
    else if (chipType=='Beadchip 24x1x2') chipFormat <- '24x1'
    else stop(glue::glue("{RET}[{funcTag}]: ERROR: Unrecognized ChipType={chipType}!{RET}{RET}"))
    tib <- tibble::tibble('ChipType'=chipType,'ChipFormat'=chipFormat)
  })
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} ChipType={chipType}, ChipFormat={chipFormat}.{RET}{RET}"))
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  
  tib # %>% gather()
}

getIdatTimeStampTib = function(idat, method='Extract', sherlockID='sherlockID', order='latest', 
                               verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'getIdatTimeStampTib'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  stime <- system.time({
    time_tib  <- idat$RunInfo %>% tibble::as_tibble()
    method_idxs <- grep(method, time_tib$BlockType)
    machine_idxs <- grep(sherlockID, time_tib$BlockPars)
    stopifnot(length(method_idxs)>0)
    if (verbose>=vt+5) print(method_idxs)
    
    stopifnot(order=='latest')
    if (order=='latest')  {
      metod_idx <- method_idxs %>% tail(n=1)
      machine_idx <- machine_idxs %>% tail(n=1)
    } else {
      metod_idx <- method_idxs %>% head(n=1)
      machine_idx <- machine_idxs %>% head(n=1)
    }
    mach_vec <- time_tib$BlockPars[machine_idx] %>% str_split('\\|', simplify=TRUE)
    mach_var <- mach_vec[2] %>% stringr::str_split('=', simplify=TRUE)
    date_str <- time_tib$RunTime[metod_idx] %>% as.POSIXct(format="%m/%d/%Y %H:%M:%S")
    name_str <- paste('iscan',method, sep='_')
    
    mach_key <- mach_var[1]
    mach_val <- mach_var[2]
    
    time_tib <- date_str %>% stringr::str_split('[\\-\\/ \\:]') %>% BiocGenerics::unlist() %>% 
      purrr::set_names('Year','Mon','Day','Hour','Min','Sec') %>% 
      tibble::enframe() %>% spread(name, value) %>% dplyr::mutate_all(.funs = (as.integer)) %>%
      tibble::add_column('Date' = date_str) %>%
      tibble::add_column(!!mach_key := !!mach_val) %>%
      dplyr::select(!!mach_key, Date, Year, Mon, Day, Hour, Min, Sec) %>%
      purrr::set_names(paste(name_str, names(.), sep='_'))
  })
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done. Date({method})={date_str}.{RET}{RET}"))
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  
  time_tib # %>% gather()
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                    Auto Sample Detection Functions::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
deltaMatrix = function(mat, beg=NULL, end=NULL, minDelta,
                       verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'deltaMatrix'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  stime <- system.time({
    ncols  <- base::ncol(mat)
    nrows  <- base::nrow(mat)
    cnames <- colnames(mat)
    
    stopifnot(nrows>0)
    
    del_mat <- matrix(0, nrow=ncols ,ncol=ncols)
    colnames(del_mat) <- cnames
    rownames(del_mat) <- cnames
    
    if (is.null(beg)) beg <- 1
    if (is.null(end)) end <- ncols
    
    max_cnt <- 0
    max_idx <- 0
    max_per <- 0
    max_key <- 'Unknown'
    for (ii in c(beg:end)) {
      for (jj in c(beg:end)) {
        if (ii<jj) {
          cur_cnt <- length(which(abs(rowDiffs(mat[,c(ii,jj)])) <= minDelta))
          del_mat[ii,jj] <- cur_cnt
          del_mat[jj,ii] <- cur_cnt
        } else if (ii==jj) {
          del_mat[ii,jj] <- length(mat[,ii])
        }
      }
    }
  })
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}{RET}"))
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  
  del_mat
}

rsquaredMatrix = function(mat, verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'rsquaredMatrix'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting.{RET}"))
  
  stime <- system.time({
    cor_val <- NULL
    
    cor_val = tryCatch({
      cor(mat, method="pearson", use="pairwise.complete.obs")
    }, warning = function(w) {
      'warning-SBM-Cor'
    }, error = function(e) {
      'error-SBM-Cor'
    }, finally = {
      'cleanup-SBM-Cor'
    })
  })
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}{RET}"))
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  
  cor_val
}

sampleDetect = function(can, ref, minPval, minDelta, dname, pname, ptype=NULL,
                        join, pval, field, suffix, del='_',
                        outDir=NULL, sname=NULL, plotMatrix=FALSE, writeMatrix=FALSE,
                        dpi=120, format='png', datIdx=4, roundData=TRUE,
                        verbose=0,vt=4,tc=1,tt=NULL) {
  funcTag <- 'sampleDetect'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting.{RET}"))
  
  stime <- system.time({
    sss <- NULL
    
    dname <- dname %>% rlang::sym()
    pname <- pname %>% rlang::sym()
    
    if (!is.null(outDir) && !dir.exists(outDir)) dir.create(outDir, recursive=TRUE)
    
    # 1. Format and Filter Data by Probe Type
    if (!is.null(ptype) && !is.null(pname)) {
      can <- can %>% dplyr::filter(!!pname==ptype)
    }
    
    # 2. Format Candidate Data
    can_nrows <- base::nrow(can)
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} can_nrows={can_nrows}.{RET}"))
    
    tib <- can %>% dplyr::select(join,!!dname,!!pname,pval,field) %>%
      purrr::set_names(join,dname,pname, 'Sample_pval','Sample_beta') %>%
      dplyr::left_join(ref, by=join)
    if (verbose>=vt+2) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Joined.{RET}"))
    
    # 3. Filter on P-value
    tib <- maskTibs(tib, field=field, pval=pval, minPval=minPval, del=del,
                    verbose=verbose,vt=vt+1,tc=tc+1,tt=tTracker)
    
    # 4. Build matrix
    mat <- tib %>% dplyr::select(-c(join,!!dname,!!pname)) %>%
      dplyr::select(ends_with(paste0(del,suffix)) ) %>% as.matrix()
    
    # 5. Plot Matrix
    if (plotMatrix && !is.null(outDir) && !is.null(sname))
      gg <- plotPairs(tib=tib, sample=sname, nameA='Sample', nameB='CanonicalReference',
                      field='beta', field_str='beta', detp='pval', outDir=outDir, minPval=minPval, minDelta=minDelta,
                      dpi=dpi, format=format, datIdx=datIdx,
                      verbose=verbose,vt=vt+1,tc=tc+1,tt=tTracker)
    
    # 5. Build R-Squared Matrix
    r2m <- rsquaredMatrix(mat, verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    r2m_tib <- r2m %>% tibble::as_tibble(rownames='Sample')
    
    # 6. Build Deleta Matrix
    dbm <- deltaMatrix(mat, minDelta=minDelta, verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    dbm_tib <- dbm %>% tibble::as_tibble(rownames='Sample')
    
    cat("Building sss r2m/dbm:\n")
    sss <- dplyr::bind_cols(
      r2m_tib %>% head(n=1) %>% dplyr::select(-Sample) %>% 
        tidyr::gather(key='Sample', value='r2') %>% 
        dplyr::mutate(Sample=stringr::str_remove(Sample, paste0('_',suffix)) ) %>%
        dplyr::filter(Sample!='Sample') %>% 
        dplyr::arrange(-r2) %>% head(n=1) %>%
        purrr::set_names('AutoSample_R2_Key', 'AutoSample_R2_Val'),
      
      dbm_tib %>% head(n=1) %>% dplyr::select(-Sample) %>%
        tidyr::gather(key='Sample', value='passCount') %>%
        dplyr::mutate(Sample=stringr::str_remove(Sample, paste0('_',suffix)),
                      failCount=can_nrows-passCount,
                      passPerc=round(100*passCount/can_nrows,1) ) %>%
        dplyr::filter(Sample!='Sample') %>% 
        dplyr::arrange(-passPerc) %>% head(n=1) %>%
        purrr::set_names('AutoSample_dB_Key', 'AutoSample_Total_Cnt', 'AutoSample_dB_Cnt', 'AutoSample_dB_Val')
    ) %>% dplyr::select(AutoSample_Total_Cnt, AutoSample_R2_Key, AutoSample_R2_Val, 
                        AutoSample_dB_Key, AutoSample_dB_Cnt, AutoSample_dB_Val)
    
    # Write matricies::
    if (writeMatrix && !is.null(outDir) && !is.null(sname)) {
      if (is.null(ptype)) ptype <- 'ALL'
      fname <- paste(sname,'Sample','VS','CanonicalReference',ptype,field,
                     paste0('pval-',minPval),paste0('delta-',minDelta), sep='_')
      r2m_csv <- file.path(outDir, paste0(fname,'.rsquaredMatrix.csv.gz'))
      dbm_csv <- file.path(outDir, paste0(fname,'.deltaMatrix.csv.gz'))
      ord_csv <- file.path(outDir, paste0(fname,'.orderedAutoMatrix.csv.gz'))
      
      if (roundData) r2m_tib <- r2m_tib %>% dplyr::mutate_if(is.numeric, round, 6)
      if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Writing RSquared(CSV)={r2m_csv}.{RET}"))
      readr::write_csv(r2m_tib, r2m_csv)
      if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Writing DeltaBeta(CSV)={dbm_csv}.{RET}"))
      readr::write_csv(dbm_tib, dbm_csv)
      
      # 7. Build r2/delta Table
      ord_tib <- dplyr::inner_join(
        tibble::enframe(r2m[,1], value='r2'), 
        tibble::enframe(dbm[,1], value='delta_PassCount'), by="name") %>%
        dplyr::mutate(delta_PassPerc=round(100*delta_PassCount/can_nrows,1),
                      delta_FailCount=can_nrows-delta_PassCount,
                      name=stringr::str_remove(name, paste0('_',suffix) )
        ) %>%
        dplyr::filter(name!='Sample') %>% dplyr::select(-delta_PassCount)
      
      if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Writing OrderAuto(CSV)={ord_csv}.{RET}"))
      readr::write_csv(ord_tib, ord_csv)
    }
  })
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}{RET}"))
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  
  sss
}

# End of file
