
suppressWarnings(suppressPackageStartupMessages( base::require("tidyverse") ))
suppressWarnings(suppressPackageStartupMessages( base::require("stringr") ))
suppressWarnings(suppressPackageStartupMessages( base::require("glue") ))

# Parallel Computing Packages
suppressWarnings(suppressPackageStartupMessages( base::require("doParallel") ))

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
sesamizeSingleSample = function(prefix, man, add, ref, opt, workflows, retData=FALSE, 
                                del='_', vt=3,tc=0) {
  funcTag <- 'sesamizeSingleSample'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  verbose <- opt$verbosity
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} prefix={prefix}.{RET}"))
  
  # TBD:: Validate all options are present in opt!!!
  
  # Retrun Value::
  ret <- NULL
  
  # tTracker <- timeTracker$new(verbose=verbose)
  tTracker <- timeTracker$new()
  stime <- system.time({
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                              Define Outputs::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    basecode <- basename(prefix)
    out_name <- paste(basecode, opt$platform, opt$manifest, sep=del)
    out_prefix <- file.path(opt$outDir, out_name)
    
    stampBeg_txt <- file.path(opt$outDir, paste(out_name, 'timestamp.beg.txt', sep=del) )
    stampEnd_txt <- file.path(opt$outDir, paste(out_name, 'timestamp.end.txt', sep=del) )
    # system(glue::glue('touch {stampBeg_txt}') )
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                           Extract Raw idat::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    if (!dir.exists(opt$outDir)) dir.create(opt$outDir, recursive=TRUE)
    idat_rds  <- paste(out_prefix, 'idat.rds', sep=del)
    idat_list <- prefixToIdat(prefix=prefix, load=opt$loadIdat, save=opt$saveIdat, rds=idat_rds,
                              verbose=verbose,vt=vt+1,tc=tc+1,tt=tTracker)

    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                     Get Sesame Manifest/Address Tables::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    # ses_man_tib <- man
    # if (opt$subManifest) ses_man_tib <- getSesameManifest(man, idat_list$sig, verbose=verbose,vt=1,tc=tc+1,tt=tTracker)
    ses_man_tib <- getSesameManifest(man, idat_list$sig, verbose=verbose,vt=1,tc=tc+1,tt=tTracker)
    ses_add_tib <- ses_man_tib %>% dplyr::select(Probe_ID) %>% dplyr::inner_join(add, by='Probe_ID')
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                     Initialize Data (Structures) Tables::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                            Build Sample Sheet::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Building Auto-Sample-Sheet: [idat/bead/pheno].{RET}"))
    
    requ_sum_tib <- NULL
    call_sum_tib <- NULL
    sset_sum_tib <- NULL
    phen_sum_tib <- NULL
    pred_sum_tib <- NULL
    
    bead_sum_tib <- ses_add_tib %>% dplyr::filter(Probe_Type=='cg') %>% 
      dplyr::select(Address) %>% dplyr::left_join(idat_list$sig, by="Address") %>% 
      summarise(CG_Bead_Count=n(), CG_Bead_Total=sum(Bead_Grn,Bead_Red, na.rm=TRUE), CG_Bead_AvgRep=CG_Bead_Total/CG_Bead_Count/2)
    pool_sum_tib <- ses_man_tib %>% dplyr::filter(Probe_Type=='cg' | Probe_Type=='ch' | Probe_Type=='rs') %>% 
      dplyr::mutate(Probe_Type=stringr::str_to_upper(Probe_Type)) %>%
      dplyr::group_by(Probe_Type) %>%
      dplyr::summarise(Count=n()) %>% tidyr::spread(Probe_Type,Count) %>% 
      purrr::set_names(paste(names(.),'Manifest_Count',sep='_') ) %>% 
      addBeadPoolToSampleSheet(field='CG_Manifest_Count') %>% dplyr::ungroup()
    
    chipFormat <- idat_list$ann %>% dplyr::select(Chip_Format) %>% dplyr::pull()
    beadPool   <- pool_sum_tib %>% dplyr::select(Bead_Pool) %>% dplyr::pull()
    
    ssheet_tib <- idat_list$ann
    if (!is.null(pool_sum_tib)) ssheet_tib <- ssheet_tib %>% dplyr::bind_cols(pool_sum_tib)
    if (!is.null(bead_sum_tib)) ssheet_tib <- ssheet_tib %>% dplyr::bind_cols(bead_sum_tib)
    # ssheet_tib <- ssheet_tib %>% add_column(minPvalUsed  = opt$minPval)
    # ssheet_tib <- ssheet_tib %>% add_column(minDeltaUsed = opt$minDelta)
    
    if (retData) {
      ret$idat <- idat_list
      ret$sman <- ses_man_tib
      ret$sadd <- ses_add_tib
      ret$bead <- bead_sum_tib
      ret$pool <- pool_sum_tib
      ret$ssheet_tib <- ssheet_tib
    }
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                            Output Files::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    if (opt$buildSubDir) opt$outDir <- file.path(opt$outDir, chipFormat, beadPool)
    if (!dir.exists(opt$outDir)) dir.create(opt$outDir, recursive=TRUE)
    
    basecode   <- idat_list$ann$Sentrix_Name[1]
    out_name   <- paste(basecode, opt$platform, opt$manifest, sep=del)
    out_prefix <- file.path(opt$outDir, out_name)
    
    ssheet_csv <- paste(out_prefix, 'AutoSampleSheet.csv.gz', sep=del)
    calls_csv  <- paste(out_prefix, 'calls.csv.gz', sep=del)
    times_csv  <- paste(out_prefix, 'run-times.csv.gz', sep=del)
    
    raw_sset_csv <- paste(out_prefix, 'raw-sset.csv.gz', sep=del)
    raw_sset_rds <- paste(out_prefix, 'raw-sset.rds', sep=del)
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                          Initialize RAW SSET::
    #
    #  - Required for initial phenotype calcs before channel inference
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    raw_sset <- NULL
    raw_sset <- raw_sset <- initSesameRaw(prefix=prefix, platform=opt$platform, manifest=ses_man_tib,
                                          load=opt$loadSsets, save=opt$saveRawSset, rds=raw_sset_rds,
                                          verbose=verbose,vt=vt+1,tc=tc+1,tt=tTracker)
    raw_call_tib <- sset2calls(sset=raw_sset, workflow='raw', percisionBeta=opt$percisionBeta, percisionPval=opt$percisionPval,
                               verbose=verbose,vt=vt+1,tc=tc+1,tt=tTracker)
    call_sum_tib <- callToSSheet(call=raw_call_tib, idx=0, key='raw', pre=call_sum_tib, 
                                 minNegPval=opt$minNegPval, minOobPval=opt$minOobPval, 
                                 percisionBeta=opt$percisionBeta, percisionPval=opt$percisionPval,
                                 verbose=verbose,vt=vt+1,tc=tc+1,tt=tTracker)
    phen_sum_tib <- ssetToInferences(sset=raw_sset, idx=0, key='raw', pre=phen_sum_tib, verbose=verbose,vt=vt+1,tc=tc+1,tt=tTracker)
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                 SSET to Calls by Order of Operations:: workflows
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    all_call_tib <- NULL
    auto_beta_key <- NULL
    auto_negs_key <- NULL
    workflow_cnt <- length(workflows)
    for (ww in seq(1,workflow_cnt)) {
      cur_workflow <- workflows[ww]
      cur_sset_rds <- paste(out_prefix, paste0(cur_workflow,'-sset.rds'), sep=del)
      cur_sset_csv <- paste(out_prefix, paste0(cur_workflow,'-sset.csv.gz'), sep=del)
      cur_ssum_csv <- paste(out_prefix, paste0(cur_workflow,'-ssum.csv.gz'), sep=del)
      
      if (is.null(auto_beta_key)) auto_beta_key <- paste(cur_workflow,'beta', sep=del)
      if (is.null(auto_negs_key)) auto_negs_key <- paste(cur_workflow,'negs', sep=del)
      
      cur_sset <- NULL
      if (opt$loadSsets && file.exists(cur_sset_rds)) {
        if (verbose>=vt+1) cat(glue::glue("[{funcTag}]:{tabsStr} Loading RDS={cur_sset_rds}.{RET}"))
        cur_sset <- readr::read_rds(cur_sset_rds)
      } else {
        stopifnot(!is.null(raw_sset))
        cur_sset <- mutateSSET_workflow(sset=raw_sset, workflow=cur_workflow, save=opt$saveSsets, rds=cur_sset_rds,
                                        verbose=verbose,vt=vt+1,tc=tc+1,tt=tTracker)
      }
      stopifnot(!is.null(cur_sset))
      
      cur_sset_tib <- sset2tib(sset=cur_sset, by="Probe_ID", des="Probe_Design",  
                               percision=opt$percisionSigs, sort=FALSE, save=opt$writeSset, csv=cur_sset_csv, 
                               verbose=verbose,vt=vt+1,tc=tc+1,tt=tTracker)
      sset_dat_tib <- sigTibToSSheet(sigs=cur_sset_tib, man=ses_man_tib, save=opt$writeSsum, csv=cur_ssum_csv,
                                     verbose=verbose,vt=vt+1,tc=tc+1,tt=tTracker)
      sset_sum_tib <- sigsSumToSSheet2(tib=sset_dat_tib, metric='avg', 
                                       verbose=verbose,vt=vt+1,tc=tc+1,tt=tTracker)
      cur_call_tib <- sset2calls(sset=cur_sset, workflow=cur_workflow, 
                                 percisionBeta=opt$percisionBeta, percisionPval=opt$percisionPval,
                                 verbose=verbose,vt=vt+1,tc=tc+1,tt=tTracker)
      call_sum_tib <- callToSSheet(call=cur_call_tib, idx=ww, key=cur_workflow, pre=call_sum_tib, 
                                   minNegPval=opt$minNegPval, minOobPval=opt$minOobPval, 
                                   percisionBeta=opt$percisionBeta, percisionPval=opt$percisionPval,
                                   verbose=verbose,vt=vt+1,tc=tc+1,tt=tTracker)
      phen_sum_tib <- ssetToInferences(sset=cur_sset, idx=ww, key=cur_workflow, pre=phen_sum_tib, 
                                       verbose=verbose,vt=vt+1,tc=tc+1,tt=tTracker)
      
      all_call_tib <- joinTibbles(all_call_tib, cur_call_tib, by='Probe_ID', side='full')
      
      if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} cur_workflow={cur_workflow}. Next...{RET}{RET}"))
    }
    readr::write_csv(all_call_tib, calls_csv)
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                       Add Summarize Phenotype Inferences::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    ssheet_tib <- ssheet_tib %>% dplyr::bind_cols(phen_sum_tib)
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                  Add Summarize Call Stats (p-vales/beta values)::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    ssheet_tib <- ssheet_tib %>% 
      add_column(MinNegPval=opt$minNegPval) %>%
      add_column(MinOobPval=opt$minOobPval) %>%
      dplyr::bind_cols(call_sum_tib)
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                             Auto-Detect Sample::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    if (opt$autoDetect) {
      if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Detecting Auto-Sample-Sheet: [SampleName, rsquared,deltaBeta].{RET}"))
      
      # Prep all_call_tib to match previous format::
      #  Probe_ID       M        U DESIGN COLOR_CHANNEL col   Probe_Type Probe_Source Next_Base
      auto_call <- ses_man_tib %>% dplyr::rename(Design_Type=DESIGN) %>% 
        dplyr::select(Probe_ID, Design_Type, Probe_Type) %>%
        dplyr::left_join(all_call_tib, by="Probe_ID") # %>% dplyr::select(!!auto_negs_key, !!auto_beta_key)
      
      auto_ref <- ref
      if (grep("_beta", names(ref)) %>% length() == 0) {
        auto_ref <- ref %>% purrr::set_names(paste(names(.),'beta', sep=del) )
        names(auto_ref)[1] <- 'Probe_ID'
      }
      
      auto_data <- sampleDetect(can=auto_call, ref=auto_ref, minPval=opt$minNegPval, minDelta=opt$minDeltaBeta,
                                dname='Design_Type', pname='Probe_Type', ptype='cg',
                                jval='Probe_ID', field=auto_beta_key, pval=auto_negs_key, suffix='beta', del=del,
                                outDir=opt$outDir, sname=out_name, plotMatrix=opt$plotAuto, writeMatrix=opt$writeAuto,
                                dpi=opt$dpi, format=opt$plotFormat, datIdx=4,
                                verbose=verbose,vt=vt+1,tc=tc+1,tt=tTracker)

      ssheet_tib <- ssheet_tib %>% dplyr::bind_cols(auto_data)
    }

    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                       Add Summarize Sample Intensities::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    ssheet_tib <- ssheet_tib %>% dplyr::bind_cols(sset_sum_tib)
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                             Add Requeue Flags::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    requ_sum_tib <- tibble::tibble(
      Requeue_Flag_Oob=case_when(ssheet_tib$Poob_Pass_0_Perc < opt$minOobPerc ~ FALSE, TRUE ~ TRUE),
      Requeue_Flag_Neg=case_when(ssheet_tib$Negs_Pass_0_Perc < opt$minNegPerc ~ FALSE, TRUE ~ TRUE),
      Requeue_Pass_Perc_Oob=opt$minOobPerc,
      Requeue_Pass_Perc_Neg=opt$minNegPerc
    )
    ssheet_tib <- ssheet_tib %>% dplyr::bind_cols(requ_sum_tib) %>% 
      dplyr::select(Requeue_Flag_Oob, Requeue_Flag_Neg, Requeue_Pass_Perc_Oob, Requeue_Pass_Perc_Neg, everything())
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                  Add to Full Return Values:: R&D Testing ONLY::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    # 
    if (retData) {
      ret$sset_dat <- sset_dat_tib
      ret$sset_sum <- sset_sum_tib
      ret$raw_sset <- raw_sset
      ret$cur_sset <- cur_sset
      ret$sset_tib <- cur_sset_tib
      ret$all_call <- all_call_tib
    }
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                               Write Outputs::
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Writing Outputs: [call,sigs,signalSummary,sampleSheet,time].{RET}"))
    
    ssheet_tib <- ssheet_tib %>% dplyr::select(-starts_with('Iscan_'), starts_with('Iscan_')) %>% 
      dplyr::mutate_if(is.numeric, list(round), 4)
    if (opt$writeSsheet && !is.null(ssheet_tib)) readr::write_csv(ssheet_tib, ssheet_csv)
    
    time_tib <- tTracker$time
    if (opt$percisionBeta!=0) time_tib <- time_tib %>% dplyr::mutate_if(base::is.numeric, round, opt$percisionBeta)
    if (!is.null(time_tib)) readr::write_csv(time_tib, times_csv)
    
    # system(glue::glue('touch {stampEnd_txt}') )
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (verbose>=vt) tTracker %>% print()
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done. elasped={etime}.{RET}{RET}"))
  
  if (retData) {
    ret$ssheet_tib <- ssheet_tib
    ret$time       <- tTracker
    return(ret)
  }
  
  ssheet_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                          Call (Beta/Pval) Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

maskTibs = function(tib, field, pval, minPval, del='_', 
                    rm.suffix=FALSE, only.field=FALSE, no.mask=FALSE,
                    verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'maskTibs'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting; field={field}, pval={pval}, minPval={minPval}.{RET}"))
  if (!no.mask) {
    pval_del <- paste0(del,pval)
    snames <- tib %>% dplyr::select(ends_with(pval_del)) %>% names() %>% stringr::str_remove(pval_del)
    for (sname in snames) {
      fstr <- paste(sname,field, sep=del)
      pstr <- paste(sname,pval, sep=del)
      tib <- maskTib(tib, field=fstr, pval=pstr, minPval=minPval, verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    }
  }
  if (only.field) tib <- tib %>% dplyr::select(1,ends_with(field))
  if (rm.suffix) tib  <- tib %>% purrr::set_names(stringr::str_remove(names(.), paste0(del,field) ) )
  
  tib
}

maskTib = function(tib, field, pval, minPval, rm.pval=FALSE, only.field=FALSE,
                   verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'maskTib'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting; field={field}, pval={pval}, minPval={minPval}.{RET}"))
  
  field <- field %>% rlang::sym()  
  pval  <- pval %>% rlang::sym()
  # minPval <- minPval %>% rlang::sym()
  
  stime <- system.time({
    tib <- tib %>% dplyr::mutate(!!field := case_when(!!pval >= !!minPval ~ NA_real_, TRUE ~ !!field))
    if (rm.pval) tib <- tib %>% dplyr::select(-!!pval)
    if (only.field) tib <- tib %>% dplyr::select(!!field)
  })
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}{RET}"))
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  
  tib
}

maskCall = function(tib, field, minKey, minVal, verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'maskCall'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Masking {field} with {minKey} at {minVal}.{RET}"))
  
  field <- field %>% rlang::sym()  
  minKey  <- minKey %>% rlang::sym()
  # minVal <- minVal %>% rlang::sym()
  
  stime <- system.time({
    tib <- tib %>% dplyr::mutate(!!field := case_when(!!minKey >= !!minVal ~ NA_real_, TRUE ~ !!field))
  })
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}{RET}"))
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  
  tib
}

loadCallFile = function(file, selKey, datKey=NULL, minKey=NULL, minVal=NULL, prefix=NULL, 
                        retMin=FALSE, retRaw=FALSE, del='_',
                        verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'loadCallFile'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) {
    if (!is.null(datKey)) {
      cat(glue::glue("[{funcTag}]:{tabsStr} selKey={selKey}, datKey={datKey} file={file}.{RET}")) 
    } else { 
      cat(glue::glue("[{funcTag}]:{tabsStr} selKey={selKey}, datKey=LOAD_ALL, file={file}.{RET}"))
    }
  }
  
  stime <- system.time({
    selKey <- selKey %>% rlang::sym()
    if (!is.null(datKey)) datKey <- datKey %>% rlang::sym()
    if (!is.null(minKey)) minKey <- minKey %>% rlang::sym()
    
    if (stringr::str_ends(file,'.rds')) {
      tib <- suppressMessages(suppressWarnings(readr::read_rds(file)) )
    } else {
      tib <- suppressMessages(suppressWarnings(readr::read_csv(file)) )
    }
    if (verbose>=vt+4) print(tib)
    
    if (retRaw && !is.null(datKey)) {
      if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr} IS_RAW retRaw={retRaw}, datKey={datKey}.{RET}"))
      tib <- tib %>% dplyr::select(!!selKey,ends_with(paste0(del,datKey)) )
      if (verbose>=vt+4) print(tib)
    } else {
      if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr} NOT_RAW retRaw={retRaw}, datKey={datKey}.{RET}"))
      
      if (!is.null(minVal) && !is.null(datKey)) {
        if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr} USE_PVAL minVal={minVal}, SELECT=[selKey={selKey}, minKey={minKey}, datKey={datKey}].{RET}"))
        tib <- tib %>% dplyr::select(!!selKey,!!minKey, !!datKey)
        if (verbose>=vt+4) print(tib)
        
        tot_cnt  <- tib %>% base::nrow()
        pre_cnt  <- tib %>% dplyr::filter(is.na(!!datKey)) %>% base::nrow()
        
        tib <- maskTib(tib=tib, field=datKey, pval=minKey, minPval=minVal, verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
        # tib <- maskCall(tib=tib, field=datKey, minKey=minKey, minVal=minVal, verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
        
        pos_cnt  <- tib %>% dplyr::filter(is.na(!!datKey)) %>% base::nrow()
        pos_per  <- round(100*pos_cnt/tot_cnt,3)
        
        if (retMin) tib <- tib %>% dplyr::select(!!selKey,!!datKey)
        if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Masked ({pos_per}%) Beta Values {pre_cnt} -> {pos_cnt}/{tot_cnt}.{RET}"))
      } else {
        if (verbose>=vt+4) cat(glue::glue("[{funcTag}]:{tabsStr} NO_PVAL.{RET}"))
        if (!is.null(datKey)) tib <- tib %>% dplyr::select(!!selKey,!!datKey)
        if (verbose>=vt+4) print(tib)
      }
      if (!is.null(prefix)) tib <- tib %>% purrr::set_names(paste(prefix,names(tib), sep='.') ) %>% dplyr::rename(!!selKey := 1)
    }
  })
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}{RET}"))
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  
  tib
}

loadCallFiles = function(files, selKey, datKey=NULL, minKey=NULL, minVal=NULL, prefix=NULL, 
                         addRep=TRUE, retMin=FALSE, retRaw=FALSE, max=NULL, del='_',
                         verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'loadCallFiles'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  stopifnot(is.vector(files))
  
  join_tibs <- NULL
  files_cnt <- length(files)
  if (verbose>=vt) {
    cat(glue::glue("[{funcTag}]:{tabsStr} Total number of call files={files_cnt}, selKey={selKey}, datKey={datKey}.{RET}"))
  } else {
    cat(glue::glue("[{funcTag}]:{tabsStr} Total number of call files={files_cnt}, selKey={selKey}, datKey=ALL_DATA.{RET}"))
  }
  
  stime <- system.time({
    for (ii in c(1:files_cnt)) {
      if (verbose>=vt+2) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} ii={ii}, file={files[ii]}.{RET}"))
      tib <- loadCallFile(file=files[ii], selKey=selKey, datKey=datKey, minKey=minKey, minVal=minVal,
                          retMin=retMin, retRaw=retRaw,
                          verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
      rep_str <- paste0('Rep',ii)
      if (!is.null(prefix)) rep_str <- paste(prefix,rep_str, sep=del)
      if (addRep) tib <- tib %>% purrr::set_names(paste(rep_str,names(tib), sep='.') ) %>% dplyr::rename(!!selKey := 1)
      
      if (is.null(join_tibs)) { join_tibs <- tib
      } else { join_tibs <- join_tibs %>% dplyr::full_join(tib, by=selKey) }
      
      if (!is.null(max) && ii >= max) break
    }
  })
  nrows <- join_tibs %>% base::nrow()
  ncols <- join_tibs %>% base::ncol()
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} nrows={nrows}, ncols={ncols}.{RET}{RET}"))
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  
  join_tibs
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
                        jval, pval, field, suffix, del='_',
                        outDir=NULL, sname=NULL, plotMatrix=FALSE, writeMatrix=FALSE,
                        dpi=120, format='png', datIdx=4, roundData=TRUE,
                        verbose=0,vt=4,tc=1,tt=NULL) {
  funcTag <- 'sampleDetect'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting. dname={dname}, pname={pname},",
                                  "field={field}, jval={jval}, pval={pval}.{RET}"))
  
  if (verbose>=vt+2) cat(glue::glue("[{funcTag}]:{tabsStr} CAN::{RET}") )
  if (verbose>=vt+2) print(can)
  if (verbose>=vt+2) cat(glue::glue("[{funcTag}]:{tabsStr} REF::{RET}") )
  if (verbose>=vt+2) print(ref)
  
  stime <- system.time({
    sss <- NULL
    can_names <- c(jval,dname,pname, 'Sample_pval','Sample_beta')
    
    # String Values::
    jval_char <- jval
    
    dname <- dname %>% rlang::sym()
    pname <- pname %>% rlang::sym()
    
    # Sym jval to see if it fixes the warning...
    field <- field %>% rlang::sym()
    jval  <- jval %>% rlang::sym()
    pval  <- pval %>% rlang::sym()
    
    if (!is.null(outDir) && !dir.exists(outDir)) dir.create(outDir, recursive=TRUE)
    
    # 1. Format and Filter Data by Probe Type
    if (!is.null(ptype) && !is.null(pname)) {
      can <- can %>% dplyr::filter(!!pname==ptype)
    }
    
    # 2. Format Candidate Data
    can_nrows <- base::nrow(can)
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} can_nrows={can_nrows}.{RET}"))
    
    if (verbose>=vt+2) {
      cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Pre-join. dname={dname}, pname={pname}, pval={pval}, field={field}, jval={jval}.{RET}"))
      cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} ref=={RET}"))
      print(ref)
      cat("\n\ntib1==\n")
      can %>% dplyr::select(!!jval,!!dname,!!pname,!!pval,!!field) %>% print()
      cat("\n\ntib2==\n")
      can %>% dplyr::select(!!jval,!!dname,!!pname,!!pval,!!field) %>%
        purrr::set_names(can_names) %>% print()
      cat("\n\nVerbose Done!\n\n")
    }
    tib <- can %>% dplyr::select(!!jval,!!dname,!!pname,!!pval,!!field) %>%
      purrr::set_names(can_names) %>%
      dplyr::left_join(ref, by=jval_char)
    if (verbose>=vt+2) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Joined.{RET}"))
    
    # 3. Filter on P-value
    tib <- maskTibs(tib, field=field, pval=pval, minPval=minPval, del=del,
                    verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    
    # 4. Build matrix
    mat <- tib %>% dplyr::select(-c(!!jval,!!dname,!!pname)) %>%
      dplyr::select(ends_with(paste0(del,suffix)) ) %>% as.matrix()
    
    # 5. Plot Matrix
    if (plotMatrix && !is.null(outDir) && !is.null(sname))
      gg <- plotPairs(tib=tib, sample=sname, nameA='Sample', nameB='CanonicalReference',
                      field='beta', field_str='beta', detp='pval', outDir=outDir, minPval=minPval, minDelta=minDelta,
                      dpi=dpi, format=format, datIdx=datIdx,
                      verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    
    # 5. Build R-Squared Matrix
    r2m <- rsquaredMatrix(mat, verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    r2m_tib <- r2m %>% tibble::as_tibble(rownames='Sample')
    
    # 6. Build Deleta Matrix
    dbm <- deltaMatrix(mat, minDelta=minDelta, verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    dbm_tib <- dbm %>% tibble::as_tibble(rownames='Sample')
    
    # cat("Building sss r2m/dbm:\n")
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
      
      if (roundData) r2m_tib <- r2m_tib %>% dplyr::mutate_if(base::is.numeric, round, 6)
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
