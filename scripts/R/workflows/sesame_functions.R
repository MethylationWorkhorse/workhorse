
suppressWarnings(suppressPackageStartupMessages( base::require("tidyverse") ))
suppressWarnings(suppressPackageStartupMessages( base::require("stringr") ))
suppressWarnings(suppressPackageStartupMessages( base::require("glue") ))

suppressWarnings(suppressPackageStartupMessages( base::require("sesame") ))

COM <- ","
TAB <- "\t"
RET <- "\n"

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#      Predefined SSET to Calls by Order of Operations Workflows::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                         Sesame Mutation Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
initSesameRaw = function(prefix, platform, manifest, load=FALSE, save=FALSE, rds=NULL, 
                         verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'initSesameRaw'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting; platform={platform}.{RET}"))
  stime <- system.time({
    if (load && !is.null(rds) && file.exists(rds)) {
      if (verbose>=vt+1) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Loading RDS={rds}.{RET}"))
      sset <- readr::read_rds(rds)
    } else {
      sset <- sesame::readIDATpair(prefix, platform=platform, manifest=manifest)
      if (save && !is.null(rds)) {
        if (verbose>=vt+1) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Writing RDS={rds}.{RET}"))
        readr::write_rds(sset, rds, compress="gz")
      }
    }
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done; elapsed={etime}.{RET}{RET}"))
  
  sset
}

# Predefined SSET to Calls by Order of Operations Workflows::
mutateSSET_workflow = function(sset, workflow, save=FALSE, rds=NULL, 
                               verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'mutateSSET'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting; workflow={workflow}.{RET}"))
  
  stime <- system.time({
    if (workflow=='din') {
      sset <- sset %>% 
        mutateSesame(method = 'dyeBiasCorrTypeINorm', verbose=opt$verbosity,vt=vt+1,tc=tc+1,tt=tt) %>%
        mutateSesame(method = 'inferTypeIChannel', verbose=opt$verbosity,vt=vt+1,tc=tc+1,tt=tt) %>%
        mutateSesame(method = 'noob', verbose=opt$verbosity,vt=vt+1,tc=tc+1,tt=tt)
    } else if (workflow=='ind') {
      sset <- sset %>% 
        mutateSesame(method = 'inferTypeIChannel', verbose=opt$verbosity,vt=vt+1,tc=tc+1,tt=tt) %>%
        mutateSesame(method = 'noob', verbose=opt$verbosity,vt=vt+1,tc=tc+1,tt=tt) %>%
        mutateSesame(method = 'dyeBiasCorrTypeINorm', verbose=opt$verbosity,vt=vt+1,tc=tc+1,tt=tt)
    } else if (workflow=='ndi') {
      sset <- sset %>% 
        mutateSesame(method = 'noob', verbose=opt$verbosity,vt=vt+1,tc=tc+1,tt=tt) %>%
        mutateSesame(method = 'dyeBiasCorrTypeINorm', verbose=opt$verbosity,vt=vt+1,tc=tc+1,tt=tt) %>%
        mutateSesame(method = 'inferTypeIChannel', verbose=opt$verbosity,vt=vt+1,tc=tc+1,tt=tt)
    } else {
      stop(glue::glue("[{funcTag}]: ERROR: Unsupported workflow={workflow}!{RET}{RET}"))
    }
    
    if (save && !is.null(rds)) {
      if (verbose>=vt+1) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Writing RDS={rds}.{RET}"))
      readr::write_rds(sset, rds, compress="gz")
    }
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done; elapsed={etime}.{RET}{RET}"))
  
  return(sset)
}

mutateSesame = function(sset, method, verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'mutateSesame'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Mutating Sesame({method}){RET}"))
  stime <- system.time({
    if (is.null(method)) stop(glue::glue("{RET}[{funcTag}]: ERROR: Missing method!!!{RET}{RET}"))
    else if (method=='open') sset <- sset %>% sesame::pOOBAH() %>% sesame::noob() %>% sesame::dyeBiasCorrTypeINorm()
    else if (method=='dyeBiasCorrTypeINorm') sset <- sset %>% sesame::dyeBiasCorrTypeINorm()
    else if (method=='detectionPnegEcdf') sset <- sset %>% sesame::detectionPnegEcdf()
    else if (method=='pOOBAH') sset <- sset %>% sesame::pOOBAH()
    else if (method=='noob') sset <- sset %>% sesame::noob()
    else if (method=='noobsb') sset <- sset %>% sesame::noobsb()
    else if (method=='inferTypeIChannel') sset <- sset %>% sesame::inferTypeIChannel(switch_failed=TRUE, verbose=FALSE)
    else if (method=='raw') { } # sset <- sset
    else stop(glue::glue("{RET}[{funcTag}]: ERROR: Unsupported method={method}!!!{RET}{RET}"))
  })
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}{RET}"))
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  
  sset
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                    Sesame SSET Conversion Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

sigsSumToSSheet2 = function(tib, metric='avg',
                            verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'sigsSumToSSheet2'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  pt_tib <- tib %>% split(.$Probe_Type)
  
  tib_names <- pt_tib %>% names()
  core_pts <- c('cg', 'ch', 'rs')
  
  ss_tib <- NULL
  for (pt in tib_names) {
    if (!c(pt) %in% core_pts) next
    
    nrows = base::nrow(pt_tib[[pt]])
    for (ii in c(1:nrows)) {
      pt_val <- pt_tib[[pt]]$Probe_Type[ii]
      pd_val <- pt_tib[[pt]]$Probe_Design[ii]
      
      key    <- paste(pt_val, pd_val, sep='_')
      valM   <- pt_tib[[pt]]$M_avg[ii];
      valU   <- pt_tib[[pt]]$U_avg[ii];
      
      if (pd_val=='II') {
        keyM <- paste(key,'G',metric, sep='_') %>% rlang::sym()
        keyU <- paste(key,'R',metric, sep='_') %>% rlang::sym()
      } else {
        keyM <- paste(key,'M',metric, sep='_') %>% rlang::sym()
        keyU <- paste(key,'U',metric, sep='_') %>% rlang::sym()
      }
      new_tib <- tibble::tibble(!!keyM := valM, !!keyU := valU) %>%
        dplyr::mutate_if(is.numeric, list(round), 0)
      
      ss_tib <- ss_tib %>% bind_cols(new_tib)
    }
    # cat(glue::glue("1st: pt={pt}.{RET}"))
  }
  # ss_tib %>% as.data.frame() %>% print()
  
  # ss_tib <- NULL
  for (pt in tib_names) {
    if (c(pt) %in% core_pts) next
    
    nrows = base::nrow(pt_tib[[pt]])
    for (ii in c(1:nrows)) {
      pt_val <- pt_tib[[pt]]$Probe_Type[ii]
      pd_val <- pt_tib[[pt]]$Probe_Design[ii]
      
      # key    <- paste(pt_val, pd_val, sep='_')
      key    <- paste(pt_val, sep='_')
      valM   <- pt_tib[[pt]]$M_avg[ii];
      valU   <- pt_tib[[pt]]$U_avg[ii];
      
      if (pd_val=='II') {
        keyM <- paste(key,'G',metric, sep='_') %>% rlang::sym()
        keyU <- paste(key,'R',metric, sep='_') %>% rlang::sym()
      } else {
        keyM <- paste(key,'M',metric, sep='_') %>% rlang::sym()
        keyU <- paste(key,'U',metric, sep='_') %>% rlang::sym()
      }
      new_tib <- tibble::tibble(!!keyM := valM, !!keyU := valU) %>%
        dplyr::mutate_if(is.numeric, list(round), 0)
      
      ss_tib <- ss_tib %>% bind_cols(new_tib)
    }
    # cat(glue::glue("2nd: pt={pt}.{RET}"))
  }
  ss_tib <- ss_tib %>% dplyr::mutate_if(is.numeric, as.integer)
  if (verbose>=vt) ss_tib %>% as.data.frame() %>% print()
  
  ss_tib
}

sigTibToSSheet = function(sigs, man=NULL, by="Probe_ID", type="Probe_Type", des="Probe_Design", 
                          percision=-1, sort=FALSE, save=FALSE, csv=NULL,
                          verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'sigTibToSSheet'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  tib <- NULL
  stime <- system.time({
    des <- des %>% rlang::sym()
    
    # tib <- sigs %>% dplyr::select(-all_of(by)) %>% dplyr::group_by(Probe_Design) %>% summarise_if(is.numeric, list(mu=mean), na.rm=TRUE)
    if (!is.null(man)) {
      type <- type %>% rlang::sym()
      sigs <- dplyr::select(man, !!by, !!type) %>% dplyr::right_join(sigs, by=by ) %>% dplyr::group_by(!!type, !!des) 
      tib <- sigs %>% summarise_if(is.numeric, list(min=min, med=median, avg=mean, sd=sd, max=max), na.rm=TRUE)
    } else {
      stop(glue::glue("\n[{func}]: ERROR: Currently only supported with the manifest!!!{RET}{RET}"))
      tib <- sigs %>% dplyr::group_by(!!des) %>% summarise_if(is.numeric, list(mu=mean), na.rm=TRUE)
    }
    
    if (save && !is.null(csv)) {
      if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Writing signal set (percision={percision}) CSV={csv}.{RET}"))
      readr::write_csv(tib, csv)
    }
  })
  nrows <- tib %>% base::nrow()
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done nrows={nrows}, elapsed={etime}.{RET}{RET}"))
  
  tib
}

sset2tib = function(sset, man=NULL, by="Probe_ID", type="Probe_Type", des="Probe_Design", 
                    percision=-1, sort=FALSE, save=FALSE, csv=NULL,
                    verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'sset2tib'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  tib <- NULL
  stime <- system.time({
    
    des <- des %>% rlang::sym()
    
    tib <- dplyr::bind_rows(
      sset@IG   %>% tibble::as_tibble(rownames='Probe_ID',.name_repair = "unique") %>% dplyr::mutate(!!des:='IG'),
      sset@oobG %>% tibble::as_tibble(rownames='Probe_ID',.name_repair = "unique") %>% dplyr::mutate(!!des:='OG'),
      
      sset@IR   %>% tibble::as_tibble(rownames='Probe_ID',.name_repair = "unique") %>% dplyr::mutate(!!des:='IR'),
      sset@oobR %>% tibble::as_tibble(rownames='Probe_ID',.name_repair = "unique") %>% dplyr::mutate(!!des:='OR'),
      
      sset@II   %>% tibble::as_tibble(rownames='Probe_ID',.name_repair = "unique") %>% dplyr::mutate(!!des:='II')
    )
    tib <- tib %>% dplyr::select(!!by, !!des, everything())
    
    if (percision!=-1) tib <- tib %>% dplyr::mutate_if(is.numeric, list(round), percision)
    if (!is.null(man)) tib <- dplyr::select(man, !!by, !!type) %>% dplyr::right_join(tib, by=by )
    
    by  <- by %>% rlang::sym()
    if (sort) tib <- tib %>% dplyr::arrange(!!by)
    
    if (save && !is.null(csv)) {
      if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Writing signal set (percision={percision}) CSV={csv}.{RET}"))
      readr::write_csv(tib, csv)
    }
  })
  nrows <- tib %>% base::nrow()
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done nrows={nrows}, elapsed={etime}.{RET}{RET}"))
  
  tib
}

sset2calls = function(sset, workflow, percisionBeta=0, percisionPval=0,
                      verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'sset2calls'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting; workflow={workflow}.{RET}"))
  
  tib <- NULL
  stime <- system.time({
    name <- paste(workflow,'beta', sep='_')
    beta <- ssetToBetaTib(sset=sset, name=name, as.enframe=FALSE, percision=percisionBeta, 
                          verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    tib <- tibble::enframe(beta, name='Probe_ID', value=name)
    
    name <- paste(workflow,'negs', sep='_')
    sset <- mutateSesame(sset=sset, method='detectionPnegEcdf', verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    pval <- ssetToPvalTib(sset=sset, name=name, percision=percisionPval, verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    tib <- tib %>% dplyr::left_join(pval, by="Probe_ID")
    
    name <- paste(workflow,'poob', sep='_')
    sset <- mutateSesame(sset=sset, method='pOOBAH', verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    pval <- ssetToPvalTib(sset=sset, name=name, percision=percisionPval, verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    tib <- tib %>% dplyr::left_join(pval, by="Probe_ID")
  })
  nrows <- tib %>% base::nrow()
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done nrows={nrows}, elapsed={etime}.{RET}{RET}"))
  
  tib
}

sset2bset = function(sset, addSigs=TRUE, addNegs=TRUE, addPoob=TRUE, addBeta=TRUE,
                     
                     # getBeta Parameters::
                     quality.mask=FALSE, nondetection.mask=FALSE, 
                     mask.use.tcga=FALSE, pval.threshold=1, sum.TypeI=FALSE,
                     as.enframe=TRUE,
                     round_dat=TRUE, round_pval=6, round_beta=4,
                     verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'sset2bset'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  bset <- NULL
  stime <- system.time({
    sigs <- dplyr::bind_rows(
      sset@IG   %>% tibble::as_tibble(rownames='Probe_ID',.name_repair = "unique") %>% dplyr::mutate(col='IG'),
      sset@oobG %>% tibble::as_tibble(rownames='Probe_ID',.name_repair = "unique") %>% dplyr::mutate(col='OG'),
      
      sset@IR   %>% tibble::as_tibble(rownames='Probe_ID',.name_repair = "unique") %>% dplyr::mutate(col='IR'),
      sset@oobR %>% tibble::as_tibble(rownames='Probe_ID',.name_repair = "unique") %>% dplyr::mutate(col='OR'),
      
      sset@II   %>% tibble::as_tibble(rownames='Probe_ID',.name_repair = "unique") %>% dplyr::mutate(col=NA)
    )
    if (round_dat) sigs <- sigs %>% dplyr::mutate_if(is.numeric, list(as.integer))
    
    if (addSigs) { bset <- sigs
    } else { bset <- dplyr::select(sigs, 'Probe_ID', 'col') }
    
    if (addNegs) {
      negs <- sset %>% sesame::detectionPnegEcdf() %>% sesame::pval() %>% tibble::enframe(name='Probe_ID', value='negs')
      if (round_dat && round_pval>0) negs <- negs %>% dplyr::mutate_if(is.numeric, list(round), round_pval)
      bset <- bset %>% dplyr::left_join(negs, by='Probe_ID')
    }
    
    if (addPoob) {
      poob <- sset %>% sesame::pOOBAH() %>% sesame::pval() %>% tibble::enframe(name='Probe_ID', value='poob')
      if (round_dat && round_pval>0) poob <- poob %>% dplyr::mutate_if(is.numeric, list(round), round_pval)
      bset <- bset %>% dplyr::left_join(poob, by='Probe_ID')
    }
    
    if (addBeta) {
      beta <- sesame::getBetas(sset, quality.mask=quality.mask, nondetection.mask=nondetection.mask, 
                               mask.use.tcga=mask.use.tcga, pval.threshold=pval.threshold, sum.TypeI=sum.TypeI)
      if (as.enframe) beta <- beta %>% tibble::enframe(name='Probe_ID', value='beta')
      if (round_dat && round_beta>0) beta <- beta %>% dplyr::mutate_if(is.numeric, list(round), round_beta)
      bset <- bset %>% dplyr::left_join(beta, by='Probe_ID')
    }
    
    # bset <- bset %>% 
    #   dplyr::left_join(negs, by='Probe_ID') %>%
    #   dplyr::left_join(poob, by='Probe_ID') %>%
    #   dplyr::left_join(beta, by='Probe_ID')
  })
  bset_nrows <- bset %>% base::nrow()
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done bset_nrows={bset_nrows}, elapsed={etime}.{RET}{RET}"))
  
  bset
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                      Sesame SSET To Inference Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# TBD:: If manifest is provided then seperate by Probe_Type and Probe_Design
callToSSheet = function(call, idx, key, pre=NULL, minNegPval, minOobPval, 
                        percisionBeta=0, percisionPval=0, del='_', onlyCG=TRUE, id='Probe_ID',
                        verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'callToSSheet'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  tib <- NULL
  stime <- system.time({

    id <- id %>% rlang::sym()
    if (onlyCG) call <- call %>% dplyr::filter(stringr::str_starts(!!id, 'cg'))

    beta_tib <- NULL
    inf_key  <- paste('Beta',idx,'Method', sep=del) %>% rlang::sym()
    inf_val  <- paste('Beta',idx,'Mean', sep=del) %>% rlang::sym()
    beta_tib <- tibble::tibble(!!inf_key := key, !!inf_val := dplyr::select(call, ends_with('_beta')) %>% 
                                 dplyr::summarise_all(list(mean=mean), na.rm=TRUE) %>% 
                                 dplyr::pull() %>% round(percisionBeta) )
    
    negs_tib <- NULL
    inf_key  <- paste('Negs_Pass',idx,'Method', sep=del) %>% rlang::sym()
    inf_val  <- paste('Negs_Pass',idx,'Perc', sep=del) %>% rlang::sym()
    negs_tib <- tibble::tibble(!!inf_key := key, !!inf_val := dplyr::select(call, ends_with('_negs')) %>% 
                                 dplyr::summarise_all(list(pass_perc=cntPer_lte), min=minNegPval) %>%
                                 dplyr::pull() %>% round(percisionPval) )
    
    poob_tib <- NULL
    inf_key  <- paste('Poob_Pass',idx,'Method', sep=del) %>% rlang::sym()
    inf_val  <- paste('Poob_Pass',idx,'Perc', sep=del) %>% rlang::sym()
    poob_tib <- tibble::tibble(!!inf_key := key, !!inf_val := dplyr::select(call, ends_with('_poob')) %>% 
                                 dplyr::summarise_all(list(pass_perc=cntPer_lte), min=minOobPval) %>%
                                 dplyr::pull() %>% round(percisionPval) )
    
    tib <- dplyr::bind_cols(pre, beta_tib, negs_tib, poob_tib)
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done; elapsed={etime}.{RET}{RET}"))
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  
  tib
}

ssetToInferences = function(sset, idx, key, pre=NULL, del='_', verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'ssetToInferences'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  tib <- NULL
  stime <- system.time({
    
    gct_tib <- NULL
    inf_key <- paste('GCT',idx,'Method', sep=del) %>% rlang::sym()
    inf_val <- paste('GCT',idx,'Score', sep=del) %>% rlang::sym()
    gct_tib <- tibble::tibble(!!inf_key := key, !!inf_val := safeGCT(sset, verbose=verbose,vt=vt+1,tc=tc+1,tt=tt) )
    
    sex_tib <- NULL
    inf_key <- paste('Sex',idx,'Method', sep=del) %>% rlang::sym()
    inf_val <- paste('Sex',idx,'Call', sep=del) %>% rlang::sym()
    sex_tib <- tibble::tibble(!!inf_key := key, !!inf_val := safeSex(sset, verbose=verbose,vt=vt+1,tc=tc+1,tt=tt) )
    
    kar_tib <- NULL
    inf_key <- paste('Karyotype',idx,'Method', sep=del) %>% rlang::sym()
    inf_val <- paste('Karyotype',idx,'Call', sep=del) %>% rlang::sym()
    kar_tib <- tibble::tibble(!!inf_key := key, !!inf_val := safeSexKaryo(sset, verbose=verbose,vt=vt+1,tc=tc+1,tt=tt) )
    
    eth_tib <- NULL
    inf_key <- paste('Ethnicity',idx,'Method', sep=del) %>% rlang::sym()
    inf_val <- paste('Ethnicity',idx,'Call', sep=del) %>% rlang::sym()
    eth_tib <- tibble::tibble(!!inf_key := key, !!inf_val := safeEthnicity(sset, verbose=verbose,vt=vt+1,tc=tc+1,tt=tt) )
    
    tib <- dplyr::bind_cols(pre, gct_tib, sex_tib, kar_tib, eth_tib)
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done; elapsed={etime}.{RET}{RET}"))
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  
  tib
}

safeGCT_org = function(sset, verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'safeGCT_org'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  stime <- system.time({
    val <- NA
    try_str <- 'pass'
    
    try(val <- sesame::bisConversionControl(sset), silent = TRUE)
    if (is.na(val)) try_str <- 'fail'
  })
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Value={val}, try_str={try_str}. Done.{RET}{RET}"))
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  
  val
}

safeGCT = function(sset, verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'safeGCT'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  stime <- system.time({
    val <- NULL
    
    try_str <- ''
    val = tryCatch({
      try_str <- 'Pass'
      suppressWarnings(sesame::bisConversionControl(sset) )
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
  })
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Value={val}, try_str={try_str}. Done.{RET}{RET}"))
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  
  val
}

safeEthnicity = function(sset, verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'safeEthnicity'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  stime <- system.time({
    val <- NULL
    
    try_str <- ''
    val = tryCatch({
      try_str <- 'Pass'
      suppressWarnings(sesame::inferEthnicity(sset) )
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
  })
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Value={val}, try_str={try_str}. Done.{RET}{RET}"))
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  
  val
}

safeSexKaryo = function(sset, verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'safeSexKaryo'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  stime <- system.time({
    val <- NULL
    
    try_str <- ''
    val = tryCatch({
      try_str <- 'Pass'
      suppressWarnings(sesame::inferSexKaryotypes(sset) )
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
  })
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Value={val}, try_str={try_str}. Done.{RET}{RET}"))
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  
  val
}

safeSex = function(sset, verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'safeSex'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  stime <- system.time({
    val <- NULL
    
    try_str <- ''
    val = tryCatch({
      try_str <- 'Pass'
      suppressWarnings(sesame::inferSex(sset) )
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
  })
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Value={val}, try_str={try_str}. Done.{RET}{RET}"))
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  
  val
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                      Sesame Beta To Predictions Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
safeSkinAge = function(beta, verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'safeSkinAge'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  stime <- system.time({
    val <- NULL
    
    try_str <- ''
    val = tryCatch({
      try_str <- 'Pass'
      suppressWarnings(sesame::predictAgeSkinBlood(betas=beta) )
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
  })
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Value={val}, try_str={try_str}. Done.{RET}{RET}"))
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  
  val
}

safePhenoAge = function(beta, verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'safePhenoAge'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  stime <- system.time({
    val <- NULL
    
    try_str <- ''
    val = tryCatch({
      try_str <- 'Pass'
      suppressWarnings(sesame::predictAgePheno(betas=beta) )
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
  })
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Value={val}, try_str={try_str}. Done.{RET}{RET}"))
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  
  val
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                        Sesame SSET to Tibs Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
ssetToPvalTib = function(sset, name, percision=0, verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'ssetToPvalTib'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} name={name}.{RET}"))
  stime <- system.time({
    dat <- sset@pval %>% tibble::enframe(name='Probe_ID', value=name)
    if (percision!=0) dat <- dat %>% dplyr::mutate_if(purrr::is_double, round, percision)
  })
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}{RET}"))
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  
  dat
}

ssetToBetaTib = function(sset, name, quality.mask=FALSE, nondetection.mask=FALSE, 
                         mask.use.tcga=FALSE, pval.threshold=1, sum.TypeI=FALSE,
                         as.enframe=FALSE, percision=0,
                         verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'ssetToBetaTib'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} name={name}.{RET}"))
  stime <- system.time({
    dat <- sesame::getBetas(sset, quality.mask=quality.mask, nondetection.mask=nondetection.mask, 
                            mask.use.tcga=mask.use.tcga, pval.threshold=pval.threshold, sum.TypeI=sum.TypeI)
    if (percision!=0) dat <- round(dat, percision)
    if (as.enframe) dat <- dat %>% tibble::enframe(name='Probe_ID', value=name)
  })
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}{RET}"))
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  
  dat
}

ssetToSigsTib = function(sset, add, name, del='_', rmAdd=TRUE, verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'ssetToSigsTib'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting, name={name}!{RET}"))
  
  stime <- system.time({
    dat <- NULL
    
    datTag <- 'sig'
    grnTag <- paste(name,'G', sep=del) %>% rlang::sym()
    redTag <- paste(name,'R', sep=del) %>% rlang::sym()
    inbTag <- paste(name,'Inb_Col', sep=del) %>% rlang::sym()
    swpTag <- paste(name,'Swap', sep=del) %>% rlang::sym()
    # grnTag <- 'G' %>% rlang::sym()
    # redTag <- 'R' %>% rlang::sym()
    
    IG_tib <- sset@IG %>% tibble::as_tibble(rownames='Probe_ID',.name_repair = "unique") %>%
      tidyr::gather(Design_Type, !!grnTag, -Probe_ID) %>%
      dplyr::mutate(Inb_Col='G')
    # dplyr::mutate(Design_Type=paste0(Design_Type,'G'))
    
    OR_tib <- sset@oobR %>% tibble::as_tibble(rownames='Probe_ID',.name_repair = "unique") %>%
      tidyr::gather(Design_Type, !!redTag, -Probe_ID) %>%
      dplyr::mutate(Inb_Col='G')
    # dplyr::mutate(Design_Type=paste0(Design_Type,'G'))
    
    IR_tib <- sset@IR %>% tibble::as_tibble(rownames='Probe_ID',.name_repair = "unique") %>%
      tidyr::gather(Design_Type, !!redTag, -Probe_ID) %>%
      dplyr::mutate(Inb_Col='R')
    # dplyr::mutate(Design_Type=paste0(Design_Type,'R'))
    
    OG_tib <- sset@oobG %>% tibble::as_tibble(rownames='Probe_ID',.name_repair = "unique") %>%
      tidyr::gather(Design_Type, !!grnTag, -Probe_ID) %>%
      dplyr::mutate(Inb_Col='R')
    # dplyr::mutate(Design_Type=paste0(Design_Type,'R'))
    
    I1 <- dplyr::bind_rows(dplyr::inner_join(IG_tib,OR_tib, by=c('Probe_ID', 'Inb_Col', 'Design_Type')),
                           dplyr::inner_join(OG_tib,IR_tib, by=c('Probe_ID', 'Inb_Col', 'Design_Type')) ) %>%
      dplyr::mutate(Design_Type=paste0(Design_Type,'I')) %>%
      dplyr::select(Probe_ID,Design_Type,Inb_Col,everything())
    
    I1 <- add %>% dplyr::filter(Design_Type!='II') %>% 
      dplyr::left_join(I1, by=c("Probe_ID", "Design_Type") ) %>%
      dplyr::mutate(Swap=case_when(Man_Col!=Inb_Col ~ TRUE,
                                   Man_Col==Inb_Col ~ FALSE,
                                   TRUE ~ NA)) %>%
      dplyr::select(-Inb_Col) %>%
      dplyr::select(Probe_ID, Man_Col, Design_Type, Probe_Type, Swap, everything())
    # dplyr::select(Probe_ID, Address, Man_Col, Design_Type, Probe_Type, Inb_Col, Swap, everything())
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #  Infinium II::
    I2 <- sset@II %>% tibble::as_tibble(rownames='Probe_ID',.name_repair = "unique") %>%
      dplyr::rename(!!redTag :=U, !!grnTag :=M) %>%
      dplyr::mutate(Inb_Col=NA, Swap=NA, Design_Type='II')
    
    I2 <- add %>% dplyr::filter(Design_Type=='II') %>% 
      dplyr::left_join(I2, by=c("Probe_ID", "Design_Type") ) %>%
      dplyr::select(-Inb_Col) %>%
      dplyr::select(Probe_ID, Man_Col, Design_Type, Probe_Type, Swap, everything())
    # dplyr::select(Probe_ID, Address, Man_Col, Design_Type, Probe_Type, Inb_Col, Swap, everything())
    
    dat <- dplyr::bind_rows(I1, I2) %>%
      dplyr::rename(!!swpTag := Swap) %>%
      dplyr::arrange(Probe_ID)
    
    if (rmAdd && "Address" %in% names(dat)) dat <- dat %>% dplyr::select(-Address)
  })
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}{RET}"))
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  
  dat
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                        Sesame Tibs Summary Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
sigsToSummary = function(tib, name, del='_', verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'sigsToSummary'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  stime <- system.time({
    dat <- NULL
    
    cntTag <- 'Count' %>% rlang::sym()
    swapTag <- 'Swap' %>% rlang::sym()
    scntTag <- 'Swap_Count' %>% rlang::sym()
    
    if (! is.null(name)) {
      cntTag  <- paste(name,'Count', sep=del) %>% rlang::sym()
      swapTag <- paste(name,'Swap', sep=del) %>% rlang::sym()
      scntTag <- paste(name,'Swap_Count', sep=del) %>% rlang::sym()
    }
    
    # Remove Address from summary if address exists::
    if ("Address" %in% names(tib)) tib <- tib %>% dplyr::select(-Address) 
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #   Summarize Data::
    cnts <- tib %>% dplyr::group_by(Design_Type, Probe_Type, Man_Col) %>%
      summarise(!!cntTag :=n(), !!scntTag := count(!!swapTag==TRUE))
    sums <- tib %>% dplyr::group_by(Design_Type, Probe_Type, Man_Col) %>%
      summarise_if(is.numeric, list(min=min, avg=mean, med=median, max=max), na.rm=TRUE )
    
    dat <- cnts %>% dplyr::full_join(sums, by=c("Design_Type","Probe_Type", "Man_Col")) %>%
      dplyr::arrange(Probe_Type)
  })
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}{RET}"))
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  
  dat
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                          Sesame Manifest Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
getSesameManifest = function(man, sig, verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'getSesameManifest'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  stime <- system.time({
    tib <- man %>% dplyr::right_join(sig, by=c("U"="Address")) %>% 
      dplyr::filter(!is.na(Probe_ID)) %>%
      dplyr::select(Probe_ID, M, U, DESIGN, COLOR_CHANNEL, col, Probe_Type, Probe_Source, Next_Base) %>%
      dplyr::arrange(Probe_ID)
  })
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}{RET}"))
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  
  tib
}



# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#             Every thing below is not used and can be removed::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #



# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                   OLD Sesame SSET Manipulation Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
sesameStepAbbreviation = function(x, verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'sesameStepAbbreviation'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  if (x=='raw') return('R')
  if (x=='dyeBiasCorrTypeINorm') return('D')
  if (x=='detectionPnegEcdf') return('N')
  if (x=='pOOBAH') return('P')
  if (x=='noob') return('N')
  if (x=='noobsb') return('S')
  if (x=='inferTypeIChannel') return('I')
  stop(glue::glue("{RET}[{funcTag}]: ERROR: Unsupported Sesame Abbreviation={x}!!!{RET}{RET}") )
  
  return('U')
}

sesameWorkflow = function(sset=NULL, add, call, sigs, swap, pheno, beadPool=NULL,
                          prefix=NULL, platform=NULL, manifest=NULL, # This is only used if sset is not present
                          stepCalls=NULL,
                          negsCalls=NULL,
                          poobCalls=NULL,
                          betaCalls=NULL, 
                          intsCalls=NULL,
                          swapCalls=NULL, 
                          fenoCalls=NULL, 
                          del='_', verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'sesameWorkflow'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  stime <- system.time({
    
    if (is.null(sset)) {
      stopifnot(is.null(prefix))
      stopifnot(is.null(platform))
      stopifnot(is.null(manifest))
      sset <- initSesameRaw(prefix=prefix, platform=platform, manifest=manifest,
                            verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    }
    
    sabr <- ''
    step_cnt <- length(stepCalls)
    for (ii in seq(1:step_cnt)) {
      step <- stepCalls[ii]
      sabr <- paste0(sabr, sesameStepAbbreviation(step) )
      if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} step={step}, sabr={sabr}.{RET}"))
      sset <- sset %>% mutateSesame(method=step, verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
      if (verbose>=vt+3) print(sset)
      
      # Make Negative Detection-Pval Call::
      if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Starting pval.{RET}"))
      
      if (!is.null(negsCalls) && !is.null(negsCalls[ii]) && negsCalls[ii]==TRUE) {
        nabr <- paste(sabr,'negs_pval', sep=del)
        call <- call %>% dplyr::left_join(
          sset %>% mutateSesame(method='detectionPnegEcdf', verbose=verbose,vt=vt+1,tc=tc+1,tt=tt) %>%
            ssetToPvalTib(name=nabr, verbose=verbose,vt=vt+1,tc=tc+1,tt=tt), by="Probe_ID")
      }
      # Make Poob Detection-Pval Call::
      if (!is.null(poobCalls) && !is.null(poobCalls[ii]) && poobCalls[ii]==TRUE) {
        pabr <- paste(sabr,'poob_pval', sep=del)
        call <- call %>% dplyr::left_join(
          sset %>% mutateSesame(method='pOOBAH', verbose=verbose,vt=vt+1,tc=tc+1,tt=tt) %>%
            ssetToPvalTib(name=pabr, verbose=verbose,vt=vt+1,tc=tc+1,tt=tt), by="Probe_ID")
      }
      
      # Make Beta Call::
      if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Starting beta.{RET}"))
      beta <- NULL
      if (!is.null(betaCalls) && !is.null(betaCalls[ii]) && betaCalls[ii]==TRUE) {
        babr <- paste(sabr,'beta', sep=del)
        beta <- ssetToBetaTib(sset=sset, name=babr, as.enframe=FALSE, verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
        call <- call %>% dplyr::left_join(tibble::enframe(beta, name='Probe_ID', value=babr), by="Probe_ID")
      }
      
      # Make Signal Intesnisty Table::
      if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Starting signal{RET}"))
      cur_sigs <- NULL
      if (!is.null(intsCalls) && !is.null(intsCalls[ii]) && intsCalls[ii]==TRUE) {
        # iabr <- paste(sabr,'sigs', sep=del)
        iabr <- sabr
        cur_sigs <- ssetToSigsTib(sset=sset, add=add, name=iabr, verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
        
        sigs <- sigs %>% dplyr::left_join(cur_sigs, by=c("Probe_ID", "Man_Col", "Design_Type", "Probe_Type") )
      }
      
      # Record Swapping Changes::
      if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Starting swapping.{RET}"))
      if (!is.null(intsCalls) && !is.null(intsCalls[ii]) && intsCalls[ii]==TRUE) {
        # wabr <- paste(sabr,'swap', sep=del)
        wabr <- sabr
        
        # Need to make the signal call before if it wasn't already called::
        if (is.null(cur_sigs))
          cur_sigs <- ssetToSigsTib(sset=sset, add=add, name=iabr, verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
        
        if (verbose>=vt+2) cat(glue::glue("{RET}{RET}CURRENT SIGS({wabr}){RET}"))
        if (verbose>=vt+2) print(cur_sigs)
        
        # TBD:: Fix the directly indexed field[5] below...
        cur_swap <- cur_sigs %>% dplyr::filter(Design_Type!='II') %>% 
          dplyr::filter(.[5]==TRUE) %>% dplyr::select(Probe_ID, Man_Col, Probe_Type, 5)
        
        if (verbose>=vt+2) cat(glue::glue("{RET}{RET}CURRENT SWAPS({wabr}){RET}"))
        if (verbose>=vt+2) print(cur_swap)
        if (is.null(swap)) {
          swap <- cur_swap
        } else {
          swap <- swap %>% dplyr::full_join(cur_swap, by=c("Probe_ID", "Man_Col", "Probe_Type")) %>% dplyr::distinct()
        }
      }
      
      # Make Sesame Inferred Phenotype Calls::
      if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Starting Inferred Phenotype calls.{RET}"))
      if (!is.null(fenoCalls) && !is.null(fenoCalls[ii]) && fenoCalls[ii]==TRUE) {
        if (verbose>=vt+1) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Calculating Pheno Type Predictions(beadPool={beadPool}).{RET}"))
        
        # fabr <- paste(sabr,'pheno', sep=del)
        fabr <- sabr
        if (is.null(beta)) beta <- ssetToBetaTib(sset=sset, name=fabr, as.enframe=FALSE, verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
        
        gct_str  <- paste('gct',fabr, sep=del)
        age_str  <- paste('agePheno',fabr, sep=del)
        skin_str <- paste('ageSkin',fabr, sep=del)
        sex_str  <- paste('sex',fabr, sep=del)
        sexKaryo_str  <- paste('sexKaryo',fabr, sep=del)
        ethnicity_str <- paste('ethnicity',fabr, sep=del)
        
        if (!is.null(beadPool) && beadPool=='EPIC') {
          gct  <- safeGCT(sset=sset, verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
          age  <- safePhenoAge(beta, verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
          skin <- safeSkinAge(beta, verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
          sex  <- safeSex(sset, verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
          sexKaryo  <- safeSexKaryo(sset, verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
          ethnicity <- safeEthnicity(sset, verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
        } else {
          gct  <- NA
          age  <- NA
          skin <- NA
          sex  <- NA
          sexKaryo  <- NA
          ethnicity <- NA
        }
        
        phen_tib <- tibble::tibble(!!gct_str  := gct,
                                   !!age_str  := age,
                                   !!skin_str := skin,
                                   !!sex_str  := sex,
                                   !!sexKaryo_str  := sexKaryo,
                                   !!ethnicity_str := ethnicity
        )
        
        pheno <- pheno %>% dplyr::bind_cols(phen_tib)
      }
    }
    
  })
  nrows <- call %>% base::nrow()
  etime <- stime[3] %>% as.double() %>% round(2)
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done. nrows={nrows}, elapsed={etime}.{RET}{RET}"))
  
  list(call, sigs, swap, pheno, sset)
}


# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                     SSET to BeadSET Conversion Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
ssetBeta2bset = function(sset, bset, nkey, del='_', 
                         quality.mask=FALSE, nondetection.mask=FALSE,
                         mask.use.tcga=FALSE, pval.threshold=1, sum.TypeI=FALSE,
                         verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'ssetBeta2bset'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  datTag  <- 'beta'
  betaTag <- paste(nkey,datTag, sep=del) # %>% rlang::sym()
  
  stime <- system.time({
    dat <- NULL
    
    betas <- sset %>% sesame::getBetas(quality.mask=quality.mask, nondetection.mask=nondetection.mask, 
                                       mask.use.tcga=mask.use.tcga, pval.threshold=pval.threshold, sum.TypeI=sum.TypeI) %>%
      tibble::enframe(name='Probe_ID', value=betaTag)
    
    dat <- add2bset(bset=bset, inf1=betas, keyA='Probe_ID', verbose=verbose,vt=1,tc=0,tt=tt)
  })
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}{RET}"))
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  
  dat
}

ssetPval2bset = function(sset, bset, nkey, pkey, del='_', verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'ssetPval2bset'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  datTag  <- 'pval'
  pvalTag <- paste(nkey,pkey,datTag, sep=del) # %>% rlang::sym()
  
  stime <- system.time({
    dat <- NULL
    
    pvals <- sset@pval %>% tibble::enframe(name='Probe_ID', value=pvalTag)
    dat <- add2bset(bset=bset, inf1=pvals, keyA='Probe_ID', verbose=verbose,vt=1,tc=0,tt=tt)
  })
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}{RET}"))
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  
  dat
}

ssetSigs2bset = function(sset, bset, nkey, del='_', verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'ssetSig2bset'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  stime <- system.time({
    dat <- NULL
    
    datTag <- 'sig'
    grnTag <- paste(nkey,'Grn',datTag, sep=del) %>% rlang::sym()
    redTag <- paste(nkey,'Red',datTag, sep=del) %>% rlang::sym()
    
    IG_tib <- sset@IG %>% tibble::as_tibble(rownames='Probe_ID',.name_repair = "unique") %>% 
      tidyr::gather(Design_Type, !!grnTag, -Probe_ID) # %>%
    # dplyr::mutate(Design_Type=as.factor(Design_Type))
    
    IR_tib <- sset@IR %>% tibble::as_tibble(rownames='Probe_ID',.name_repair = "unique") %>% 
      tidyr::gather(Design_Type, !!redTag, -Probe_ID) # %>%
    # dplyr::mutate(Design_Type=as.factor(Design_Type))
    
    OG_tib <- sset@oobG %>% tibble::as_tibble(rownames='Probe_ID',.name_repair = "unique") %>% 
      tidyr::gather(Design_Type, !!grnTag, -Probe_ID) # %>%
    # dplyr::mutate(Design_Type=as.factor(Design_Type))
    
    OR_tib <- sset@oobR %>% tibble::as_tibble(rownames='Probe_ID',.name_repair = "unique") %>% 
      tidyr::gather(Design_Type, !!redTag, -Probe_ID) # %>%
    # dplyr::mutate(Design_Type=as.factor(Design_Type))
    
    I1 <- dplyr::bind_rows(dplyr::inner_join(IG_tib,OR_tib, by=c('Probe_ID', 'Design_Type')),
                           dplyr::inner_join(OG_tib,IR_tib, by=c('Probe_ID', 'Design_Type')) )
    
    I2 <- sset@II %>% tibble::as_tibble(rownames='Probe_ID',.name_repair = "unique") %>% 
      dplyr::rename(!!redTag :=U, !!grnTag :=M)
    
    dat <- add2bset(bset=bset, inf1=I1, inf2=I2, keyA='Probe_ID',keyB='Design_Type', 
                    verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
    
  })
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}{RET}"))
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  
  dat
}

add2bset = function(bset, inf1, inf2=NULL, keyA=NULL,keyB=NULL,
                    verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'add2bset'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  stime <- system.time({
    dat <- NULL
    if (is.null(inf2)) inf2 <- inf1
    if (is.null(keyB)) {
      dat[['I1']] <- bset[[1]] %>% dplyr::inner_join(inf1, by=c(keyA))
    } else {
      dat[['I1']] <- bset[[1]] %>% dplyr::inner_join(inf1, by=c(keyA,keyB))
    }
    dat[['I2']] <- bset[[2]] %>% dplyr::inner_join(inf2, by=c(keyA))
  })
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}{RET}"))
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  
  dat
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                          Function Template::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
templateFuncT = function(x, verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'templateFuncT'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  # tabsStr <- paste0(rep(' ', length(funcTag)),tabStr, rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  stime <- system.time({
    ## PUT CODE HERE
  })
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}{RET}"))
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  
  x
}



# End of file
