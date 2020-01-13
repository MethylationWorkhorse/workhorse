
suppressWarnings(suppressPackageStartupMessages(require("tidyverse")) )
suppressWarnings(suppressPackageStartupMessages(require("stringr")) )
suppressWarnings(suppressPackageStartupMessages(require("glue")) )

COM <- ","
TAB <- "\t"
RET <- "\n"


# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                    Sesame SSET Manipulation Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
sesameStepAbbreviation = function(x, verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'sesameStepAbbreviation'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  if (x=='Raw') return('R')
  if (x=='dyeBiasCorrTypeINorm') return('D')
  if (x=='detectionPnegEcdf') return('N')
  if (x=='pOOBAH') return('P')
  if (x=='noob') return('N')
  if (x=='noobsb') return('S')
  if (x=='inferTypeIChannel') return('I')
  stop(glue::glue("{RET}[{funcTag}]: ERROR: Unsupported Sesame Abbreviation={x}!!!{RET}{RET}") )
  
  return('U')
}

sesameWorkflow = function(sset=NULL, add, call, sigs, pheno, beadPool=NULL,
                          prefix=NULL, platform=NULL, manifest=NULL, # This is only used if sset is not present
                          stepCalls=NULL,
                          negsCalls=NULL,
                          poobCalls=NULL,
                          betaCalls=NULL, 
                          intsCalls=NULL,
                          fenoCalls=NULL, del='_',
                          verbose=0,vt=3,tc=1,tt=NULL) {
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
      
      sset <- sset %>% mutateSesame(method=step, verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
      
      # Make Negative Detection-Pval Call::
      if (!is.null(negsCalls) && !is.null(negsCalls[ii]) && negsCalls[ii]==TRUE) {
        nabr <- paste(sabr,'negs_pval', sep=del)
        call <- call %>% dplyr::left_join(
          sset %>% mutateSesame(method='detectionPnegEcdf', verbose=opt$verbosity,vt=vt+1,tc=tc+1,tt=tt) %>%
            ssetToPvalTib(name=nabr, verbose=opt$verbosity,vt=vt+1,tc=tc+1,tt=tt), by="Probe_ID")
      }
      # Make Poob Detection-Pval Call::
      if (!is.null(poobCalls) && !is.null(poobCalls[ii]) && poobCalls[ii]==TRUE) {
        pabr <- paste(sabr,'poob_pval', sep=del)
        call <- call %>% dplyr::left_join(
          sset %>% mutateSesame(method='pOOBAH', verbose=opt$verbosity,vt=vt+1,tc=tc+1,tt=tt) %>%
            ssetToPvalTib(name=pabr, verbose=opt$verbosity,vt=vt+1,tc=tc+1,tt=tt), by="Probe_ID")
      }
      
      # Make Beta Call::
      beta <- NULL
      if (!is.null(betaCalls) && !is.null(betaCalls[ii]) && betaCalls[ii]==TRUE) {
        babr <- paste(sabr,'beta', sep=del)
        beta <- ssetToBetaTib(sset=sset, name=babr, as.enframe=FALSE, verbose=opt$verbosity,vt=vt+1,tc=tc+1,tt=tt)
        call <- call %>% dplyr::left_join(tibble::enframe(beta, name='Probe_ID', value=babr), by="Probe_ID")
      }
      
      # Make Intesnisty Call::
      if (!is.null(intsCalls) && !is.null(intsCalls[ii]) && intsCalls[ii]==TRUE) {
        # iabr <- paste(sabr,'sigs', sep=del)
        iabr <- sabr
        sigs <- sigs %>% dplyr::left_join(
          ssetToSigsTib(sset=sset, add=add, name=iabr, verbose=opt$verbosity,vt=vt+1,tc=tc+1,tt=tt), 
          by=c("Probe_ID", "Man_Col", "Design_Type", "Probe_Type") )
      }
      
      # Make Sesame Inference Calls::
      if (!is.null(fenoCalls) && !is.null(fenoCalls[ii]) && fenoCalls[ii]==TRUE) {
        if (verbose>=vt+1) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Calculating Pheno Type Predictions(beadPool={beadPool}).{RET}"))
        
        # fabr <- paste(sabr,'pheno', sep=del)
        fabr <- sabr
        if (is.null(beta)) beta <- ssetToBetaTib(sset=sset, name=fabr, as.enframe=FALSE, verbose=opt$verbosity,vt=vt+1,tc=tc+1,tt=tt)
        
        gct_str  <- paste('gct',fabr, sep=del)
        age_str  <- paste('agePheno',fabr, sep=del)
        skin_str <- paste('ageSkin',fabr, sep=del)
        sex_str  <- paste('sex',fabr, sep=del)
        sexKaryo_str  <- paste('sexKaryo',fabr, sep=del)
        ethnicity_str <- paste('ethnicity',fabr, sep=del)

        if (!is.null(beadPool) && beadPool=='EPIC') {
          gct  <- safeGCT(sset=sset, verbose=opt$verbosity,vt=vt+1,tc=tc+1,tt=tt)
          age  <- safePhenoAge(beta, verbose=opt$verbosity,vt=vt+1,tc=tc+1,tt=tt)
          skin <- safeSkinAge(beta, verbose=opt$verbosity,vt=vt+1,tc=tc+1,tt=tt)
          sex  <- safeSex(sset, verbose=opt$verbosity,vt=vt+1,tc=tc+1,tt=tt)
          sexKaryo  <- safeSexKaryo(sset, verbose=opt$verbosity,vt=vt+1,tc=tc+1,tt=tt)
          ethnicity <- safeEthnicity(sset, verbose=opt$verbosity,vt=vt+1,tc=tc+1,tt=tt)
        } else {
          cat("\n\n\nBEAD_POOL(NULL)=",beadPool,"\n\n\n",sep='')
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
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}{RET}"))
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  
  # list(sset, call, sigs)
  # list(call, sigs, pheno)
  list(call, sigs, pheno)
}

initSesameRaw = function(prefix, platform, manifest, verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'initSesameRaw'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  stime <- system.time({
    sset <- sesame::readIDATpair(prefix, platform=platform, manifest=manifest)
  })
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}{RET}"))
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  
  sset
}

mutateSesame = function(sset, method, verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'mutateSesame'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Mutating Sesame({method}){RET}"))
  stime <- system.time({
    if (is.null(method)) stop(glue::glue("{RET}[{funcTag}]: ERROR: Missing method!!!{RET}{RET}"))
    else if (method=='dyeBiasCorrTypeINorm') sset <- sset %>% sesame::dyeBiasCorrTypeINorm()
    else if (method=='detectionPnegEcdf') sset <- sset %>% sesame::detectionPnegEcdf()
    else if (method=='pOOBAH') sset <- sset %>% sesame::pOOBAH()
    else if (method=='noob') sset <- sset %>% sesame::noob()
    else if (method=='noobsb') sset <- sset %>% sesame::noobsb()
    else if (method=='inferTypeIChannel') sset <- sset %>% sesame::inferTypeIChannel(switch_failed=TRUE, verbose=FALSE)
    else if (method=='Raw') { } # sset <- sset
    else stop(glue::glue("{RET}[{funcTag}]: ERROR: Unsupported method={method}!!!{RET}{RET}"))
  })
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}{RET}"))
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  
  sset
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                         Sesame SSET Call Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
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



# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                        Sesame SSET to Tibs Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
ssetToPvalTib = function(sset, name, verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'ssetToPvalTib'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} name={name}.{RET}"))
  stime <- system.time({
    dat <- sset@pval %>% tibble::enframe(name='Probe_ID', value=name)
  })
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}{RET}"))
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  
  dat
}

ssetToBetaTib = function(sset, name, quality.mask=FALSE, nondetection.mask=FALSE, 
                         mask.use.tcga=FALSE, pval.threshold=1, sum.TypeI=FALSE,
                         as.enframe=FALSE,
                         verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'ssetToBetaTib'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} name={name}.{RET}"))
  stime <- system.time({
    dat <- sesame::getBetas(sset, quality.mask=quality.mask, nondetection.mask=nondetection.mask, 
                            mask.use.tcga=mask.use.tcga, pval.threshold=pval.threshold, sum.TypeI=sum.TypeI)
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
    
    # CODE:: To add Color String not sure if its needed for calculating Swap Count...
    # ses_man_tib %>% dplyr::mutate(PDesign_Key=case_when(DESIGN=='II' ~ DESIGN,
    #                                                     TRUE ~ paste0(col,DESIGN)))
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
    
    dat <- add2bset(bset=bset, inf1=betas, keyA='Probe_ID', verbose=opt$verbosity,vt=1,tc=0,tt=tTracker)
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
    dat <- add2bset(bset=bset, inf1=pvals, keyA='Probe_ID', verbose=opt$verbosity,vt=1,tc=0,tt=tTracker)
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
templateT = function(x, verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'templateT'
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
