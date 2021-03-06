
suppressWarnings(suppressPackageStartupMessages( base::require("tidyverse") ))
suppressWarnings(suppressPackageStartupMessages( base::require("stringr") ))
suppressWarnings(suppressPackageStartupMessages( base::require("glue") ))

suppressWarnings(suppressPackageStartupMessages( base::require("matrixStats") ))
suppressWarnings(suppressPackageStartupMessages( base::require("scales") ))
suppressWarnings(suppressPackageStartupMessages( base::require("grid") ))

# Plotting Extras
suppressWarnings(suppressPackageStartupMessages( base::require("GGally") ))
suppressWarnings(suppressPackageStartupMessages( base::require("hexbin") ))

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                   Per CpG Screening/Summary Functions::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
screenTitration = function(tib, summary=FALSE,
                           verbose=0,vt=5,tc=1,tt=NULL) {
  funcTag <- 'flagLoci'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  # Single Sample Filtering:: beta variation
  flag_tib <- NULL
  flag_tib <- flagLoci(tib=tib, field='beta_31', cutoff=0.10, min=FALSE, cmp='gt', pre=flag_tib)
  flag_tib <- flagLoci(tib=tib, field='beta_sd', cutoff=0.10, min=FALSE, cmp='gt', pre=flag_tib)
  
  # Single Sample Filtering:: deltaBeta variation
  flag_tib <- flagLoci(tib=tib, field='dB_q2',  cutoff=0.05, min=FALSE, cmp='gt', pre=flag_tib)
  flag_tib <- flagLoci(tib=tib, field='dB_mu',  cutoff=0.05, min=FALSE, cmp='gt', pre=flag_tib)
  
  # Cross Sample Filtering:: CSS Scores
  flag_tib <- flagLoci(tib=tib, field='CSS_per',  cutoff=90, min=TRUE, cmp='lt', pre=flag_tib)
  flag_tib <- flagLoci(tib=tib, field='CSS_q2', cutoff=0.10, min=TRUE, cmp='lt', pre=flag_tib)
  flag_tib <- flagLoci(tib=tib, field='CSS_mu', cutoff=0.10, min=TRUE, cmp='lt', pre=flag_tib)
  # TBD:: Add cutoff for CSS_mu and CSS_q2 for each combinaiton directly (i.e. T50BZ_T99BZ_CSS_mu, cutoff=0.1, min=TRUE, cmp=lte)
  # flag_tib <- flagLoci(tib=tib, field='T00BZ_T50BZ_CSS_mu', cutoff=90, min=TRUE, cmp='lt', pre=flag_tib)
  # flag_tib <- flagLoci(tib=tib, field='T50BZ_T99BZ_CSS_mu', cutoff=90, min=TRUE, cmp='lt', pre=flag_tib)
  
  flag_tib <- dplyr::inner_join(dplyr::select(tib, Probe_ID, RSquared), flag_tib, by="Probe_ID")
  if (verbose>=vt+1) print(flag_tib)

  # Summary of all flags::
  if (summary) {
    flag_tib <- flag_tib %>% dplyr::select(RSquared, ends_with('_flag')) %>% group_by_if(is_logical) %>% 
      dplyr::summarise(Count=n(), R2_mu=mean(RSquared, na.rm=TRUE), R2_q2=median(RSquared, na.rm=TRUE)) %>% 
      dplyr::arrange(-Count) # %>% as.data.frame()
    # flag_tib <- flag_tib %>% dplyr::select(ends_with('_flag')) %>% group_by_all() %>% dplyr::summarise(Count=n()) %>% dplyr::arrange(-Count) # %>% as.data.frame()
  }
  
  flag_tib
}

flagLoci = function(tib, field, cutoff, min=FALSE, cmp='gt', ids='Probe_ID', del='_', pre=NULL,
                    verbose=0,vt=5,tc=1,tt=NULL) {
  funcTag <- 'flagLoci'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  end_col <- paste0(del,field)
  ids_str <- ids
  ids_col <- ids_str %>% rlang::sym()
  rnames  <- tib %>% dplyr::pull(ids_col)
  val_col <- paste(field,'min', sep=del) %>% rlang::sym()
  flg_col <- paste(field,'min_flag', sep=del) %>% rlang::sym()
  if (!min) {
    val_col <- paste(field,'max', sep=del) %>% rlang::sym()
    flg_col <- paste(field,'max_flag', sep=del) %>% rlang::sym()
  }

  cur_tib <- tib %>% dplyr::select(ends_with(!!end_col)) %>% as.matrix()
  if (min)  cur_tib <- cur_tib %>% rowMins(na.rm=TRUE)
  if (!min) cur_tib <- cur_tib %>% rowMaxs(na.rm=TRUE)
  
  cur_tib <- cur_tib %>%
    as.data.frame(row.names = rnames) %>% 
    tibble::as_tibble(rownames=ids_str, .name_repair='unique') %>% 
    purrr::set_names(c(ids_col, val_col) ) %>% 
    dplyr::mutate(!!flg_col:=TRUE)
  if (cmp=='lte') {
    cur_tib <- cur_tib %>% dplyr::mutate(
      !!flg_col:=case_when(is.na(!!val_col) ~ FALSE, !!val_col==-Inf ~ FALSE, !!val_col==Inf ~ FALSE, !!val_col <= !!cutoff ~ FALSE, TRUE ~ TRUE))
  } else if (cmp=='lt') {
    cur_tib <- cur_tib %>% dplyr::mutate(
      !!flg_col:=case_when(is.na(!!val_col) ~ FALSE, !!val_col==-Inf ~ FALSE, !!val_col==Inf ~ FALSE, !!val_col < !!cutoff ~ FALSE, TRUE ~ TRUE))
  } else if (cmp=='gt') {
    cur_tib <- cur_tib %>% dplyr::mutate(
      !!flg_col:=case_when(is.na(!!val_col) ~ FALSE, !!val_col==-Inf ~ FALSE, !!val_col==Inf ~ FALSE, !!val_col > !!cutoff ~ FALSE, TRUE ~ TRUE))
  } else {
    stop(glue::glue("{RET}[{funcTag}]: ERROR: Unsupported comparison type={cmp}!!!{RET}{RET}"))
    return(NULL)
  }
  
  if (!is.null(pre)) cur_tib <- dplyr::left_join(pre, cur_tib, by=ids_str)
  
  cur_tib
}

plotEndsWith = function(tib, ending, del='_', xmin=NULL, xmax=NULL) {
  
  end_str <- paste0(del,ending)
  ending  <- ending %>% rlang::sym()
  
  stack <- tib %>% dplyr::select(ends_with(!!end_str)) %>% 
    purrr::set_names(stringr::str_remove(names(.), end_str) ) %>% 
    tidyr::gather(Sample, !!ending)
  # print(stack)

  gg <- ggplot2::ggplot(data=stack, aes(x=!!ending, color=Sample)) +
    ggplot2::geom_density(alpha=0.8, na.rm=TRUE)
  
  if (!is.null(xmin) && !is.null(xmax)) gg <- gg + xlim(c(xmin,xmax))

  gg
}

crossSampleLociRSquared = function(tib, id='Probe_ID', max=NULL, retData=FALSE,
                                   cpp.verbose=0,verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'crossSampleLociRSquared'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  if (!is.null(max)) tib <- tib %>% head(n=max)
  
  stime <- system.time({
    # Builds Samle Counts and Sample Names vectors::
    # sam_tib <- tib %>% dplyr::select(-id) %>% names() %>% stringr::str_remove('_.*$') %>% 
    sam_tib <- tib %>% dplyr::select(-all_of(id)) %>% names() %>% stringr::str_remove('_.*$') %>% 
      tibble::enframe() %>% dplyr::group_by(value) %>% summarise(Count=n())
    sam_vec <- sam_tib %>% dplyr::pull(1) %>% as.vector()
    cnt_vec <- sam_tib %>% dplyr::pull(2) %>% as.vector()
    
    # Builds Summary Matrix (c++) and converts to tibble::
    sam_mat <- tib %>% column_to_rownames(id) %>% as.matrix()
    sum_mat <- C_crossSampleLociRSquared(sam_mat, cnt_vec, sam_vec, verbose=cpp.verbose)
    
    # TBD:: Need to fix this::
    sum_tib <- sum_mat %>% tibble::as_tibble(rownames=id)
    # sum_tib <- sum_mat
    
    # cpg_vec <- tib %>% dplyr::pull(1) %>% as.vector()
  })
  etime <- stime[3] %>% as.double() %>% round(2)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done; elapsed={etime}.{RET}{RET}"))
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  
  sum_tib
}

getCrossSampleLociVariation = function(tib, id='Probe_ID', max=NULL, retData=FALSE,
                                       verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'getLociSampleSummary'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  if (!is.null(max)) tib <- tib %>% head(n=max)
  
  stime <- system.time({
    # Builds Samle Counts and Sample Names vectors::
    # sam_tib <- tib %>% dplyr::select(-id) %>% names() %>% stringr::str_remove('_.*$') %>% 
    sam_tib <- tib %>% dplyr::select(-all_of(id)) %>% names() %>% stringr::str_remove('_.*$') %>% 
      tibble::enframe() %>% dplyr::group_by(value) %>% summarise(Count=n())
    sam_vec <- sam_tib %>% dplyr::pull(1) %>% as.vector()
    cnt_vec <- sam_tib %>% dplyr::pull(2) %>% as.vector()
    
    # HACK for properly labeled titration samples i.e. TXX (where X is a digit)
    dif_vec <- sam_vec %>% stringr::str_remove('T') %>% as.numeric() %>% cbind() %>% colDiffs() /100 %>% as.vector()
    # print(dif_vec)

    # Builds Summary Matrix (c++) and converts to tibble::
    sam_mat <- tib %>% column_to_rownames(id) %>% as.matrix()
    sum_mat <- C_crossSampleLociVariation(sam_mat, cnt_vec, sam_vec, dif_vec, verbose=0)
    sum_tib <- sum_mat %>% tibble::as_tibble(rownames=id)
    
    # TBD:: Old Code that should be deleted::
    cpg_vec <- tib %>% dplyr::pull(1) %>% as.vector()
    # sum_mat <- C_lociRSquared(sam_mat, cnt_vec, sam_vec, cpg_vec)
    # id <- id %>% rlang::sym()
    # sum_tib <- sum_mat %>% tibble::as_tibble() %>% tibble::add_column(!!id := cpg_vec) %>% dplyr::select(!!id, everything())
  })
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}{RET}"))
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  
  if (retData) {
    rdat <- NULL
    rdat$cntVec <- cnt_vec
    rdat$samVec <- sam_vec
    rdat$cpgVec <- cpg_vec
    rdat$samMat <- sam_mat
    rdat$sumMat <- sum_mat
    rdat$sumTib <- sum_tib
    
    return(rdat)
  }
  
  sum_tib
}

getSingleSampleLociVariation = function(tib, id='Probe_ID', max=NULL, retData=FALSE,
                                        verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'getLociSampleSummary'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  if (!is.null(max)) tib <- tib %>% head(n=max)
  
  stime <- system.time({
    # Builds Samle Counts and Sample Names vectors::
    # sam_tib <- tib %>% dplyr::select(-id) %>% names() %>% stringr::str_remove('_.*$') %>% 
    sam_tib <- tib %>% dplyr::select(-all_of(id)) %>% names() %>% stringr::str_remove('_.*$') %>% 
      tibble::enframe() %>% dplyr::group_by(value) %>% summarise(Count=n())
    sam_vec <- sam_tib %>% dplyr::pull(1) %>% as.vector()
    cnt_vec <- sam_tib %>% dplyr::pull(2) %>% as.vector()
    
    # Builds Summary Matrix (c++) and converts to tibble::
    #  This was the old name of this function::
    #      sum_mat <- C_lociRSquared(sam_mat, cnt_vec, sam_vec)
    sam_mat <- tib %>% column_to_rownames(id) %>% as.matrix()
    sum_mat <- C_singleSampleLociVariation(sam_mat, cnt_vec, sam_vec)
    sum_tib <- sum_mat %>% tibble::as_tibble(rownames=id)
    
    # TBD:: Old Code that should be deleted::
    cpg_vec <- tib %>% dplyr::pull(1) %>% as.vector()
    # sum_mat <- C_lociRSquared(sam_mat, cnt_vec, sam_vec, cpg_vec)
    # id <- id %>% rlang::sym()
    # sum_tib <- sum_mat %>% tibble::as_tibble() %>% tibble::add_column(!!id := cpg_vec) %>% dplyr::select(!!id, everything())
  })
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}{RET}"))
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  
  if (retData) {
    rdat <- NULL
    rdat$cntVec <- cnt_vec
    rdat$samVec <- sam_vec
    rdat$cpgVec <- cpg_vec
    rdat$samMat <- sam_mat
    rdat$sumMat <- sum_mat
    rdat$sumTib <- sum_tib
    
    return(rdat)
  }
  
  sum_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#  Below is the current code from plotByExperiments Jan-7-2019
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                            Clean Functions::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

selectSamples = function(ss, max=3, spread='mid', sort,
                         verbose=0,vt=4,tc=1) {
  funcTag <- 'selectSamples'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} max={max}, sort={sort}{RET}"),sep='')
  
  nrows <- ss %>% base::nrow()
  # Arrange by variable
  sort <- rlang::sym(sort)
  if (!is.null(sort)) ss <- ss %>% dplyr::arrange(desc(!!sort))
  ss <- ss %>% dplyr::mutate(Sample_Rep_Num=row_number())
  
  if (nrows==1 || nrows<=max) return(ss)
  
  if (verbose>=vt)
    cat(glue::glue("[{funcTag}]:{tabsStr} nrows={nrows}, max={max}, spread={spread}..."),"\n", sep='')
  
  # Select Top,Mid,Last
  sel_idxs <- NULL
  if (spread=='top') {
    sel_idxs <- c(1:max)
  } else if (spread=='mid') {
    mid_idx <- as.integer(nrows/2)+1
    sel_idxs <- c(1,mid_idx,nrows)
  } else if (spread=='bot') {
    end_idx  <- nrows-max
    if (end_idx<1) end_idx<- 1
    sel_idxs <- c(end_idx:nrows)
  } else {
    stop(glue::glue("[{funcTag}]: ERROR: Unsupported Sample Quality spread={spread}, nrows={nrows}"),"\n\n", sep='')
  }
  ss_sel <- ss %>% dplyr::slice(sel_idxs)
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}{RET}"),sep='')
  
  ss_sel
}

loadSampleRDS_byRep = function(ss, exp=NULL, minPval, pvalKey, betaKey,
                               parallel=FALSE, verbose=0,vt=3,tc=1) {
  funcTag <- 'loadSampleRDS_byRep'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  ss_split <- ss %>% split(.$Sample_Rep_Num)
  if (verbose>=vt)
    cat(glue::glue("[{funcTag}]:{tabsStr} Starting(isParallel={parallel}), minPval={minPval}..."),"\n", sep='')
  
  cur_tibs <- list()
  for (prefix in names(ss_split)) {
    exp_name <- prefix
    if (!is.null(exp)) exp_name <- paste(exp,exp_name, sep='.')
    prb_tib <- loadProbeRDS(ss_split[[prefix]]$Calls_Path, minPval=minPval, 
                            pvalKey=pvalKey, betaKey=betaKey, prefix=exp_name, 
                            verbose=verbose, vt=vt, tc=tc+1)
    cur_tibs[[prefix]] <- prb_tib
  }
  
  join_tib <- cur_tibs[[1]]
  join_ids <- names(cur_tibs)
  join_len <- join_ids %>% length()
  for (sIdx in c(2:join_len)) {
    join_tib <- join_tib %>% dplyr::left_join(cur_tibs[[ join_ids[sIdx] ]], by="Probe_ID")
  }
  if (verbose>=vt)
    cat(glue::glue("[{funcTag}]:{tabsStr} SampleCount={join_len}, Done."),"\n\n", sep='')
  
  join_tib
}

loadProbeRDS = function(file, minPval=NULL, prefix=NULL, pvalKey, betaKey,
                        verbose=0,vt=3,tc=1) {
  funcTag <- 'loadProbeRDS'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} minPva={minPval}, file={file}"),"\n", sep='')
  if (is.null(minPval)) stop(glue::glue("[{funcTag}]: ERROR: MUST PROVIDE minPval!!!"),"\n\n", sep='')
  
  pvalKey <- pvalKey %>% rlang::sym()
  betaKey <- betaKey %>% rlang::sym()
  
  if (stringr::str_ends(file,'.rds')) {
    tib <- suppressMessages(suppressWarnings(readr::read_rds(file)) ) %>% dplyr::bind_rows()
  } else {
    tib <- suppressMessages(suppressWarnings(readr::read_csv(file)) ) %>% dplyr::bind_rows()
  }
  # print(tib)
  # tib <- tib %>% dplyr::select(1,2,3,!!pvalKey, !!betaKey)
  tib <- tib %>% dplyr::select(1,!!pvalKey, !!betaKey)
  # print(tib)
  
  if (!is.null(minPval)) {
    
    tot_cnt  <- tib %>% base::nrow()
    pre_cnt  <- tib %>% dplyr::filter(is.na(!!betaKey)) %>% base::nrow()
    
    tib <- maskTib(tib, field=betaKey, pval=pvalKey, minPval=minPval, verbose=verbose) %>%
      dplyr::rename(Beta=!!betaKey)
    
    #  dplyr::rename(Pval=!!pvalKey, Beta=!!betaKey)
    # print(tib)
    # mask_idx <- which(tib$NegsDetP>minPval)
    # tib$Beta[mask_idx] <- NA
    
    # pos_cnt  <- tib %>% dplyr::filter(is.na(!!betaKey)) %>% base::nrow()
    pos_cnt  <- tib %>% dplyr::filter(is.na(Beta)) %>% base::nrow()
    pos_per  <- round(100*pos_cnt/tot_cnt,3)
    if (verbose>=vt)
      cat(glue::glue("[{funcTag}]:{tabsStr} Masked ({pos_per}%) Beta Values {pre_cnt} -> {pos_cnt}.{RET}"))
  }
  cname1 <- names(tib)[1]
  if (!is.null(prefix)) tib <- tib %>% purrr::set_names(paste(prefix,names(tib), sep='.') ) %>% dplyr::rename(!!cname1 :=1)
  
  tib
}


# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                Plotting Functions for Pairwise Comparison::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

plotBetaMatrix_bySample = function(tib, sample, manifest, minPval, pvalKey, betaKey,
                                   outDir, field='Beta', join='Probe_ID', sort, joinFull=FALSE,
                                   max=3, spread='mid', outType='auto', dpi=72, format='png', maxCnt=30000,
                                   retTibs=FALSE,
                                   verbose=0,vt=4,tc=1) {
  funcTag <- 'plotBetaMatrix_bySample'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  if (verbose>=vt)
    cat(glue::glue("[{funcTag}]:{tabsStr} Starting sample={sample}: minPval={minPval}, outType={outType}"),"\n", sep='')
  
  sam_dir <- file.path(outDir, 'samples', sample, spread)
  if (!dir.exists(sam_dir)) dir.create(sam_dir, recursive=TRUE)
  ss_exp_tibs  <- tib[[sample]] %>% split(.$Experiment_Key)
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #  Load Sample Probe Data from RDS
  prbs_exp_tibs <- list()
  idat_exp_tibs <- list()
  for (exp_key in names(ss_exp_tibs)) {
    exp_nrows <- ss_exp_tibs[[exp_key]] %>% base::nrow()
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} experiment={exp_key}, exp_nrows={exp_nrows}.{RET}"))
    
    loadRDS  <- TRUE
    loadIDAT <- FALSE
    if (exp_nrows>1) {
      ss_exp_tibs[[exp_key]]  <- selectSamples(ss_exp_tibs[[exp_key]], max=max, spread=spread, sort=sort,
                                               verbose=verbose,vt=vt,tc=tc+1)
      if (loadRDS)
        prbs_exp_tibs[[exp_key]] <- loadSampleRDS_byRep(ss=ss_exp_tibs[[exp_key]], exp=exp_key, 
                                                        minPval=minPval, pvalKey=pvalKey, betaKey=betaKey,
                                                        verbose=verbose+9,vt=vt+1,tc=tc+1)
      
      if (loadIDAT)
        idat_exp_tibs[[exp_key]] <- loadSampleIDAT_byRep(ss=ss_exp_tibs[[exp_key]], exp=exp_key, verbose=verbose+9,vt=vt+1,tc=tc+1)
    }
  }
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Loaded ALL Experimental Probe Data!{RET}{RET}"))
  
  if (retTibs) {
    ret <- NULL
    ret$ss_exp_tibs   <- ss_exp_tibs
    ret$prbs_exp_tibs <- prbs_exp_tibs
    # ret$idat_exp_tibs <- idat_exp_tibs
    return(ret)
  }
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #  Plot Sample Probe Data
  exp_prb_cnts <- length(prbs_exp_tibs)
  if (verbose>=vt)
    cat(glue::glue("[{funcTag}]:{tabsStr} sample={sample}, exp_prb_cnts={exp_prb_cnts}"),"\n", sep='')
  
  gg <- NULL
  exp_names <- NULL
  if (exp_prb_cnts>0) {
    exp_names <- names(prbs_exp_tibs)
    for (idxA in seq(exp_names)) {
      for (idxB in seq(exp_names)) {
        if (idxA < idxB) {
          nameA <- exp_names[idxA] %>% head(n=1)
          nameB <- exp_names[idxB] %>% head(n=1)
          
          prbsA <- prbs_exp_tibs[[nameA]]
          prbsB <- prbs_exp_tibs[[nameB]]
          
          idatA <- idat_exp_tibs[[nameA]]
          idatB <- idat_exp_tibs[[nameB]]
          
          if (joinFull) {
            join_prbs <- manifest %>% dplyr::right_join(dplyr::full_join(prbsA,prbsB, by="Probe_ID"), by="Probe_ID") %>%
              rename(Design_Type=DESIGN)
          } else {
            join_prbs <- manifest %>% dplyr::right_join(dplyr::inner_join(prbsA,prbsB, by="Probe_ID"), by="Probe_ID") %>%
              rename(Design_Type=DESIGN)
          }
          
          # join_prbs <- joinTibsAndManifest(datA=prbsA,keyA=nameA, datB=prbsB,keyB=nameB, 
          #                                  manifest=manifest, field=field, join=join,
          #                                  verbose=verbose, vt=vt+1, tc=tc+1)
          
          # join_idat <- joinTibsAndManifest(datA=prbsA,keyA=nameA, datB=prbsB,keyB=nameB, 
          #                                  manifest=manifest, field=field, join=join,
          #                                  verbose=verbose, vt=vt+1, tc=tc+1)
          
          # TBD::Ensure removing these and replacing with parameter variables is OK???
          # field_str <- 'NDI_beta' # betaKey
          # detP_str <- 'NDI_negs_pval' # pvalKey
          gg <- plotPairs(join_prbs, sample=sample, nameA=nameA, nameB=nameB, outDir=sam_dir,
                          probeType='cg', field=field, field_str=betaKey, detp=pvalKey, maxCnt=maxCnt, minPval=minPval,
                          # probeType='cg', field=field, field_str=field_str, detp=detP_str, maxCnt=maxCnt, minPval=minPval,
                          spread=spread, outType=outType, dpi=dpi, format=format,
                          verbose=verbose, tc=tc+1)
          
          # ret <- NULL
          # ret$gg <- gg
          # ret$sample <- sample
          # ret$nameA <- nameA
          # ret$nameB <- nameB
          # ret$sam_dir <- sam_dir
          # ret$prbsA <- prbsA
          # ret$prbsB <- prbsB
          # ret$idatA <- idatA
          # ret$idatB <- idatB
          # ret$join_prbs <- join_prbs
          # return(ret)
          
          # gg <- plotPvalViolins()
          
          # return(gg)
          # break
        }
      }
      # break
    }
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done."),"\n\n", sep='')
  }
  
  # exp_names
  # gg
  join_prbs
}

plotPairs = function(tib, sample, nameA, nameB, 
                     probeType='cg', groupA='Probe_Type', groupB='Design_Type',
                     field='Beta', field_str, detp,
                     outDir, 
                     maxCnt=30000, wsize=2.8, minPval, minDelta=0.2,
                     spread='mid', outType='auto', dpi=72, format='png',
                     datIdx=4,
                     verbose=0,vt=4,tc=1,tt=NULL) {
  funcTag <- 'plotPairs'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting sample={sample}: {nameA} vs. {nameB}{RET}"), sep='')
  
  plotDir  <- file.path(outDir, paste(nameA,'VS',nameB, sep='_'))
  plotName <- paste(sample,nameA,'VS',nameB,probeType,field_str, sep='_')
  if (!dir.exists(plotDir)) dir.create(plotDir, recursive=TRUE)
  
  # Filter out Controls::
  tib <- tib %>% dplyr::filter(Probe_Type=='cg'|Probe_Type=='ch'|Probe_Type=='rs')
  
  # Combine Probe_Type+Design_Type
  group <- 'Group'
  # group  <- rlang::sym(group)
  groupA <- rlang::sym(groupA)
  groupB <- rlang::sym(groupB)
  tib <- tib %>% dplyr::mutate(!!group := paste0(!!groupA,'.',!!groupB) )
  
  pdat <- tib %>% dplyr::select(!!group, dplyr::ends_with(field)) %>%
    purrr::set_names(stringr::str_remove(names(.), paste0('.',field,'$')) )
  ddat <- tib %>% dplyr::select(!!group, dplyr::ends_with(detp)) %>%
    purrr::set_names(stringr::str_remove(names(.), paste0('.',detp,'$')) )
  
  ncols_org <- base::ncol(pdat)
  nrows_org <- base::nrow(pdat)
  tarColIdxes <- c(2:ncols_org)
  
  if (nrows_org>maxCnt) sdat <- pdat %>% sampleGroup_n(n=maxCnt, field='Group')
  nrows_sub  <- base::nrow(sdat)
  nrows_per  <- round(100*nrows_sub/nrows_org,1)
  nrows_orgK <- number_as_commaK(nrows_org)
  nrows_subK <- number_as_commaK(nrows_sub)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Total={nrows_orgK}, Subset={nrows_subK}, Percent={nrows_per}{RET}"))
  
  gg_type  <- paste0('RSquared',field)
  name_pdf <- glue::glue("{plotName}.{format}")
  gg_pdf   <- file.path(plotDir, name_pdf)
  
  ptime <- Sys.time()
  gg_mtitle <- glue::glue("{nameA} VS {nameB}; sample={sample}")
  gg_stitle <- glue::glue("PlottedProbes={probeType}, Metric={field_str}, DetP={detp}, MinPval={minPval}")
  gg_ctitle <- glue::glue("TimeStamp={ptime}, DPI={dpi}, Plot displays downsampled percent = {nrows_per}% ",
                          "({nrows_subK}/{nrows_orgK})")
  
  alpha_hih <- 0.3
  alpha_mid <- 0.3
  alpha_low <- 0.9
  alpha_lab <- 0.3
  
  dsize_low <- 0.1
  wsize_low <- 2.4
  wsize_low <- 1.5
  
  if (ncols_org>6) wsize_low <- 1.0
  
  # Calculate Upper Righer Delta x/y buffers::
  bufs <- NULL
  bufs <- getDeltaMaxXY(data=sdat,datIdx=2, verbose=verbose,vt=vt+1,tc=tc+1)
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Building ggpairs(group={group}, field={field}, gg_mtitle={gg_mtitle})"),"\n", sep='')
  print(sdat)
  
  # group  <- rlang::sym(group)
  gg <- GGally::ggpairs(
    data=sdat, columns=tarColIdxes,
    mapping = ggplot2::aes_(color=as.name(group) ),
    
    # upper = 'blank',
    upper = list(combo = "box_no_facet",
                 continuous = wrap(geomDensity1d_Delta, fdat=pdat, bufs=bufs, group=group, minDelta=minDelta, only=probeType,
                                   alpha=alpha_hih,alpha_lab=alpha_lab, size=dsize_low, wsize=wsize_low,
                                   verbose=verbose,vt=vt+1,tc=tc+1) ),
    
    # diag  = list(continuous='blankDiag'),
    diag  = list(continuous=wrap(diag_density, fdat=ddat, group=group, minPval=minPval, only=probeType, field=field,
                                 alpha=alpha_mid,alpha_lab=alpha_lab, size=dsize_low, wsize=wsize_low,
                                 x_min=0,x_max=1,y_min=0,ticks=3,
                                 verbose=verbose,vt=vt+1,tc=tc+1) ),
    
    # lower = 'blank'
    lower = list(combo = "box_no_facet",
                 continuous = wrap(geomDensity2d_RSquared, fdat=pdat, group=group, only=probeType,
                                   alpha=alpha_low,alpha_lab=alpha_lab, size=dsize_low, wsize=wsize_low,
                                   x_min=0,x_max=1,y_min=0,y_max=1,ticks=3,
                                   verbose=verbose,vt=vt+1,tc=tc+1) )
  )
  gg <- gg +
    theme(panel.grid.major = element_blank()) +
    labs(title=gg_mtitle, subtitle=gg_stitle, caption=gg_ctitle)
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} GGally::ggpairs writing image={gg_pdf}.{RET}"))
  suppressMessages(ggplot2::ggsave(gg_pdf, gg, dpi=dpi) )
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}{RET}"))
  
  gg
}

diag_density <- function(data, mapping, fdat=NULL, group, minPval, only, field='Beta',
                         alpha=0.8,alpha_lab=0.3,size=0.5,wsize=3,
                         x_min=NULL,x_max=NULL,y_min=NULL,y_max=NULL,ticks=3,
                         verbose=0,vt=4,tc=1, ...) {
  funcTag <- 'diag_density'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting; group={group}, minPval={minPval}, only={only}, field={field}.{RET}"))
  # print(mapping)
  # cat("fdat::\n")
  # print(fdat)
  # cat("data::\n")
  # print(head(data))
  
  if (is.null(fdat)) fdat <- data
  # Canonical Reference Samples don't have pvals, so we need to mock an fdat with all zeros::
  if (c(mapping$x) %in% names(fdat)) { fdat <- fdat %>% dplyr::mutate(Pval=GGally::eval_data_col(fdat, mapping$x))
  } else { fdat <- fdat %>% dplyr::mutate(Pval=0) }
  data <- data %>% dplyr::mutate(!!field :=GGally::eval_data_col(data, mapping$x))
  # data <- data %>% dplyr::mutate(Beta=GGally::eval_data_col(data, mapping$x))
  
  sum_str <- getSummaryLabel(fdat, group=group, field='Pval', minVal=minPval)
  
  if (!is.null(only)) {
    group <- rlang::sym(group)
    data <- data %>% dplyr::filter(str_starts(!!group, only))
  }
  
  field <- rlang::sym(field)
  if (is.null(x_min)) x_min <- min(data[[field]], na.rm=TRUE)
  if (is.null(x_max)) x_max <- max(data[[field]], na.rm=TRUE)
  if (is.null(y_min)) y_min <- 0
  if (is.null(y_max)) {
    y_idx  <- which.max(density(data[[field]], na.rm=TRUE)$y)
    y_max  <- density(data[[field]], na.rm=TRUE)$y[y_idx]
    y_max  <- y_max + (y_max*0.3)
  }
  
  if (verbose>=vt+1) cat(glue::glue("[{funcTag}]:{tabsStr} x_min={x_min}, x_max={x_max}.{RET}"))
  gg <- ggplot(data=data, mapping=mapping) +
    geom_density(..., alpha=alpha, aes(fill=!!group) ) +
    scale_x_continuous(breaks = round(seq(x_min,x_max, length.out=ticks),1), limits=c(x_min, x_max)) +
    scale_y_continuous(breaks = round(seq(y_min,y_max, length.out=ticks),1), limits=c(y_min, y_max)) +
    ggplot2::geom_label(
      data = data.frame(
        xlabel = x_min,
        ylabel = y_max,
        lab = sum_str),
      mapping = ggplot2::aes(x = xlabel,
                             y = ylabel,
                             label = lab),
      hjust = 'left', vjust = 'top',
      size = wsize, fontface = "bold",
      alpha=0.3, family='mono',
      show.legend=TRUE,
      inherit.aes = FALSE # do not inherit anything from the ...
    )
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done...{RET}{RET}"))
  
  gg
}

geomDensity1d_Delta <- function(data, mapping, fdat=NULL, bufs=NULL, group, minDelta, only=NULL,
                                alpha=0.8,alpha_lab=0.3,size=0.5,wsize=3, ticks=3,
                                verbose=0,vt=4,tc=1) {
  funcTag <- 'geomDensity1d_Delta'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  data <- data %>% dplyr::mutate(Delta=GGally::eval_data_col(data, mapping$x) - GGally::eval_data_col(data, mapping$y)) %>% 
    tibble::as_tibble()
  if (is.null(fdat)) { 
    fdat <- data
  } else {
    fdat <- fdat %>% dplyr::mutate(Delta=GGally::eval_data_col(fdat, mapping$x) - GGally::eval_data_col(fdat, mapping$y))
  }
  
  sum_str <- getSummaryLabel(fdat, group=group, field='Delta', minVal=minDelta, rm.tot=TRUE)
  
  x_max <- 0
  x_min <- 0
  y_max <- 0
  y_min <- 0
  if (is.null(bufs)) {
    x_min <- min(data$Delta, na.rm=TRUE)
    x_max <- max(data$Delta, na.rm=TRUE)
    
    y_idx  <- which.max(density(data$Delta, na.rm=TRUE)$y)
    y_max  <- density(data$Delta, na.rm=TRUE)$y[y_idx]
  } else {
    x_max <- bufs$x_buf
    y_max <- bufs$y_buf
    x_min <- -x_max
    y_min <- 0
  }
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} x=[{x_min}, {x_max}], y=[{y_min},{y_max}]{RET}"))
  
  str_df <- data.frame(
    xlabel = as.double(x_min),
    ylabel = as.double(y_max),
    lab = sum_str)
  
  if (!is.null(only)) {
    group <- rlang::sym(group)
    data <- data %>% dplyr::filter(str_starts(!!group, only))
  }
  
  gg <- ggplot(data=data) +
    geom_density(aes(x=Delta, color=!!group, fill=!!group) ) +
    geom_vline(xintercept = -minDelta, color="red",  linetype="dotted") +
    geom_vline(xintercept =  minDelta, color="blue", linetype="dotted") +
    scale_x_continuous(breaks = round(seq(x_min, x_max, length.out=ticks),1), limits=c(x_min, x_max)) +
    scale_y_continuous(breaks = round(seq(y_min, y_max, length.out=ticks),1), limits=c(y_min, y_max)) +
    
    ggplot2::geom_label(data=str_df,
                        mapping = ggplot2::aes(x = xlabel, y = ylabel, label = lab),
                        hjust = "left", vjust = "top",
                        size = wsize, fontface = "bold",
                        alpha=alpha_lab, family='mono',
                        inherit.aes = FALSE # do not inherit anything from the ...
    )
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}{RET}"))
  
  gg
}

geomDensity2d_RSquared <- function(data, mapping, fdat=NULL, group, only=NULL,
                                   alpha=0.8,alpha_lab=0.3,size=0.5,wsize=3, 
                                   x_min=0,x_max=1,y_min=0,y_max=1,ticks=3,
                                   verbose=0,vt=4,tc=1) {
  funcTag <- 'geomDensity2d_RSquared'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  if (is.null(fdat)) fdat <- data
  
  cor_tib <- NULL
  sp_tib <- fdat %>% split(data[[group]])
  for (sp in names(sp_tib)) {
    x <- GGally::eval_data_col(sp_tib[[sp]], mapping$x)
    y <- GGally::eval_data_col(sp_tib[[sp]], mapping$y)
    cor <- round(cor(x, y,  method='pearson', use='pairwise.complete.obs'),3)
    cor_tib <- cor_tib %>% dplyr::bind_rows(tibble(!!group := !!sp, R2 = cor) )
  }
  sum_str <- getSummaryLabel(cor_tib, group=group, field='R2', minVal=NULL)
  
  if (!is.null(only)) {
    group <- rlang::sym(group)
    data <- data %>% dplyr::filter(str_starts(!!group, only))
  }
  
  gg <- ggplot(data=data, mapping=mapping) +
    geom_point(size=size) +
    geom_density2d(alpha=alpha) +
    scale_x_continuous(breaks = round(seq(x_min, x_max, length.out=ticks),1), limits=c(x_min, x_max)) +
    scale_y_continuous(breaks = round(seq(y_min, y_max, length.out=ticks),1), limits=c(y_min, y_max)) +
    
    ggplot2::geom_label(
      data = data.frame(
        xlabel = min(x, na.rm = TRUE),
        ylabel = max(y, na.rm = TRUE),
        lab = sum_str),
      mapping = ggplot2::aes(x = xlabel,
                             y = ylabel, 
                             label = lab),
      hjust = 'left', vjust = 'top',
      size = wsize, fontface = "bold",
      alpha=alpha_lab, family='mono',
      inherit.aes = FALSE # do not inherit anything from the ...
    )
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}{RET}"))
  
  gg
}

getDeltaMaxXY = function(data, datIdx=2, verbose=0,vt=4,tc=1) {
  funcTag <- 'getDeltaMaxXY'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  x_buf_fin <- NULL
  y_buf_fin <- NULL
  
  ncols <- base::ncol(data)
  for (idxA in c(datIdx:ncols)) {
    for (idxB in c(datIdx:ncols)) {
      if (idxA<idxB) {
        # Normal Method without 'all_of(idxA)'
        # d <- data %>% dplyr::select(idxA, idxB) %>% purrr::set_names('V1', 'V2')
        # d <- data %>% dplyr::select(!!idxA, !!idxB) %>% purrr::set_names('V1', 'V2')
        d <- data %>% dplyr::select(all_of(idxA), all_of(idxB)) %>% purrr::set_names('V1', 'V2')
        d <- d %>% mutate(Delta=V1-V2)
        
        x_min <- min(d$Delta, na.rm=TRUE)
        x_max <- max(d$Delta, na.rm=TRUE)
        x_buf <- max(abs(x_min), abs(x_max) )
        x_mid <- 0
        
        y_idx <- which.max(density(d$Delta, na.rm=TRUE)$y)
        y_max <- density(d$Delta, na.rm=TRUE)$y[y_idx]
        
        if (is.null(x_buf_fin) || x_buf_fin<x_buf) x_buf_fin <- x_buf
        if (is.null(y_buf_fin) || y_buf_fin<y_max) y_buf_fin <- y_max
        
        if (verbose>=vt)
          cat(glue::glue("[{funcTag}]:{tabsStr} idxA={idxA}, idxB={idxB}, x_buf={x_buf}, y_idx/y_max={y_idx}/{y_max},",
                         " x_buf_fin={x_buf_fin}, y_buf_fin={y_buf_fin}{RET}") )
      }
    }
  }
  ret <- NULL
  ret$x_buf <- x_buf_fin
  ret$y_buf <- y_buf_fin
  if (verbose>=vt)
    cat(glue::glue("[{funcTag}]:{tabsStr} Done. x_buf_fin={x_buf_fin}, y_buf_fin={y_buf_fin}{RET}{RET}") )
  
  ret
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                             General Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
sampleGroup_n = function(tib, n, field,
                         verbose=0,vt=4,tc=1) {
  funcTag <- 'sampleGroup_n'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  tib_sp <- tib %>% split(.[[field]])
  
  tib_sub <- NULL
  for (sp in names(tib_sp)) {
    min_n <- tib_sp[[sp]] %>% base::nrow() %>% min(n)
    
    tib_sub <- tib_sp[[sp]] %>% dplyr::sample_n(min_n) %>% dplyr::bind_rows(tib_sub)
  }
  
  tib_sub
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                    General Text Formatting Functions::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

number_as_commaK <- function(value, cutoff=10000) {
  if (!is.numeric(value)) {
    cat("[ERROR]: Not-Numeric=",value,"\n\n", sep='')
  }
  stopifnot(is.numeric(value))
  
  s <- ''
  if (value>=cutoff) {
    s <- 'k'
    value <- as.integer(value/1000)
  }
  paste0(format(value, big.mark = ",", scientific = FALSE, digits = 22),s)
}

getSummaryLabel = function(tib, pad=' ', group, field, minVal=NULL, rm.tot=FALSE,
                           verbose=0,vt=4,tc=1) {
  funcTag <- 'getSummaryLabel'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  group <- rlang::sym(group)
  field <- rlang::sym(field)
  sum_tib <- tib %>% dplyr::select(!!group,!!field) 
  
  if (!is.null(minVal)) {
    sum_tib <- sum_tib %>% group_by(!!group) %>% 
      dplyr::summarise(TotalCnt=n(),
                       Pass_Cnt=count(!!field <= minVal, na.rm=TRUE),
                       Fail_Cnt=count(!!field >  minVal, na.rm=TRUE),
                       Nana_Cnt=count(is.na(!!field)),
                       Pass=round(100*Pass_Cnt/TotalCnt,1),
                       Fail=number_as_commaK(count(!!field >  minVal, na.rm=TRUE)),
                       Total=number_as_commaK(n())) %>% 
      dplyr::select(!!group,Pass,Fail,Total) %>%
      # dplyr::mutate(Pass=paste0(Pass,'%')) %>%
      dplyr::rename('Pass%'=Pass)
    if (rm.tot) sum_tib <- sum_tib %>% dplyr::select(-Total)
  }
  sum_tib <- sum_tib %>% dplyr::mutate_all(as.character) %>%
    dplyr::rename(type=!!group)
  key_ord <- names(sum_tib)
  key_tib <- tibble(key_ord) %>% dplyr::mutate(val=key_ord) %>% spread(key_ord, val)
  sum_tib <- key_tib %>% dplyr::bind_rows(sum_tib) %>% dplyr::select(key_ord)
  cnt_tib <- sum_tib %>% lapply(function(x) max(nchar(x))) %>% dplyr::bind_cols() %>% 
    dplyr::select(key_ord)
  
  del <- ','
  lab_str <- ''
  nrows <- base::nrow(sum_tib)
  for (ii in c(1:nrows)) {
    ncols <- sum_tib[ii,] %>% base::ncol()
    
    cur_str <- ''
    for (jj in c(1:ncols)) {
      len <- as.integer(cnt_tib[jj])
      str <- as.character(head(sum_tib[ii,jj],n=1) )
      if (str=='100' || str=='R2')
        val <- stringr::str_pad(string=str, width=len, side='right', pad=pad)
      else
        val <- stringr::str_pad(string=str, width=len, side='left', pad=pad)
      
      cur_str <- paste(cur_str,val, sep=del)
    }
    cur_str <- cur_str %>% 
      stringr::str_remove(del) %>%
      stringr::str_replace(del,'=')
    
    if (ncols>2) {
      cur_str <- cur_str %>%
        stringr::str_replace(del,' [') %>%
        stringr::str_replace(del,', ')
      cur_str <- paste0(cur_str,']')
    }
    lab_str <- paste(lab_str,cur_str,'\n', sep='')
  }
  lab_str <- lab_str %>% stringr::str_remove('\\n$')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}{RET}"))
  
  lab_str
}



# End of file
