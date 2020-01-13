
suppressWarnings(suppressPackageStartupMessages(require("tidyverse")) )
suppressWarnings(suppressPackageStartupMessages(require("stringr")) )
suppressWarnings(suppressPackageStartupMessages(require("glue")) )
suppressWarnings(suppressPackageStartupMessages(require("scales")) )
suppressWarnings(suppressPackageStartupMessages(require("matrixStats")) )
suppressWarnings(suppressPackageStartupMessages(require("grid")) )

# Plotting Extras
suppressWarnings(suppressPackageStartupMessages(require("GGally")) )
suppressWarnings(suppressPackageStartupMessages(require("hexbin")) )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                        Sample Sheet I/O Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

loadAutoSampleSheets = function(dir, platform, manifest, suffix='AutoSampleSheet.csv.gz', 
                                addSampleName=FALSE, addPaths=FALSE,
                                
                                pvalDetectFlag=TRUE, pvalDetectMinKey='CG_NDI_negs_pval_PassPerc', pvalDetectMinVal=96,
                                
                                flagSampleDetect=TRUE, filterRef=FALSE,
                                dbMin=90, r2Min=0.9,
                                dbKey='AutoSample_dB_Key', r2Key='AutoSample_R2_Key',
                                dbVal='AutoSample_dB_Val', r2Val='AutoSample_R2_Val',
                                verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'loadAutoSampleSheets'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting, dir={dir}.{RET}"))

  stime <- system.time({
    dat <- NULL
    
    dbKey <- dbKey %>% rlang::sym()
    dbVal <- dbVal %>% rlang::sym()
    r2Key <- r2Key %>% rlang::sym()
    r2Val <- r2Val %>% rlang::sym()
    
    pvalDetectMinKey <- pvalDetectMinKey %>% rlang::sym()
    # pvalDetectMinVal <- pvalDetectMinVal %>% rlang::sym()

    pattern <- paste(platform,manifest,suffix, sep='.')
    auto_ss_list <- list.files(dir, pattern=pattern, recursive=TRUE, full.names=TRUE)
    auto_ss_llen <- auto_ss_list %>% length()
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting, SampleSheetCount={auto_ss_llen}{RET}"))
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #  Load Samples::
    auto_ss_tibs <- suppressMessages(suppressWarnings(lapply(auto_ss_list, readr::read_csv) )) %>% 
      dplyr::bind_rows() # %>% 
    #  addBeadPoolToSampleSheet(field='CG_Loci_Count', verbose=verbose,vt=vt+1,tc=tc+1)
    auto_ss_tlen <- base::nrow(auto_ss_tibs)
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} SampleSheetNrows(bPool)={auto_ss_tlen}{RET}"))
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #  Add General Sample Name Field::
    if (addSampleName) {
      auto_ss_tibs <- auto_ss_tibs %>% dplyr::mutate(Sample_Name=!!dbKey) %>%
        dplyr::select(Sample_Name, dplyr::everything())
    }
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #  Flag Probe Detected (pval)::
    if (pvalDetectFlag && !is.null(pvalDetectMinKey) && !is.null(pvalDetectMinVal)) {
      fail_tag <- paste0("Failed<",pvalDetectMinVal)
      pass_tag <- paste0("Passed>=",pvalDetectMinVal)
      
      auto_ss_tibs <- auto_ss_tibs %>% dplyr::mutate(detectPval=case_when(
        !!pvalDetectMinKey < !!pvalDetectMinVal ~ fail_tag,
        TRUE ~ pass_tag)
      )
    }
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #  Flag Auto-Detected Samples::
    if (flagSampleDetect) {
      fail_r2  <- paste0('Failed_r2<',r2Min)
      fail_db  <- paste0('Failed_db<',dbMin)
      fail_mt  <- paste0('Failed_r21-db')
      pass_tag <- 'Passed'
      
      auto_ss_tibs <- auto_ss_tibs %>% dplyr::mutate(detectedSample=case_when(
        !!r2Val < r2Min ~ fail_r2,
        !!dbVal < dbMin ~ fail_db,
        !!dbKey != !!r2Key ~ fail_mt,
        TRUE ~ pass_tag)
      )
    }
    
    if (filterRef) {
      auto_ss_tibs <- auto_ss_tibs %>% dplyr::filter() %>% dplyr::filter(!!dbVal >= dbMin) %>% dplyr::filter(!!r2Val >= r2Min)
      auto_ss_flen <- base::nrow(auto_ss_tibs)
      if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} SampleSheetNrows(paths)={auto_ss_flen}{RET}"))

      # Remove Unidentifiable
      # if (rmOdd) auto_ss_tibs <- auto_ss_tibs %>% dplyr::filter(Bead_Pool!='Odd')
    }
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #  Add Paths::
    if (addPaths) {
      auto_ss_tibs <- auto_ss_tibs %>%
        addPathsToSampleSheet(dir=dir, platform=platform, manifest=manifest, 
                              field='Calls_Path', suffix='call.csv.gz$', verbose=verbose)
      
      #  # addPathsToSampleSheet(dir=dir, platform=platform, field='ProbesI_CSV_Path',  suffix='both.ProbesI.tib.csv.gz$', verbose=verbose) %>%
      #  # addPathsToSampleSheet(dir=dir, platform=platform, field='ProbesII_CSV_Path', suffix='both.ProbesII.tib.csv.gz$', verbose=verbose) %>%
      #   addPathsToSampleSheet(dir=dir, platform=platform, field='Probes_RDS_Path',   suffix='both.Probes.tib.rds$', verbose=verbose)
      auto_ss_tlen <- base::nrow(auto_ss_tibs)
      if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} SampleSheetNrows(paths)={auto_ss_tlen}{RET}"))
    }
    
    dat <- auto_ss_tibs
  })
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}{RET}"))
  if (!is.null(tt)) tt$addTime(stime,funcTag)

  dat
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                     Sample Sheet Manipulation Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

getUniqueFields = function(ss, keys, verbose=0,vt=1,tc=1) {
  funcTag <- 'getUniqueFields'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  # Remove Variables with a single value
  uniq_exp_tib <- ss %>% dplyr::select(keys) %>% 
    # dplyr::summarise_each(list(n_distinct) ) %>% 
    dplyr::summarise_all(list(n_distinct) ) %>% 
    tidyr::gather() %>% dplyr::filter(value>1) %>%
    dplyr::arrange(value) %>% 
    tidyr::spread(key, value)
  
  uniq_exp_keys <- uniq_exp_tib %>% names()
  
  # Special Treatment for ChipFormat (put in front)
  field <- 'ChipFormat'  
  isPresent <- grep(field, uniq_exp_keys) %>% length() > 0
  if (isPresent)
    uniq_exp_tib <- uniq_exp_tib %>% dplyr::select(field, everything())
  
  # Special Treatment for Bead_Pool (put in last)
  field <- 'Bead_Pool'  
  isPresent <- grep(field, uniq_exp_keys) %>% length() > 0
  if (isPresent)
    uniq_exp_tib <- uniq_exp_tib %>% dplyr::select(-field, everything())
  
  uniq_exp_keys <- uniq_exp_tib %>% names()
  
  uniq_exp_keys
}

addPathsToSampleSheet = function(ss, dir, platform, manifest, field, suffix, del='.',
                                 verbose=0,vt=1,tc=1) {
  funcTag <- 'addRdsPathsToSampleSheet'
  
  pattern=paste(platform,manifest,suffix, sep=del)
  file_list <- NULL
  file_list <- list.files(dir, pattern=pattern, recursive=TRUE, full.names=TRUE)
  if (length(file_list)==0) file_list <- NULL
  
  # TBD:: Clean up code below!!!
  ss_path_tibs <- tibble::tibble(Sentrix_Name=basename(file_list) %>% 
                                   str_remove(paste('',opt$platform,'.*$', sep=del)), 
                                 !!field := file_list)
  
  tib <- ss_path_tibs %>%
    dplyr::distinct(Sentrix_Name, .keep_all=TRUE) %>% 
    dplyr::inner_join(ss, by="Sentrix_Name")
  
  tib
}

addBeadPoolToSampleSheet = function(ss, field, verbose=0, vt=1, tc=1) {
  funcTag <- 'addBeadPoolToSampleSheet'
  
  field <- field %>% rlang::sym()
  ss <- ss %>% dplyr::mutate(
    Bead_Pool=case_when(!!field >862926+1000 ~ 'EPIC_Plus',
                        !!field >862926-1000 ~ 'EPIC',
                        !!field >612329-1000 ~ 'BP1234',
                        !!field >448696-1000 ~ 'BP123',
                        !!field >148977-1000 ~ 'BP2',
                        !!field >103113-1000 ~ 'Odd',
                        !!field<=103113-1000 ~ 'UNK',
                        TRUE ~ NA_character_) ) %>%
    dplyr::select(Bead_Pool, everything())
  
  ss
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#  Below is the current code from plotByExperiments Jan-7-2019
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #


# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                            Clean Functions::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

selectSamples = function(ss, max=3, spread='mid', sort='CG_inf_negs_pval_PassPerc',
                         verbose=0,vt=4,tc=1,tabsStr='') {
  funcTag <- 'selectSamples'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} max={max}, sort={sort}{RET}"),sep='')
  max=3
  
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

loadSampleRDS_byRep = function(ss, exp=NULL, minPval, parallel=FALSE, verbose=0,vt=3,tc=1) {
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
                            prefix=exp_name, verbose=verbose, vt=vt, tc=tc+1)
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

loadProbeRDS = function(file, minPval=NULL, prefix=NULL, 
                        pvalKey='inf_negs_pval', betaKey='inf_noob_dye_beta',
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

plotBetaMatrix_bySample = function(tib, sample, manifest, minPval, outDir, field='Beta', join='Probe_ID',
                                   max=3, spread='mid', outType='auto', dpi=72, format='png', maxCnt=30000,
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
      ss_exp_tibs[[exp_key]]  <- selectSamples(ss_exp_tibs[[exp_key]], max=max, spread=spread, sort='CG_inf_negs_pval_PassPerc',
                                               verbose=verbose,vt=vt,tc=tc+1)
      if (loadRDS)
        prbs_exp_tibs[[exp_key]] <- loadSampleRDS_byRep(ss=ss_exp_tibs[[exp_key]], exp=exp_key, minPval=minPval, verbose=verbose+9,vt=vt+1,tc=tc+1)
      
      if (loadIDAT)
        idat_exp_tibs[[exp_key]] <- loadSampleIDAT_byRep(ss=ss_exp_tibs[[exp_key]], exp=exp_key, verbose=verbose+9,vt=vt+1,tc=tc+1)
    }
  }
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Loaded ALL Experimental Probe Data!{RET}{RET}"))
  
  # ret <- NULL
  # ret$ss_exp_tibs   <- ss_exp_tibs
  # ret$prbs_exp_tibs <- prbs_exp_tibs
  # # ret$idat_exp_tibs <- idat_exp_tibs
  # return(ret)
  
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
          
          join_prbs <- manifest %>% dplyr::right_join(dplyr::full_join(prbsA,prbsB, by="Probe_ID"), by="Probe_ID") %>%
            rename(Design_Type=DESIGN)
          
          # join_prbs <- joinTibsAndManifest(datA=prbsA,keyA=nameA, datB=prbsB,keyB=nameB, 
          #                                  manifest=manifest, field=field, join=join,
          #                                  verbose=verbose, vt=vt+1, tc=tc+1)
          
          # join_idat <- joinTibsAndManifest(datA=prbsA,keyA=nameA, datB=prbsB,keyB=nameB, 
          #                                  manifest=manifest, field=field, join=join,
          #                                  verbose=verbose, vt=vt+1, tc=tc+1)
          
          field_str <- 'inf_noob_dye_beta'
          detP_str <- 'inf_negs_pval'
          gg <- plotPairs(join_prbs, sample=sample, nameA=nameA, nameB=nameB, outDir=sam_dir,
                          probeType='cg', field=field, field_str=field_str, detp=detP_str, maxCnt=maxCnt, minPval=minPval,
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
  
  # Calculate Upper Righer Delta x/y buffers::
  bufs <- NULL
  bufs <- getDeltaMaxXY(data=sdat,datIdx=2, verbose=verbose,vt=vt+1,tc=tc+1)
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Building ggpairs(group={group}, field={field}, gg_mtitle={gg_mtitle})"),"\n", sep='')
  
  # group  <- rlang::sym(group)
  gg <- GGally::ggpairs(sdat, columns=tarColIdxes,
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
  
  suppressMessages(ggplot2::ggsave(gg_pdf, gg, dpi=dpi) )
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}{RET}"), sep='')
  
  gg
}

diag_density <- function(data, mapping, fdat=NULL, group, minPval, only, field='Beta',
                         alpha=0.8,alpha_lab=0.3,size=0.5,wsize=3,
                         x_min=NULL,x_max=NULL,y_min=NULL,y_max=NULL,ticks=3,
                         verbose=0,vt=4,tc=1, ...) {
  funcTag <- 'diag_density'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  if (is.null(fdat)) fdat <- data
  fdat <- fdat %>% dplyr::mutate(Pval=GGally::eval_data_col(fdat, mapping$x))
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
  
  gg
}

geomDensity1d_Delta <- function(data, mapping, fdat=NULL, bufs=NULL, group, minDelta, only=NULL,
                                alpha=0.8,alpha_lab=0.3,size=0.5,wsize=3, ticks=3,
                                verbose=0,vt=4,tc=1,tabsStr='') {
  funcTag <- 'geomDensity1d_Delta'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"), sep='')
  
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
  if (verbose>=vt+1) cat(glue::glue("[{funcTag}]:{tabsStr} x=[{x_min}, {x_max}], y=[{y_min},{y_max}]{RET}"), sep='')
  
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
  
  gg
}

geomDensity2d_RSquared <- function(data, mapping, fdat=NULL, group, only=NULL,
                                   alpha=0.8,alpha_lab=0.3,size=0.5,wsize=3, 
                                   x_min=0,x_max=1,y_min=0,y_max=1,ticks=3,
                                   verbose=0,vt=4,tc=1,tabsStr='') {
  funcTag <- 'geomDensity2d_RSquared'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
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
  
  gg
}

getDeltaMaxXY = function(data, datIdx=2, verbose=0,vt=4,tc=1,tabsStr='') {
  funcTag <- 'getDeltaMaxXY'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"), sep='')
  
  x_buf_fin <- NULL
  y_buf_fin <- NULL
  
  ncols <- base::ncol(data)
  for (idxA in c(datIdx:ncols)) {
    for (idxB in c(datIdx:ncols)) {
      if (idxA<idxB) {
        d <- data %>% dplyr::select(idxA, idxB) %>% purrr::set_names('V1', 'V2')
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
                         verbose=0,vt=4,tc=1,tabsStr='') {
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
                           verbose=0,vt=4,tc=1,tabsStr='') {
  funcTag <- 'getSummaryLabel'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
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
  
  lab_str
}



# End of file
