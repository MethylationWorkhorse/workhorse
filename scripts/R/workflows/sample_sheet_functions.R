

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                          Sample Sheet Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

suppressWarnings(suppressPackageStartupMessages( base::require("tidyverse") ))
suppressWarnings(suppressPackageStartupMessages( base::require("stringr") ))
suppressWarnings(suppressPackageStartupMessages( base::require("glue") ))

COM <- ","
TAB <- "\t"
RET <- "\n"

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                        Sample Sheet I/O Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

listChipIds = function(dir, verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'listChipIds'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting, dir={dir}.{RET}"))

  tib <- NULL
  dirs <- list.dirs(opt$buildDir, recursive = FALSE)
  for (d in dirs) {
    exp_name <- base::basename(d)
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Experiment={exp_name}, dir={d}.{RET}"))

    chip_paths <- list.dirs(d, recursive = FALSE)
    chip_names <- chip_paths %>% base::basename()
    
    cur_tib <- tibble::tibble(Build_Paths=chip_paths, 
                              Sentrix_Barcode=chip_names,
                              Experiment_Name=exp_name) %>%
      dplyr::filter(Sentrix_Barcode!='shells') %>%
      dplyr::mutate(Sentrix_Barcode=as.double(Sentrix_Barcode))
    
    # cur_tib %>% print()
    tib <- tib %>% dplyr::bind_rows(cur_tib)
  }
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}{RET}"))
  
  tib <- tib %>% dplyr::select(Sentrix_Barcode, Experiment_Name, everything())

  tib
}

loadAutoSampleSheets = function(dir, platform, manifest, suffix='AutoSampleSheet.csv.gz', 
                                addSampleName=FALSE, addPathsCall=FALSE, addPathsSigs=FALSE,
                                pvalDetectFlag=TRUE, pvalDetectMinKey, pvalDetectMinVal=96,
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
    
    if (!is.null(pvalDetectMinKey)) pvalDetectMinKey <- pvalDetectMinKey %>% rlang::sym()
    # pvalDetectMinVal <- pvalDetectMinVal %>% rlang::sym()
    
    pattern <- paste(platform,manifest,suffix, sep='_')
    auto_ss_list <- list.files(dir, pattern=pattern, recursive=TRUE, full.names=TRUE)
    auto_ss_llen <- auto_ss_list %>% length()
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting, SampleSheetCount={auto_ss_llen}, pattern={pattern}.{RET}"))
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #  Load Samples::
    auto_ss_tibs <- suppressMessages(suppressWarnings(lapply(auto_ss_list, readr::read_csv) )) %>% 
      dplyr::bind_rows() # %>% 
    #  addBeadPoolToSampleSheet(field='CG_Loci_Count', verbose=verbose,vt=vt+1,tc=tc+1)
    auto_ss_tlen <- base::nrow(auto_ss_tibs)
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} SampleSheetNrows(bPool)={auto_ss_tlen}{RET}"))
    # print(auto_ss_tibs)
    
    if (addSampleName || pvalDetectFlag || flagSampleDetect || filterRef) {
      if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Adding flags...{RET}"))
      dbKey <- dbKey %>% rlang::sym()
      dbVal <- dbVal %>% rlang::sym()
      r2Key <- r2Key %>% rlang::sym()
      r2Val <- r2Val %>% rlang::sym()
      
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
    }
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #  Add Paths::
    if (addPathsCall) {
      auto_ss_tibs <- auto_ss_tibs %>%
        addPathsToSampleSheet(dir=dir, platform=platform, manifest=manifest, 
                              field='Calls_Path', suffix='calls.csv.gz$', verbose=verbose)
      
      auto_ss_tlen <- base::nrow(auto_ss_tibs)
      if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} SampleSheetNrows(paths)={auto_ss_tlen}{RET}"))
    }
    if (addPathsSigs) {
      auto_ss_tibs <- auto_ss_tibs %>%
        addPathsToSampleSheet(dir=dir, platform=platform, manifest=manifest, 
                              field='Sigs_Path', suffix='sigs.csv.gz$', verbose=verbose)
      
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
  
  # Special Treatment for Chip_Format (put in front)
  field <- 'Chip_Format'  
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
#                             Manifest I/O::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
loadManifestSource = function(file, verbose=0,vt=4,tc=1,tt=NULL) {
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
                             verbose=0,vt=4,tc=1,tt=NULL) {
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

# End of file
