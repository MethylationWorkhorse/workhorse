

suppressWarnings(suppressPackageStartupMessages(require("tidyverse")) )
suppressWarnings(suppressPackageStartupMessages(require("stringr")) )
suppressWarnings(suppressPackageStartupMessages(require("glue")) )
suppressWarnings(suppressPackageStartupMessages(require("scales")) )
suppressWarnings(suppressPackageStartupMessages(require("matrixStats")) )

# Plotting Extras
suppressWarnings(suppressPackageStartupMessages(require("GGally")) )
suppressWarnings(suppressPackageStartupMessages(require("hexbin")) )

# Parallel Computing Packages
suppressWarnings(suppressPackageStartupMessages(require("doParallel")) )

doParallel::registerDoParallel()
num_cores   <- detectCores()
num_workers <- getDoParWorkers()

COM <- ","
TAB <- "\t"
RET <- "\n"

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                             Time Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

formatSysTime = function(stime, key, verbose=0, vt=1, tc=1) {
  funcTag <- 'formatTime'
  
  ctime <- stime %>% tibble::enframe() %>% 
    dplyr::mutate(Method=key, name=stringr::str_replace_all(name,'\\.','_')) %>% 
    tidyr::spread('name', 'value')
  
  ctime
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                          Sample Sheet Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

loadAutoSampleSheets = function(dir, platform, suffix='both.SampleSheet.csv.gz', rmOdd=TRUE,
                                verbose=0, vt=4, tc=1) {
  funcTag <- 'loadAutoSampleSheets'
  
  pattern <- paste(platform,suffix, sep='_')
  auto_ss_list <- list.files(dir, pattern=pattern, recursive=TRUE, full.names=TRUE)
  auto_ss_tibs <- suppressMessages(suppressWarnings(lapply(auto_ss_list, readr::read_csv) )) %>% 
    dplyr::bind_rows() %>% 
    addBeadPoolToSampleSheet(verbose=verbose) %>%
    addRdsPathsToSampleSheet(dir=dir, platform=platform, verbose=verbose) %>%
    dplyr::mutate(Sample_Name=Auto_Sample_BETA_bDB_Key) %>%
    dplyr::select(Sample_Name, dplyr::everything())
  
  # Remove Unidentifiable
  if (rmOdd) auto_ss_tibs <- auto_ss_tibs %>% dplyr::filter(Bead_Pool!='Odd')
  
  if (verbose>vt) {
    auto_ss_tibs %>% dplyr::group_by(iscan_Extract_Year,Bead_Pool,ChipFormat) %>% 
      dplyr::summarise(Count=n()) %>% print()
    auto_ss_tibs %>% dplyr::group_by(iscan_Extract_Year,Bead_Pool,ChipFormat,Sentrix_Barcode) %>% 
      dplyr::summarise(Count=n()) %>% print()
  }
  
  auto_ss_tibs
}

addBeadPoolToSampleSheet = function(ss, verbose=0, vt=1, tc=1) {
  funcTag <- 'addBeadPoolToSampleSheet'
  
  ss <- ss %>% dplyr::mutate(
    Bead_Pool=case_when(Total_Count_CG >862926+1000 ~ 'EPIC_Plus',
                        Total_Count_CG >862926-1000 ~ 'EPIC',
                        Total_Count_CG >612329-1000 ~ 'BP1234',
                        Total_Count_CG >448696-1000 ~ 'BP123',
                        Total_Count_CG >148977-1000 ~ 'BP2',
                        Total_Count_CG >103113-1000 ~ 'Odd',
                        Total_Count_CG<=103113-1000 ~ 'UNK',
                        TRUE ~ NA_character_) ) %>%
    dplyr::select(Bead_Pool, everything())
  
  ss
}

addRdsPathsToSampleSheet = function(ss, dir, platform, suffix='both.Probes.tib.rds$',
                                    verbose=0, vt=1, tc=1) {
  funcTag <- 'addRdsPathsToSampleSheet'
  
  pattern=paste(platform,suffix, sep='_')
  rds_list <- list.files(dir, pattern=pattern, recursive=TRUE, full.names=TRUE)
  ss_path_tibs <- tibble::tibble(Sentrix_Name=basename(rds_list) %>% str_remove(paste('',opt$platform,'.*$', sep='_')), 
                                 Probes_Path=rds_list) %>% 
    dplyr::distinct(Sentrix_Name, .keep_all=TRUE) %>% 
    dplyr::inner_join(ss, by="Sentrix_Name")
  
  ss_path_tibs
}

loadSampleRDS = function(ss, max, parallel=TRUE, verbose=0,vt=3,tc=1,tabsStr='') {
  funcTag <- 'loadSampleRDS'
  for (i in c(1:tc)) tabsStr <- paste0(tabsStr, TAB)

  cat(glue::glue("[{funcTag}]:{tabsStr} Starting(isParallel={parallel})..."),"\n", sep='')
  
  ss_split <- ss %>% head(n=max) %>% split(.$Sentrix_Name)
  
  cur_tibs <- list()
  
  if (parallel) {
    cur_tibs <- foreach (prefix=names(ss_split), .inorder=T, .final = function(x) setNames(x, names(ss_split))) %dopar% {
      prb_tib <- loadProbeTib(ss_split[[prefix]]$Probes_Path, minPval=opt$minDetpNegs)
      match_col <- names(prb_tib)[1]
      exp_name  <- head(ss_split[[prefix]]$Experiment_Key, n=1)
      prb_names <- paste(exp_name, names(prb_tib), sep='.')
      prb_names[1] <- match_col
      prb_tib <- prb_tib %>% purrr::set_names(prb_names)
      
      prb_tib
    }
  } else {
    for (prefix in names(ss_split)) {
      prb_tib <- loadProbeTib(ss_split[[prefix]]$Probes_Path, minPval=opt$minDetpNegs)
      match_col <- names(prb_tib)[1]
      exp_name  <- head(ss_split[[prefix]]$Experiment_Key, n=1)
      prb_names <- paste(exp_name, names(prb_tib), sep='.')
      prb_names[1] <- match_col
      prb_tib <- prb_tib %>% purrr::set_names(prb_names)
      
      cur_tibs[[prefix]] <- prb_tib
    }
  }
  sample_cnt <- names(cur_tibs) %>% length()
  join_tib <- NULL
  for (sentrix_name in names(cur_tibs)) {
    join_tib <- joinTibs(cur_tibs[[sentrix_name]], join_tib, cleanPrefix=TRUE)
  }
  cat(glue::glue("[{funcTag}]:{tabsStr} SampleCount={sample_cnt}, Done."),"\n\n", sep='')
  
  
  join_tib
}












# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                       OLD Sample Sheet Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

getManifestAnnotation = function(opt, verbose=0, vt=1) {
  funcTag <- 'getManifestAnnotation'
  
  man.csv  <- '/Users/bbarnes/Documents/Projects/manifests/methylation/MethylationEPIC_v-1-0_B4.csv.gz'
  man.tib <- base::suppressMessages(base::suppressWarnings(readr::read_csv(man.csv, skip=7)))
  ctl.idx <- grep('Controls', man.tib$IlmnID)
  
  mana.tib <- man.tib %>% head(n=ctl.idx-1) %>% dplyr::rename(Probe_ID=IlmnID) %>% dplyr::arrange(Probe_ID)
  
  subset.len <- 1000
  subset.len <- 865918
  mana.core.tib <- mana.tib %>% head(n=subset.len) %>%
    dplyr::select(Probe_ID, UCSC_RefGene_Name, UCSC_RefGene_Group,
                  UCSC_CpG_Islands_Name, Relation_to_UCSC_CpG_Island,
                  Phantom5_Enhancers, DMR, 
                  DNase_Hypersensitivity_NAME, 
                  OpenChromatin_NAME, 
                  TFBS_NAME) %>%
    dplyr::mutate(UCSC_RefGene_Name=stringr::str_replace(UCSC_RefGene_Name, ';.*$', ''), 
                  UCSC_RefGene_Group=stringr::str_replace(UCSC_RefGene_Group, ';.*$', '')) %>%
    dplyr::rename(Gene_Name=UCSC_RefGene_Name, Gene_Group=UCSC_RefGene_Group,
                  Islands_Name=UCSC_CpG_Islands_Name, Islands_Group=Relation_to_UCSC_CpG_Island,
                  Phantom5=Phantom5_Enhancers, DNaseHype=DNase_Hypersensitivity_NAME,
                  OpenChrom=OpenChromatin_NAME, TFBS=TFBS_NAME) %>%
    dplyr::mutate(Gene_Name=str_replace_all(Gene_Name,'_','-')) %>% 
    dplyr::mutate(Gene_Group=str_replace_all(Gene_Group,"'",'-')) %>% 
    dplyr::mutate(Islands_Group=str_replace_all(Islands_Group,'_','-')) %>% 
    tidyr::unite(Gene, Gene_Name,Gene_Group, sep='_') %>%
    tidyr::unite(Island, Islands_Name,Islands_Group, sep='_') %>%
    dplyr::mutate(Gene=na_if(Gene, 'NA_NA')) %>%
    dplyr::mutate(Island=na_if(Island, 'NA_NA'))
  
  mana.core.tib
}

formatSampleSheet_MAGS = function(file=NULL, verbose=0, vt=1) {
  funcTag <- 'formatSampleSheet_FFPE'
  
  if (is.null(file)) file <- '/Users/bbarnes/Documents/Projects/workhorse/sampleSheets/TableControl_20190819_EPIC123_8x1_MagPrep_exp_cleanHeader.csv'
  stopifnot(file.exists(file))
  raw_tib <- suppressMessages(suppressWarnings(readr::read_csv(file) ))
  #  [1] "Sentrix_Name"          "Decoder_Pool_User"     "Bead_Pool_User"        "Chip_Format_User"      "Sample_Name_User"      "Sample_Type_User"     
  #  [7] "Sample_Prep_User"      "Sample_Plate_User"     "Sentrix_ID_User"       "Sentrix_Position_User"
  tib <- raw_tib %>% dplyr::mutate(Sentrix_Name=paste(Sentrix_Barcode,Sample_Section, sep='_'),
                                   Decoder_Pool_User='BETA', 
                                   Bead_Pool_User='EPIC-B4',
                                   Chip_Format_User='8x1',
                                   Sample_Name_User=stringr::str_to_upper(stringr::str_remove(Sample_ID, '_[0-9]+$')),
                                   Sample_Type_User=case_when(
                                     Sample_Name_User=='U' || Sample_Name_User=='H' || Sample_Name_User=='M' ~ 'REFERENCE',
                                     TRUE ~ 'CELL_LINE'
                                   ),
                                   Sample_Prep_User=stringr::str_replace_all(Workflow,' ','_'),
                                   Sample_Plate_User=NA,
                                   Sample_Workflow_User=Sample_Plate,
                                   Sentrix_ID_User=Sentrix_Barcode,
                                   Sentrix_Position_User=Sample_Section
  ) %>%
    dplyr::group_by(Sample_Name_User) %>%
    dplyr::mutate(Sample_Rep=row_number()) %>% 
    dplyr::ungroup() %>%
    dplyr::select(Sentrix_Name,Decoder_Pool_User,Bead_Pool_User,Chip_Format_User,Sample_Name_User,
                  Sample_Type_User,Sample_Prep_User,Sample_Plate_User,Sample_Workflow_User,Sentrix_ID_User,Sentrix_Position_User)
  
  tib
}

formatSampleSheet_FFPE = function(file=NULL, verbose=0, vt=1) {
  funcTag <- 'formatSampleSheet_FFPE'
  
  if (is.null(file)) file <- '/Users/bbarnes/Documents/Projects/workhorse/sampleSheets/TableControl_8X1_EPICFFPE_Samplesheet_RA_noHeader.csv'
  stopifnot(file.exists(file))
  raw_tib <- suppressMessages(suppressWarnings(readr::read_csv(file) ))
  #  [1] "Sentrix_Name"          "Decoder_Pool_User"     "Bead_Pool_User"        "Chip_Format_User"      "Sample_Name_User"      "Sample_Type_User"     
  #  [7] "Sample_Prep_User"      "Sample_Plate_User"     "Sentrix_ID_User"       "Sentrix_Position_User"
  tib <- raw_tib %>% dplyr::mutate(Sentrix_Name=paste(Sentrix_ID,Sentrix_Position, sep='_'),
                                   Decoder_Pool_User='BETA', 
                                   Bead_Pool_User='EPIC-B4',
                                   Chip_Format_User='8x1',
                                   Sample_Name_User=stringr::str_to_upper(stringr::str_remove(Sample_Name, '_[0-9]+$')),
                                   Sample_Type_User=case_when(
                                     Sample_Name_User=='U' || Sample_Name_User=='H' || Sample_Name_User=='M' ~ 'TITRATION',
                                     TRUE ~ 'CELL_LINE'
                                   ),
                                   Sample_Prep_User='FFPE',
                                   Sample_Plate_User=NA,
                                   Sample_Workflow_User=Sample_Plate,
                                   Sentrix_ID_User=Sentrix_ID,
                                   Sentrix_Position_User=Sentrix_Position
  ) %>%
    dplyr::group_by(Sample_Name_User) %>%
    dplyr::mutate(Sample_Rep=row_number()) %>% 
    dplyr::ungroup() %>%
    dplyr::select(Sentrix_Name,Decoder_Pool_User,Bead_Pool_User,Chip_Format_User,Sample_Name_User,
                  Sample_Type_User,Sample_Prep_User,Sample_Plate_User,Sample_Workflow_User,Sentrix_ID_User,Sentrix_Position_User)
  
  tib
}

addNanoGramInputs = function(tib, opt, file=NULL, verbose=0, vt=4) {
  funcTag <- 'addNanoGramInputs'
  
  # Nano-gram input is not added yet!!!
  #
  if (is.null(file)) file <- file.path(opt$sesDir, 'sampleSheets/official/nano_gram_input.SampleSheet.tsv')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]: Starting, nanoFile={file}"),"\n", sep='')
  stopifnot(file.exists(file))
  
  ss_inp_tib <- suppressMessages(suppressWarnings(readr::read_tsv(file) )) %>%
    dplyr::arrange(Sentrix_Name) %>% dplyr::rename(Sentrix_Barcode=Sentrix_Name)
  
  if (verbose>=vt) print(ss_inp_tib)
  ss_nrows <- ss_inp_tib %>% base::nrow()
  if (verbose>=vt) cat(glue::glue("[{funcTag}]: ss_nrows={ss_nrows}"),"\n", sep='')
  
  ret_tib <- ss_inp_tib %>% dplyr::full_join(tib,  by='Sentrix_Barcode') %>%
    dplyr::select(Sentrix_Name, Sentrix_Barcode, Sentrix_Poscode, everything()) %>% 
    dplyr::rename(Sample_Input=InputNg) %>% 
    dplyr::mutate(Sample_Input=as.integer(Sample_Input),
                  Bead_Pool=case_when(Total_Count_CG >862926-1000 ~ 'BP123456',
                                      Total_Count_CG >612329-1000 ~ 'BP1234',
                                      Total_Count_CG >448696-1000 ~ 'BP123',
                                      Total_Count_CG >148977-1000 ~ 'BP2',
                                      Total_Count_CG >103113-1000 ~ 'Odd',
                                      Total_Count_CG<=103113-1000 ~ 'UNK',
                                      TRUE ~ NA_character_) )
  if (verbose>=vt+4) print(ret_tib)
  
  ret_tib <- ret_tib %>% dplyr::mutate(Sample_Input=ifelse(is.na(Sample_Input), 1000, Sample_Input),
                                       Decoder_Pool_User=ifelse(is.na(Decoder_Pool_User), 'BETA', Decoder_Pool_User))
  if (verbose>=vt+4) print(ret_tib)
  ret_tib <- ret_tib %>% 
    dplyr::mutate(Chip_Format_User=case_when(is.na(Chip_Format_User) ~ ChipFormat, TRUE ~ Chip_Format_User),
                  Sample_Name_User=case_when(is.na(Sample_Name_User) ~ Auto_Sample_BETA_bDB_Key, TRUE ~ Sample_Name_User),
                  Sample_Type_User=case_when(Sample_Name_User=='U' ~ 'TITRATION',
                                             Sample_Name_User=='H' ~ 'TITRATION',
                                             Sample_Name_User=='M' ~ 'TITRATION', 
                                             Sample_Name_User=='MCF7'   ~ 'CELL_LINE', 
                                             Sample_Name_User=='RAJI'   ~ 'CELL_LINE', 
                                             Sample_Name_User=='JURKAT' ~ 'CELL_LINE', 
                                             Sample_Name_User=='HELA'   ~ 'CELL_LINE', 
                                             Sample_Name_User=='A431'   ~ 'CELL_LINE', 
                                             Sample_Name_User=='K562'   ~ 'CELL_LINE', 
                                             Sample_Name_User=='SPLEEN' ~ 'TISSUE', 
                                             TRUE ~ Sample_Type_User),
                  Sample_Input=as.integer(Sample_Input),
                  # Decoder_Pool_User=case_when(is.na(Decoder_Pool_User) ~ 'BETA', TRUE ~ Decoder_Pool_User),
                  Decoder_Pool_User=case_when(is.na(Decoder_Pool_User) ~ 'BETA', TRUE ~ Decoder_Pool_User)
    )
  if (verbose>=vt+4) print(ret_tib)
  
  ret_tib <- ret_tib %>% dplyr::filter(!is.na(Sentrix_Name)) %>%
    dplyr::select(Sentrix_Name, Sentrix_Barcode, Sentrix_Poscode, Bead_Pool, everything())
  if (verbose>=vt+4) print(ret_tib)
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]: Done."),"\n\n", sep='')
  
  ret_tib
}

getCurratedHumanSampleSheet = function(opt, file=NULL, verbose=0, vt=4) {
  funcTag <- 'getCurratedHumanSampleSheet'
  
  ss.src.list <- list.files(opt$ssSrcDir, full.names=TRUE)
  
  if (is.null(file)) 
    ss.tar.csv <- file.path(opt$sesDir, 'sampleSheets/official/Target_SampleSheet.csv')
  # else
  #   ss.tar.csv <- '/Users/bbarnes/Documents/Projects/workhorse/sampleSheets/EPIC-B4.SampleSheets.n4331.Dec1-2019.csv.gz'
  
  ss.tar.tib <- suppressMessages(suppressWarnings(readr::read_csv(ss.tar.csv) )) %>%
    dplyr::arrange(Sentrix_Name) %>%
    dplyr::mutate(Sample_Name=case_when(
      Sample_Name=='HEMI_METHYL' ~ 'H', 
      Sample_Name=='RA1' ~ 'RAJI',
      TRUE ~ Sample_Name) )
  
  ss.tar.tib %>% base::nrow()
  ss.tar.tib$Sentrix_Name %>% unique() %>% length()
  
  ss.all.tib <- NULL
  for (ii in seq(1:length(ss.src.list))) {
    cur.ss.csv <- ss.src.list[ii]
    cur.ss.tib <- suppressMessages(suppressWarnings(readr::read_csv(cur.ss.csv) ))
    
    ncols <- base::ncol(cur.ss.tib)
    nrows <- base::nrow(cur.ss.tib)
    
    # Skip Non-Standard Formats::
    if (cur.ss.csv==ss.tar.csv) next
    if (ncols!=10) next
    
    if (verbose>=vt) cat(glue::glue("[{funcTag}]: ncols={ncols}, nrows={nrows} {cur.ss.csv}"),"\n", sep='')
    # print(cur.ss.tib %>% head(n=2))
    ss.all.tib <- ss.all.tib %>% dplyr::bind_rows(cur.ss.tib)
  }
  
  ss.all.tib <- ss.all.tib %>% 
    dplyr::arrange(Sentrix_Name) %>%
    distinct() %>%
    dplyr::mutate(Sample_Name=case_when(
      Sample_Name=='HEMI_METHYL' ~ 'H', 
      Sample_Name=='RA1' ~ 'RAJI',
      TRUE ~ Sample_Name) )
  ss.all.tib %>% base::nrow()
  ss.all.tib$Sentrix_Name %>% unique() %>% length()
  
  # Resolve Overlaps and Check inconsistent Chip_Formats::
  ss.all.tib %>% dplyr::full_join(dplyr::select(ss.tar.tib,-c(Sentrix_ID,Sentrix_Position)), by="Sentrix_Name") %>% 
    dplyr::filter(!is.na(Chip_Format.x) | !is.na(Chip_Format.y)) %>%
    dplyr::filter(Chip_Format.x == Chip_Format.y)
  
  # Chip_Format Sample_Name Decoder_Pool Sample_Type Sample_Prep Sample_Plate Bead_Pool   Sentrix_ID Sentrix_Position
  ss.clean.tib <- ss.all.tib %>% dplyr::full_join(dplyr::select(ss.tar.tib,-c(Sentrix_ID,Sentrix_Position)), by="Sentrix_Name") %>%
    matchAndMergeFields(field='Sample_Plate', verbose=verbose,vt=vt) %>%
    matchAndMergeFields(field='Sample_Prep',  verbose=verbose,vt=vt) %>%
    matchAndMergeFields(field='Sample_Type',  verbose=verbose,vt=vt) %>%
    matchAndMergeFields(field='Sample_Name',  verbose=verbose,vt=vt) %>%
    matchAndMergeFields(field='Chip_Format',  verbose=verbose,vt=vt) %>%
    matchAndMergeFields(field='Bead_Pool',    verbose=verbose,vt=vt) %>%
    matchAndMergeFields(field='Decoder_Pool', verbose=verbose,vt=vt) %>%
    dplyr::mutate(Sample_Workflow=NA) %>%
    dplyr::arrange(Sentrix_Name) %>% 
    dplyr::group_by(Sentrix_ID) %>% 
    dplyr::mutate(Chip_Size=n(), Chip_Size=paste0(Chip_Size,'x1')) %>% 
    dplyr::mutate(Chip_Format=case_when(
      is.na(Chip_Format) ~ Chip_Size,
      TRUE ~ Chip_Format
    )) %>%
    dplyr::mutate(Bead_Pool=case_when(
      Bead_Pool=='EPIC_B4' ~ 'EPIC-B4',
      TRUE ~ Bead_Pool
    )) %>%
    dplyr::select(-Chip_Size) %>%
    dplyr::ungroup()
  
  ss.clean.ncols <- base::ncol(ss.clean.tib)
  ss.clean.nrows <- base::nrow(ss.clean.tib)
  colnames(ss.clean.tib)[2:ss.clean.ncols] <- paste(colnames(ss.clean.tib)[2:ss.clean.ncols],'User', sep='_')
  
  ss.clean.tib <- ss.clean.tib %>%
    dplyr::select(Sentrix_Name,Decoder_Pool_User,Bead_Pool_User,Chip_Format_User,Sample_Name_User,
                  Sample_Type_User,Sample_Prep_User,Sample_Plate_User,Sample_Workflow_User,Sentrix_ID_User,Sentrix_Position_User)
  
  ss.clean.tib
}


matchAndMergeFields = function(tib, field, verbose=0, vt=1) {
  funcTag <- 'matchAndMergeFields'
  
  strX <- paste0(field,'.x')
  strY <- paste0(field,'.y')
  
  field_x <- dplyr::sym(strX)
  field_y <- dplyr::sym(strY)
  # field_x <- base::quote(strX)
  # field_y <- base::quote(strY)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{TAB} field_x={field_x}, field_y={field_y}"),"\n", sep='')
  
  tib <- tib %>% 
    dplyr::mutate(!!field :=case_when(
      is.na(!!field_x) & is.na(!!field_y) ~ NA_character_,
      is.na(!!field_y) ~ !!field_x,
      is.na(!!field_x) ~ !!field_y,
      !!field_x==!!field_y ~ !!field_x,
      TRUE ~ 'failed'
    )) %>%
    dplyr::select(Sentrix_Name, field, !!field_x, !!field_y, everything()) %>%
    dplyr::select(-c(!!field_x, !!field_y))
  
  sum <- tib %>% dplyr::group_by(!!sym(field)) %>% dplyr::summarise(Counts=n())
  if (verbose>=vt) print(sum)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{TAB} Done."),"\n", sep='')
  
  tib
}

joinTibs = function(tib, org=NULL, cleanPrefix=FALSE, verbose=0, vt=1) {
  funcTag <- 'joinTibs'
  
  if (cleanPrefix)
    tib <- tib %>% purrr::set_names(stringr::str_remove(names(tib), '^[^.]+\\.' ) )
  
  if (is.null(org)) {
    jName <- names(tib)[1]
    tib <- tib %>% purrr::set_names(paste(names(tib),1, sep='_') ) %>%
      dplyr::rename(!!jName := 1)
    return(tib)
  }
  
  ncolsN <- base::ncol(tib)-1
  nrowsN <- base::nrow(tib)
  ncolsO <- base::ncol(org)-1
  nrowsO <- base::nrow(org)
  ncols  <- as.integer(ncolsO/ncolsN)
  
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{TAB}{TAB} new={nrowsN}x{ncolsN}, old={nrowsO}x{ncolsO},ncols={ncols}"),"\n", sep='')
  
  jName <- names(tib)[1]
  tib <- tib %>% purrr::set_names(paste(names(tib),ncols+1, sep='_') ) %>%
    dplyr::rename(!!jName := 1)
  
  tib <- org %>% dplyr::left_join(tib, by="Probe_ID")
  
  tib
}

loadBetaTib = function(file, minPval=0.05, verbose=0, vt=1) {
  funcTag <- 'loadBetaTib'
  
  tib <- suppressMessages(suppressWarnings(readr::read_rds(file)) ) %>%
    dplyr::bind_rows()
  
  mask_idx <- which(tib$NegsDetP>minPval)
  tib$Beta[mask_idx] <- NA
  
  tib <- tib %>% dplyr::select(Probe_ID, Beta)
  
  tib
}

loadProbeTib = function(file, minPval=0.05, verbose=0, vt=1) {
  funcTag <- 'loadBetaTib'
  
  tib <- suppressMessages(suppressWarnings(readr::read_rds(file)) ) %>%
    dplyr::bind_rows()
  
  mask_idx <- which(tib$NegsDetP>minPval)
  tib$Beta[mask_idx] <- NA
  # tib <- tib %>% dplyr::select(Probe_ID, Beta)
  
  tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                Plotting Functions for Pairwise Comparison::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

reduceExpNames = function(expA, expB, sep='_', verbose=0, vt=1) {
  funcTag <- 'reduceExpNames'
  
  expA_vec <- stringr::str_split(expA, sep, simplify=TRUE)
  expB_vec <- stringr::str_split(expB, sep, simplify=TRUE)
  
  exp_strA <- ''
  exp_strB <- ''
  min_len <- min(length(expA_vec), length(expB_vec))
  for (ii in c(1:min_len)) {
    if (verbose>=vt) cat(glue::glue("[{funcTag}]: ii={ii}, A={expA_vec[ii]}, B={expB_vec[ii]}"),"\n", sep='')
    if (expA_vec[ii]!=expB_vec[ii]) {
      exp_strA <- paste(exp_strA,expA_vec[ii], sep=sep)
      exp_strB <- paste(exp_strB,expB_vec[ii], sep=sep)
    }
  }
  exp_strA <- exp_strA %>% stringr::str_remove(paste0('^',sep))
  exp_strB <- exp_strB %>% stringr::str_remove(paste0('^',sep))
  
  ret_tib <- tibble::tibble(expA=exp_strA,
                            expB=exp_strB)
  
  ret_tib
}

plotTibPairs = function(tib, outType='auto', sample, expA, expB, probeType, outDir,
                        maxCnt=20000, wsize=2.8, minPval=0.02, minDelta=0.2,
                        colName, outFormat='pdf',verbose=0, vt=1) {
  funcTag <- 'plotTibPairs'
  
  if (!dir.exists(outDir)) dir.create(outDir)
  
  ncols <- base::ncol(tib)
  tarColIdxes <- c(4:ncols)
  
  nrows_org <- base::nrow(tib)
  p_tib <- tib
  if (nrows_org>maxCnt) p_tib <- tib %>% dplyr::sample_n(maxCnt)
  nrows_sub  <- base::nrow(p_tib)
  nrows_per  <- round(100*nrows_sub/nrows_org,1)
  nrows_orgK <- number_as_commaK(nrows_org)
  nrows_subK <- number_as_commaK(nrows_sub)
  
  if (nrows_sub<=10000) lowerSZ <- 0.5
  if (nrows_sub<=5000) lowerSZ <- 0.8
  if (nrows_sub<=500) lowerSZ <- 1.5
  
  upperAP <- 1.0
  lowerAP <- 1.0
  upperSZ <- 1.0
  lowerSZ <- 0.4
  if (outType=='auto') {
    if (nrows_sub<=10000) {
      plotType <- 'points'
      lowerPL <- geomPoints_Range
    } else {
      plotType <- 'density'
      lowerPL <- geomDensity2d_Range
    }
  } else if (outType=='points') {
    plotType <- outType
    lowerPL <- geomPoints_Range
  } else if  (outType=='density') {
    plotType <- outType
    lowerPL <- geomDensity2d_Range
  } else {
    stop("\n\n",glue::glue("[{funcTag}]: Unsupported OutType={outType}"),"\n\n", sep='')
  }
  cat(glue::glue("[{funcTag}]: nrows_sub={nrows_sub}"),"\n", sep='')
  
  plotName <- paste(sample,expA,'VS',expB,probeType, sep='_')
  
  gg_type  <- paste0('RSquared',colName)
  name_pdf <- glue::glue("{plotName}.{plotType}.{outFormat}")
  gg_pdf   <- file.path(outDir, name_pdf)
  
  gg_mtitle <- glue::glue("{expA} VS {expB} sample={sample}")
  gg_stitle <- glue::glue("ProbeType={probeType}, minPval={minPval}, {plotType}")
  gg_ctitle <- glue::glue("Plot displays downsampled percent = {nrows_per}% ",
                          "({nrows_subK}/{nrows_org})")
  
  gg <- GGally::ggpairs(p_tib, columns=tarColIdxes,
                        mapping = ggplot2::aes(color=Design_Type),
                        diag  = list(continuous=my_dens),
                        upper = list(combo = wrap("box_no_facet", alignPercent=0.7),
                                     continuous = wrap('blank', alpha=upperAP, size=upperSZ, wsize=wsize)),
                        lower = list(combo = "box_no_facet",
                                     continuous = wrap(lowerPL, alpha=lowerAP, size=lowerSZ, wsize=wsize) )
  )
  
  markTop <- FALSE
  for (idx_x in seq(tarColIdxes)) {
    tar.x <- tarColIdxes[idx_x]
    for (idx_y in seq(tarColIdxes)) {
      tar.y <- tarColIdxes[idx_y]
      
      if (idx_y<idx_x) {
        if (verbose>vt+3) cat(glue::glue("[{funcTag}]:{TAB}{TAB} Upper Corner Plot ",
                                         "idx_x={idx_x}, idx_y={idx_y}"),"\n", sep='')
        
        gg_yx <- deltaHist(p_tib,
                           idx.x=idx_x, idx.y=idx_y, tars=NULL, dets=NULL,
                           minPval=minPval, minDelta=minDelta, wsize=wsize,
                           maxCnt=maxCnt, verbose=verbose-3)
        if (verbose>vt+3) cat(glue::glue("[{funcTag}]:{TAB}{TAB} Adding ",
                                         "idx_x={idx_x}, idx_y={idx_y}"),"\n", sep='')
        gg <- GGally::putPlot(gg, gg_yx, idx_y,idx_x)
        if (verbose>vt+3) cat(glue::glue("[{funcTag}]:{TAB}{TAB} Done ",
                                         "idx_x={idx_x}, idx_y={idx_y}"),"\n", sep='')
        
      }
      if (markTop) break
    }
  }
  gg <- gg + 
    theme(panel.grid.major = element_blank()) +
    labs(title=gg_mtitle, subtitle=gg_stitle, caption=gg_ctitle)
  
  suppressMessages(ggplot2::ggsave(gg_pdf, gg) ) # , dpi=8)
  if (verbose>vt+2) cat(glue::glue("[{funcTag}]:{TAB} Plotted {gg_type} (pt={probeType}) PDF={gg_pdf}"),"\n\n", sep='')
  
  gg
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                        General Plotting Functions::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #


deltaHist = function(dat, idx.x, idx.y, tars=NULL, dets=NULL, minPval, minDelta, 
                     maxCnt=10000, wsize=2.8, verbose=0, vt=3) {
  funcTag <- 'deltaHist'
  
  if (verbose>vt) cat("\n\n",glue::glue("Starting {funcTag}"),"\n", sep='')
  
  # pdat <- dat %>% dplyr::select(1:3,tars[idx.x],dets[idx.x],tars[idx.y],dets[idx.y]) %>%
  pdat <- dat %>%
    dplyr::select(-Probe_ID) %>%
    dplyr::mutate(Delta=.[[3]]-.[[5]], Filter_Group=NA)
  
  pdat <- pdat %>% dplyr::select(Probe_Type,Design_Type,Filter_Group,Delta)
  pp_nrows  <- pdat %>% dplyr::filter(base::abs(Delta) <= minDelta) %>% base::nrow()
  tot.nrows <- pdat %>% base::nrow()
  pp.str <- paste0("\n",'Delta=',round(100*pp_nrows/tot.nrows, 2),'% ',"\n",
                   number_as_commaK(pp_nrows),'/',number_as_commaK(tot.nrows),' ')
  
  sumDat <- pdat %>% dplyr::group_by(Design_Type) %>% 
    dplyr::summarise(PassCount=count(abs(Delta)<=minDelta, na.rm=TRUE),
                     FailCount=count(abs(Delta)>minDelta | is.na(minDelta)),
                     Total=n(), PassPer=round(100*PassCount/Total,2)) %>%
    split(.$Design_Type)
  # print(sumDat)
  
  pp.str <- "\n"
  for (sp in names(sumDat)) {
    # pCnt <- sumDat[[sp]] %>% head(n=1) %>% .$PassCount
    # nCnt <- sumDat[[sp]] %>% 
    pCnt <- sumDat[[sp]] %>% head(n=1) %>% .$FailCount
    dVal <- sumDat[[sp]] %>% head(n=1) %>% .$PassPer
    if (is.na(pCnt)) pCnt <- 0
    if (is.na(dVal)) dVal <- 0
    
    pp.str <- paste0(pp.str,sp,'=',round(dVal,1),'% (-',number_as_commaK(pCnt),')\n')
  }
  if (verbose>vt) cat(glue::glue("pp.str={pp.str}"),"\n", sep='')
  # return(pdat)
  
  # SubSample for Plotting::
  #
  # TDB:: UnGroup => SumSample => Regroup???
  #
  if (verbose>vt) cat(glue::glue("[{funcTag}]: Pre-SubSize(pdat)={tot.nrows}"),"\n", sep='')
  if (tot.nrows>maxCnt) pdat <- pdat %>% dplyr::sample_n(maxCnt)
  tot.nrows <- pdat %>% base::nrow()
  if (verbose>vt) cat(glue::glue("[{funcTag}]: Pos-SubSize(pdat)={tot.nrows}"),"\n", sep='')
  
  xmin <- -1
  xmid <- 0
  xmax <- 1
  
  colorFilter <- FALSE
  colorFilter <- TRUE
  gg <- NULL
  if (colorFilter)
    gg <- ggplot2::ggplot(data=pdat, aes(x=Delta, color=Design_Type))
  else 
    gg <- ggplot2::ggplot(data=pdat, aes(x=Delta))
  
  gg <- gg +
    geom_histogram(binwidth=0.01, alpha=0.1, na.rm=TRUE) +
    geom_vline(xintercept = -0.2, color="red",  linetype="dotted") +
    geom_vline(xintercept =  0.2, color="blue", linetype="dotted") +
    scale_x_continuous(breaks = c(xmin, xmid, xmax)) +
    # annotate("text",  x=Inf, y = Inf, label = pp.str,
    geom_text(x=Inf, y=Inf, label=pp.str,
              vjust=1, hjust=1,
              # nudge_x = -0.1, nudge_y = -0.1,
              size=wsize, # fontface = "bold",
              family='mono',
              alpha=0.3,
              inherit.aes = FALSE)
  
  if (verbose>vt) cat(glue::glue("[{funcTag}]: Done."),"\n\n", sep='')
  
  gg
}
geomDensity2d_Range <- function(data,mapping, field='Design_Type', alpha,size, 
                                xmin=0,xmax=1,ymin=0,ymax=1, alignPercent=0.8, wsize=3, verbose=0) {
  funcTag <- 'geomDensity2d_Range'
  vt <- 4
  
  xmid=round((xmax-xmin)/2, 1)
  ymid=round((ymax-ymin)/2, 1)
  
  # Get Correlation Coeffiecent String
  cor_str <- ''
  data_split <- data %>% split(data[[field]])
  for (sp in names(data_split)) {
    x <- GGally::eval_data_col(data_split[[sp]], mapping$x)
    y <- GGally::eval_data_col(data_split[[sp]], mapping$y)
    cor <- cor(x, y,  method='pearson', use='pairwise.complete.obs')
    cor_str <- paste0(cor_str,sp,'=',round(cor,3),'\n')
  }
  cor_str <- cor_str %>% stringr::str_remove('\\n$')
  
  ggplot(data=data, mapping=mapping)+
    geom_density2d(alpha=alpha, size=size)+
    scale_x_continuous(breaks = c(xmin, xmid, xmax))+
    scale_y_continuous(breaks = c(ymin, ymid, ymax))+
    ggplot2::geom_label(
      data = data.frame(
        xlabel = min(x, na.rm = TRUE),
        ylabel = max(y, na.rm = TRUE),
        lab = cor_str),
      mapping = ggplot2::aes(x = xlabel,
                             y = ylabel, 
                             label = lab),
      hjust = 0, vjust = 1,
      size = wsize, fontface = "bold",
      alpha=0.3, family='mono',
      inherit.aes = FALSE # do not inherit anything from the ...
    )
}

geomPoints_Range <- function(data,mapping, field='Design_Type', method='pearson', alpha,size, 
                             xmin=0,xmax=1,ymin=0,ymax=1, alignPercent=0.8, wsize=3, verbose=0) {
  funcTag <- 'geomPoints_Range'
  vt <- 4
  
  xmid=round((xmax-xmin)/2, 1)
  ymid=round((ymax-ymin)/2, 1)
  
  # Get Correlation Coeffiecent String
  cor_str <- ''
  data_split <- data %>% split(data[[field]])
  for (sp in names(data_split)) {
    x <- GGally::eval_data_col(data_split[[sp]], mapping$x)
    y <- GGally::eval_data_col(data_split[[sp]], mapping$y)
    cor <- cor(x, y,  method='pearson', use='pairwise.complete.obs')
    cor_str <- paste0(cor_str,sp,'=',round(cor,3),'\n')
  }
  cor_str <- cor_str %>% stringr::str_remove('\\n$')
  
  ggplot(data=data, mapping=mapping)+
    geom_point(alpha=alpha, size=size)+
    scale_x_continuous(breaks = c(xmin, xmid, xmax))+
    scale_y_continuous(breaks = c(ymin, ymid, ymax))+
    ggplot2::geom_label(
      data = data.frame(
        xlabel = min(x, na.rm = TRUE),
        ylabel = max(y, na.rm = TRUE),
        lab = cor_str),
      mapping = ggplot2::aes(x = xlabel,
                             y = ylabel, 
                             label = lab),
      hjust = 0, vjust = 1,
      size = wsize, fontface = "bold",
      alpha=0.3, family='mono',
      inherit.aes = FALSE # do not inherit anything from the ...
    )
}

my_dens <- function(data, mapping, alpha=0.8, ..., low = "#132B43", high = "#56B1F7",
                    xmin=0,xmax=1,ymin=0,ymax=1) {
  xmid=round((xmax-xmin)/2, 1)
  ymid=round((ymax-ymin)/2, 1)
  
  ggplot(data = data, mapping=mapping) +
    geom_density(..., alpha=alpha) +
    scale_x_continuous(breaks = c(xmin, xmid, xmax))
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                    General Text Formatting Functions::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

number_as_commaK <- function(value) {
  if (!is.numeric(value)) {
    cat("[ERROR]: Not-Numeric=",value,"\n\n", sep='')
  }
  stopifnot(is.numeric(value))
  
  s <- ''
  if (value>10000 || value==10000) {
    s <- 'k'
    value <- as.integer(value/1000)
  }
  paste0(format(value, big.mark = ",", scientific = FALSE, digits = 22),s)
}

# End of file
