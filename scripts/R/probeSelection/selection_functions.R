

suppressWarnings(suppressPackageStartupMessages(require("tidyverse")) )
suppressWarnings(suppressPackageStartupMessages(require("stringr")) )
suppressWarnings(suppressPackageStartupMessages(require("glue")) )
suppressWarnings(suppressPackageStartupMessages(require("scales")) )
suppressWarnings(suppressPackageStartupMessages(require("matrixStats")) )

# Parallel Computing Packages
suppressWarnings(suppressPackageStartupMessages(require("doParallel")) )

doParallel::registerDoParallel()
num_cores   <- detectCores()
num_workers <- getDoParWorkers()

COM <- ","
TAB <- "\t"
RET <- "\n"

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                          Sample Sheet Methods::
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

getCurratedHumanSampleSheet = function(opt, verbose=0, vt=1) {
  funcTag <- 'getCurratedHumanSampleSheet'
  
  ss.src.list <- list.files(opt$ssSrcDir, full.names=TRUE)
  
  # Nano-gram input is not added yet!!!
  #
  # ss.inp.csv <- file.path(opt$sesDir, 'sampleSheets/official/nano_gram_input.SampleSheet.tsv')
  # ss.inp.tib <- suppressMessages(suppressWarnings(readr::read_tsv(ss.inp.csv) )) %>%
  #   dplyr::arrange(Sentrix_Name)
  # ss.inp.tib %>% base::nrow()
  # ss.inp.tib$Sentrix_Name %>% unique() %>% length()
  
  
  ss.tar.csv <- file.path(opt$sesDir, 'sampleSheets/official/Target_SampleSheet.csv')
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
    
    cat(glue::glue("[SampleSheet]: ncols={ncols}, nrows={nrows} {cur.ss.csv}"),"\n", sep='')
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
    matchAndMergeFields(field='Sample_Plate') %>%
    matchAndMergeFields(field='Sample_Prep') %>%
    matchAndMergeFields(field='Sample_Type') %>%
    matchAndMergeFields(field='Sample_Name') %>%
    matchAndMergeFields(field='Chip_Format') %>%
    matchAndMergeFields(field='Bead_Pool') %>%
    matchAndMergeFields(field='Decoder_Pool') %>%
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
  
  ss.clean.tib
}


matchAndMergeFields = function(tib, field, summary=TRUE, verbose=0, vt=1) {
  funcTag <- 'matchAndMergeFields'
  
  strX <- paste0(field,'.x')
  strY <- paste0(field,'.y')
  
  field_x <- dplyr::sym(strX)
  field_y <- dplyr::sym(strY)
  # field_x <- base::quote(strX)
  # field_y <- base::quote(strY)
  cat(glue::glue("[{funcTag}]:{TAB} field_x={field_x}, field_y={field_y}"),"\n", sep='')
  cat("\n")
  
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
  
  if (summary) {
    sum <- tib %>% dplyr::group_by(!!sym(field)) %>% dplyr::summarise(Counts=n())
    print(sum)
  }
  
  tib
}

joinTibs = function(tib, org=NULL) {
  funcTag <- 'joinTibs'
  
  if (is.null(org)) return(tib)
  
  ncols <- base::ncol(org)
  cat("\t",glue::glue("[{funcTag}]: ncols={ncols}"),"\n", sep='')
  
  s1 <- paste0('_',ncols-1)
  s2 <- paste0('_',ncols)
  
  tib <- org %>% dplyr::left_join(tib, by="Probe_ID", suffix=c(s1,s2))
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

# End of file
