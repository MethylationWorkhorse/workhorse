
rm(list=ls(all=TRUE))
while (dev.cur()>1) dev.off()

suppressWarnings(suppressPackageStartupMessages(require("optparse",quietly=TRUE)))

suppressWarnings(suppressPackageStartupMessages(require("tidyverse")) )
suppressWarnings(suppressPackageStartupMessages(require("stringr")) )
suppressWarnings(suppressPackageStartupMessages(require("glue")) )
suppressWarnings(suppressPackageStartupMessages(require("scales")) )
suppressWarnings(suppressPackageStartupMessages(require("matrixStats")) )

COM <- ","
TAB <- "\t"
RET <- "\n"

opt <- NULL
opt$prgmTag <- 'buildMedianReferenceSamples'
opt$runMode  <- NULL

opt$Rscript <- 'Rscript'
opt$topDir  <- '/Users/bbarnes/Documents/Projects/workhorse'
opt$srcDir  <- '/Users/bbarnes/Documents/Projects/workhorse/scripts'
opt$outDir  <- NULL
opt$idatsDir <- NULL
opt$funcScript <- 'workhorse_functions.R'
opt$sampleSheet <- NULL

opt$platform  <- 'EPIC-B4'
opt$build     <- 'hg19'
opt$method    <- 'both'

opt$poobMinPval <- 0.2
opt$negsMinPval <- 0.02

opt$decoderPool <- 'BETA'
# opt$decoderPool <- 'DELTA'



opt$verbosity <- 3

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                          Sample Sheet Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
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

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                    Build Human Reference Sample Sheet::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
opt$sesDir   <- '/Users/bbarnes/Documents/Projects/darkmatter/Projects/sesamize'
opt$ssSrcDir <- file.path(opt$sesDir, 'sampleSheets/official')

ss.src.list <- list.files(opt$ssSrcDir, full.names=TRUE)

ss.inp.csv <- file.path(opt$sesDir, 'sampleSheets/official/nano_gram_input.SampleSheet.tsv')
ss.inp.tib <- suppressMessages(suppressWarnings(readr::read_tsv(ss.inp.csv) )) %>%
  dplyr::arrange(Sentrix_Name)
ss.inp.tib %>% base::nrow()
ss.inp.tib$Sentrix_Name %>% unique() %>% length()

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

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                    Load Tibs To Build Grouped Results
#                         From Reference Tibs (Beta)
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# opt$refList <- list.files(file.path(opt$sesDir, 'idats_refs'), 
#                           pattern='EPIC-B0_both.SampleSheet.csv.gz', recursive=TRUE, full.names=TRUE)
# 
# opt$refList <- list.files(file.path(opt$sesDir, 'idats_refs'), 
#                           pattern='EPIC-B2.sset_core.Auto_SampleSheet.csv.gz', recursive=TRUE, full.names=TRUE)

opt$refList <- list.files(file.path(opt$topDir, 'builds/workhorse/References'), 
                          pattern='EPIC-B4_both.SampleSheet.csv.gz', recursive=TRUE, full.names=TRUE)

ss.ref.tib <- suppressMessages(suppressWarnings(lapply(opt$refList, readr::read_csv))) %>%
  dplyr::bind_rows() %>%
  dplyr::mutate(Auto_Sample_DB_Key=stringr::str_remove(Auto_Sample_DB_Key, '_Beta$'),
                Auto_Sample_R2_Key=stringr::str_remove(Auto_Sample_R2_Key, '_Beta$')) %>%
  dplyr::arrange(Sentrix_Name) %>%
  split(.$Auto_Sample_R2_Key)
  # dplyr::group_split()

# man.ref.tib <- foreach (sample=names(ss.ref.tib), .inorder=T, .final = function(x) setNames(x, names(ss.ref.tib))) %dopar% {
sam_tib <- NULL
for (sample in names(ss.ref.tib)) {
  # sample <- 'U'
  ss.cur.tib <- ss.clean.tib %>% dplyr::inner_join(ss.ref.tib[[sample]], by="Sentrix_Name")
  # print(ss.cur.tib)
  
  ss.tib.locs <- file.path(opt$topDir, 'builds/workhorse/References', ss.cur.tib$Barcode, opt$platform)
  ss.tib.list <- paste0(ss.tib.locs,'/',ss.cur.tib$Sentrix_Name, '_EPIC-B4_both.Probes.tib.rds')
  # cat(ss.tib.locs,"\t",ss.tib.list, sep=' ')
  tibs <- lapply(ss.tib.list, loadBetaTib) %>%
    dplyr::bind_cols() %>% 
    dplyr::select(Probe_ID, starts_with('Beta'))
  
  sam_tib[[sample]] <- tibs
  # sam_tib <- sam_tib %>% dplyr::bind_rows(.id=sample)
}









ss.clean.tib %>% dplyr::inner_join(sam_tib[[sample]], by="Probe_ID")

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                           Build Reference Tibs::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

opt$autoList <- list.files(file.path(opt$topDir, 'sampleSheets'), 
                          pattern='EPIC-B4_both.SampleSheet.csv.gz', recursive=TRUE, full.names=TRUE)

ss.auto.rds <- file.path(opt$topDir, 'sampleSheets/rds/auto_EPIC-B4_both.SampleSheet.csv.gz')
ss.auto.tib <- suppressMessages(suppressWarnings(lapply(opt$autoList, readr::read_csv))) %>%
  dplyr::bind_rows() %>%
  dplyr::mutate(Auto_Sample_DB_Key=stringr::str_remove(Auto_Sample_DB_Key, '_Beta$'),
                Auto_Sample_R2_Key=stringr::str_remove(Auto_Sample_R2_Key, '_Beta$')) %>%
  dplyr::arrange(Sentrix_Name)
readr::write_rds(ss.auto.tib, ss.auto.rds, compress='gz')

ss.auto.tib <- ss.auto.tib %>% dplyr::arrange(Sentrix_Name)

# Current:: U=202761400009_R05C01
# Current:: H=202761400009_R06C01
# Current:: M=202761400009_R07C01
ss.clean.tib %>% dplyr::inner_join(dplyr::arrange(ss.auto.tib, Sentrix_Name), by="Sentrix_Name") %>%
  dplyr::select(Sentrix_Name, Decoder_Pool_User, Sample_Name_User,SwappedCount_CG,Bead_Pool_User,
                cg_Beta_mean,cg_Beta_median,cg_Beta_median1,cg_G_median,cg_R_median,
                PassDetpNegs_Percent_CG,PassDetpPoob_Percent_CG,
                Auto_Sample_DB_Key,Auto_Sample_R2_Key,Auto_Sample_DB_Val,Auto_Sample_R2_Val,
                everything()) %>%
  dplyr::filter(Sample_Name_User=='U') %>%
  dplyr::arrange(cg_Beta_mean) %>% dplyr::filter(Decoder_Pool_User=='BETA', Bead_Pool_User=='EPIC-B4') %>%
  select(1:17) %>% as.data.frame() 

  # dplyr::filter(Sentrix_Name=='202761400009_R05C01')


# Add Chunk Index
minLen <- 100
grpIdx <- 0
man.chr.tib <- man.cpg.tib %>% dplyr::arrange(seqnames, start) %>%
  dplyr::mutate(Local_Group_Key=0, Local_Group_Cnt=0) %>% split(.$seqnames)
for (chr in names(man.chr.tib)) {
  nrows <- base::nrow(man.chr.tib[[chr]])

  prePos <- 0
  grpCnt <- 0
  stime <- system.time({
    for (idx in seq(1:nrows)) {
      curPos <- man.chr.tib[[chr]]$start[idx]
  
      if (curPos-prePos < minLen) {
        grpCnt <- grpCnt + 1
      } else {
        grpCnt <- 1
        grpIdx <- grpIdx + 1
        prePos <- curPos
      }
      man.chr.tib[[chr]]$Local_Group_Key[idx] <- grpIdx
      man.chr.tib[[chr]]$Local_Group_Cnt[idx] <- grpCnt
      
      # if (grpIdx>10) break
    }
  })
  cat("\t",glue::glue("[{funcTag}]: grpIdx={grpIdx}, prePos={prePos}"),"\n", sep='')
  print(stime)
  break
}

ses.man.tib <- sesameData::sesameDataGet('EPIC.hg19.manifest')

man.chr.tib %>% dplyr::filter(Local_Group_Cnt>5) %>% dplyr::arrange(-Local_Group_Cnt) %>% dplyr::select(c(1,3,4,10:19)) %>% dplyr::select(Local_Group_Key) %>% distinct()

sam.prb.rds <- file.path(opt$topDir, 'builds/workhorse/DeltaBetaCore/202761400009/EPIC-B4/202761400009_R01C01_EPIC-B4_both.Probes.tib.rds')
sam.prb.rds <- file.path(opt$topDir, 'workspace/workhorse/202761400009/EPIC-B4/202761400009_R01C01_EPIC-B4_both.Probes.tib.rds')
sam.prb.tib <- readr::read_rds(sam.prb.rds)

# Probe_ID   Probe_Type  Beta   Mval NegsDetP PoobDetP NegsDetP_Delta PoobDetP_Delta     G     R
# Probe_ID   Probe_Type Design_Type Swapped   Beta   Mval NegsDetP PoobDetP NegsDetP_Delta PoobDetP_Delta    GM    GU    RM    RU
bind.prb.tib <- sam.prb.tib$II %>% dplyr::mutate(Design_Type='II', Swapped=NA) %>%
  dplyr::bind_rows(sam.prb.tib$I) %>% dplyr::arrange(Probe_ID)

bind.prb.csv <- file.path(opt$topDir, 'workspace/202761400009_R01C01_EPIC-B4_both.Probes.tib.csv.gz')
bind.prb.rds <- file.path(opt$topDir, 'workspace/202761400009_R01C01_EPIC-B4_both.Probes.tib.rds')
readr::write_csv(bind.prb.tib, bind.prb.csv)
readr::write_rds(bind.prb.tib, bind.prb.rds, compress="gz")

man.grp.sorted.tib <- man.grp.tib %>% dplyr::arrange(Probe_ID)

full.grp.tib <- man.grp.sorted.tib %>% dplyr::inner_join(bind.prb.tib, by="Probe_ID") 

# full.grp.tib %>% dplyr::filter(Local_Group_Cnt>3) %>% dplyr::group_by(Local_Group_Key) %>% select(Beta) %>% summary()

# STATS ON GROUPS::
full.grp.tib %>% dplyr::filter(Local_Group_Cnt>3) %>% 
  dplyr::group_by(seqnames, Local_Group_Key,Local_Group_Cnt) %>% 
  summarise(Beta_Avg=mean(Beta, na.rm=TRUE), Beta_Med=median(Beta, na.rm=TRUE), Beta_SD=sd(Beta, na.rm=TRUE))



opt$loadRDS <- TRUE
opt$retPRBS <- TRUE
opt$retSSET <- TRUE
opt$verbosity <- 10

ssets <- singleSampleWorkflow(prefix=chipPrefixes[[prefix]],
                     opt=opt,
                     man.cpg = man.cpg.tib,
                     man.add = man.add.tib,
                     man.ses = NULL,
                     ref.tib = ref.sam.tib,
                     ctl.tib = NULL,
                     vt=1)

prbs2 <- ssetsToTibs(ssets, verbose=opt$verbosity)


prbs <- singleSampleWorkflow(prefix=chipPrefixes[[prefix]],
                             opt=opt,
                             man.cpg = man.cpg.tib,
                             man.add = man.add.tib,
                             man.ses = NULL,
                             ref.tib = ref.sam.tib,
                             ctl.tib = NULL,
                             vt=1)


addT = function(stime, ptimes=NULL, name=NULL, printTime=TRUE) {
  
  order <- 1
  if (!is.null(ptimes)) order <- base::nrow(ptimes) + 1
  if (is.null(name)) name <- paste0('method-',order)
  
  curRow <- stime %>% tibble::enframe() %>% 
    tidyr::spread(name, value) %>% 
    dplyr::mutate(Method=name, Order=order) %>% 
    dplyr::select(Method,Order, everything())
  
  ptimes <- ptimes %>% dplyr::bind_rows(curRow)
  
  if (printTime==TRUE) print(ptimes)
  ptimes
}

ss.clean.tib %>% dplyr::inner_join(dplyr::arrange(ss.auto.tib, Sentrix_Name), by="Sentrix_Name") %>%
  dplyr::select(Sentrix_Name, Decoder_Pool_User, Sample_Name_User,SwappedCount_CG,Bead_Pool_User,
                cg_Beta_mean,PassDetpNegs_Percent_CG,PassDetpPoob_Percent_CG,
                Auto_Sample_DB_Key,Auto_Sample_R2_Key,Auto_Sample_DB_Val,Auto_Sample_R2_Val,
                everything()) %>%
  dplyr::filter(Sample_Name_User=='U') %>%
  dplyr::arrange(-cg_Beta_mean)


  dplyr::arrange(-PassDetpNegs_Percent_CG)
  dplyr::arrange(-Auto_Sample_R2_Val)
  dplyr::arrange(-PassDetpPoob_Percent_CG)
  dplyr::arrange(-Auto_Sample_DB_Val)



ss.auto.tib %>% dplyr::filter(Auto_Sample_R2_Key=='U') %>% 
  dplyr::select(Sentrix_Name, 
                Auto_Sample_R2_Key,Auto_Sample_DB_Key,Auto_Sample_DB_Val,Auto_Sample_R2_Val, 
                everything()) %>% 
  dplyr::filter(Auto_Sample_R2_Key==Auto_Sample_DB_Key)


ss.full.tib <- ss.clean.tib %>% dplyr::inner_join(ss.auto.tib, by="Sentrix_Name") %>%
  dplyr::select(Sentrix_Name,Sample_Name_User,Auto_Sample_DB_Key,Auto_Sample_R2_Key, everything())

ss.mis.db.tib <- ss.full.tib %>%
  dplyr::filter(Sample_Name_User!=Auto_Sample_DB_Key)

ss.mat.r2.tib <- ss.full.tib %>%
  dplyr::filter(Sample_Name_User==Auto_Sample_R2_Key)


# End of file
