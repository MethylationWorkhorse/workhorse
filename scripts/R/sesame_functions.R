
# Safe Sesame Functions

suppressWarnings(suppressPackageStartupMessages(require("tidyverse")) )
suppressWarnings(suppressPackageStartupMessages(require("stringr")) )
suppressWarnings(suppressPackageStartupMessages(require("glue")) )
suppressWarnings(suppressPackageStartupMessages(require("scales")) )
suppressWarnings(suppressPackageStartupMessages(require("matrixStats")) )

COM <- ","
TAB <- "\t"
RET <- "\n"

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                      Sesame Sample Methods::Try-Catched
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

safeGCT = function(sset, verbose=0,vt=4,tc=1,tabsStr='') {
  funcTag <- 'safeGCT'
  for (i in c(1:tc)) tabsStr <- paste0(tabsStr, TAB)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting..."),"\n", sep='')
  
  GCT = tryCatch({
    suppressWarnings(sesame::bisConversionControl(sset) )
  }, warning = function(w) {
    'warning-GCT'
  }, error = function(e) {
    'error-GCT'
  }, finally = {
    'cleanup-GCT'
  })
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done."),"\n\n", sep='')  
  
  GCT
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                          Sesame Manifest Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

getManifest = function(man, dat, verbose=0,vt=4,tc=1,tabsStr='') {
  funcTag <- 'getManifest'
  for (i in c(1:tc)) tabsStr <- paste0(tabsStr, TAB)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting..."),"\n", sep='')
  if (verbose>=vt+5) print(man)
  
  tib <- man %>% dplyr::right_join(dat, by=c("U"="Address")) %>% 
    dplyr::filter(!is.na(Probe_ID)) %>%
    dplyr::select(Probe_ID, M, U, DESIGN, COLOR_CHANNEL, col, Probe_Type, Probe_Source, Next_Base) %>%
    dplyr::arrange(Probe_ID)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done."),"\n\n", sep='')
  
  tib
}

addManifestGroups = function(platform='EPIC',build='hg19', manCPG=NULL, maxCnt=NULL, 
                             verbose=0,vt=4,tc=1,tabsStr='') {
  funcTag <- 'addManifestGroups'
  for (i in c(1:tc)) tabsStr <- paste0(tabsStr, TAB)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting platform={platform}, build={build}..."),"\n", sep='')
  
  minLen <- 100
  grpIdx <- 0
  if (is.null(manCPG))
    manCPG <- sesameManifest(platform='EPIC', build=opt$build, verbose=opt$verbosity)
  
  man.chr.tib <- manCPG %>% dplyr::arrange(seqnames, start) %>%
    dplyr::filter(!str_starts(Probe_ID, 'ctl')) %>%
    dplyr::filter(str_starts(seqnames, 'chr')) %>%
    dplyr::mutate(Local_Group_Key=0, Local_Group_Idx=0) %>% 
    dplyr::mutate(Local_Group_Key=as.integer(Local_Group_Key), 
                  Local_Group_Idx=as.integer(Local_Group_Idx)) %>%
    split(.$seqnames)
  
  if (length(man.chr.tib[['*']])>0) man.chr.tib[['*']] <- NULL
  if (verbose>=vt+2) print(names(man.chr.tib))
  
  man.ret.tib <- NULL
  # for (chr in names(man.chr.tib)) {
  man.ret.tib <- foreach (chr=names(man.chr.tib), .inorder=T, .final = function(x) setNames(x, names(man.chr.tib))) %dopar% {
    man.tib <- man.chr.tib[[chr]]
    if (!is.null(maxCnt)) man.tib <- man.chr.tib[[chr]] %>% head(n=maxCnt)
    
    nrows <- base::nrow(man.tib)
    if (verbose>=vt)
      cat("\t",glue::glue("[{funcTag}]:{TAB} chr={chr}, grpIdx={grpIdx}, nrows={nrows}"),"\n", sep='')
    
    prePos <- 0
    grpCnt <- 0
    for (idx in seq(1:nrows)) {
      curPos <- man.tib$start[idx]
      
      if (curPos-prePos < minLen) {
        grpCnt <- grpCnt + 1
      } else {
        grpCnt <- 1
        grpIdx <- grpIdx + 1
        prePos <- curPos
      }
      man.tib$Local_Group_Key[idx] <- grpIdx
      man.tib$Local_Group_Idx[idx] <- grpCnt
      
      # if (grpIdx>10) break
    }
    man.tib <- man.tib %>% dplyr::add_count(Local_Group_Key, name="Local_Group_Cnt") %>%
      dplyr::mutate(Local_Group_Key=as.integer(Local_Group_Key), 
                    Local_Group_Idx=as.integer(Local_Group_Idx))
    
    if (verbose>=vt)
      cat("\t",glue::glue("[{funcTag}]:{TAB} grpIdx={grpIdx}, prePos={prePos}"),"\n", sep='')
    # break
    man.tib
  }
  man.ret.tib <- dplyr::bind_rows(man.ret.tib)
  if (verbose>=vt+2) print(man.ret.tib)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done."),"\n\n", sep='')
  
  man.ret.tib
}

sesameManifest = function(platform='EPIC', build='hg19', verbose=0, vt=4) {
  funcTag <- 'sesameManifest'
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{TAB} Starting..."),"\n", sep='')
  name <- paste(platform,build,'manifest', sep='.')
  # tib <- suppressMessages(suppressWarnings(sesameData::sesameDataGet(name) )) %>% 
  tib <- sesameData::sesameDataGet(name) %>% 
    as.data.frame() %>% 
    rownames_to_column('Probe_ID') %>% 
    tibble::as_tibble() %>% 
    dplyr::mutate(col=NA, Masked=NA) %>%
    dplyr::mutate(
      col=case_when(
        stringr::str_starts(channel, 'G') ~ 'G',
        stringr::str_starts(channel, 'R') ~ 'R'),
      Masked=case_when(
        MASK_mapping==TRUE ~ 'Mapping',
        MASK_general==TRUE ~ 'General'
      ) ) %>%
    dplyr::mutate(M=address_B,
                  U=address_A,
                  COLOR_CHANNEL=channel,
                  Probe_Type=probeType,
                  Probe_Source='EPIC-B4',
                  Next_Base=nextBaseRef,
                  Probe_CpG_Cnt=probeCpGcnt) %>%
    dplyr::rename(DESIGN=designType) %>%
    dplyr::select(Probe_ID, M, U, DESIGN, COLOR_CHANNEL, col, Probe_Type, Probe_Source, Next_Base, 
                  seqnames, start, end, width, strand, Probe_CpG_Cnt, Masked) %>%
    dplyr::arrange(Probe_ID)
  
  # return(tib)
  
  name <- paste(platform,'address', sep='.')
  ses.all.que <- sesameData::sesameDataGet(name)
  
  # Probe_ID        M        U DESIGN COLOR_CHANNEL  col
  # ses.org.tib <- ses.all.que$ordering %>% tibble::as_tibble()
  ses.ctl.tib <- ses.all.que$controls
  ses.ctl.tib <- ses.ctl.tib %>% 
    dplyr::rename(U=Address) %>%
    dplyr::mutate(U=as.character(U)) %>%
    dplyr::mutate(U=as.integer(U),
                  M=NA,
                  M=as.integer(M),
                  Probe_ID=str_replace_all(Name, "[-\\(\\) ]+",'_'),
                  Probe_ID=str_replace_all(Probe_ID, "_+$",''),
                  Probe_ID=paste('ctl',Probe_ID,  sep='_'),
                  M=NA, 
                  COLOR_CHANNEL=as.character(Color_Channel),
                  col=NA,
                  Probe_Type=str_replace_all(Type, "\\s", '_'),
                  Probe_Source='EPIC-B4',
                  Next_Base=NA) %>% 
    dplyr::mutate(col=as.character(col),
                  Next_Base=as.character(Next_Base)) %>%
    tibble::as_tibble() %>%
    dplyr::mutate(DESIGN='II') %>%
    dplyr::select(Probe_ID, M, U, DESIGN, COLOR_CHANNEL, col, Probe_Type, Probe_Source, Next_Base)
  
  ses.all.man <- ses.ctl.tib %>% dplyr::bind_rows(tib) %>%
    dplyr::arrange(Probe_ID)
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{TAB} Done."),"\n", sep='')
  
  ses.all.man
}



# End of file
