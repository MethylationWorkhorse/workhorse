
suppressWarnings(suppressPackageStartupMessages(require("tidyverse")) )
suppressWarnings(suppressPackageStartupMessages(require("stringr")) )
suppressWarnings(suppressPackageStartupMessages(require("glue")) )
suppressWarnings(suppressPackageStartupMessages(require("scales")) )
suppressWarnings(suppressPackageStartupMessages(require("matrixStats")) )

suppressWarnings(suppressPackageStartupMessages(require("Biostrings")) )


COM <- ","
TAB <- "\t"
RET <- "\n"

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                Infinium Methylation Probe toStrings::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

srdsToBrac = function(tib, 
                      beg1=1, end1=60, mid1=61,
                      beg2=63,end2=122,mid2=62) {
  # TBD:: Calculate alll data points based on sequence length
  
  tib %>% 
    dplyr::mutate(StrandFR=case_when(FR ~ 'F', !FR ~ 'R', TRUE ~ NA_character_),
                  StrandCO=case_when(CO ~ 'C', !CO ~ 'O', TRUE ~ NA_character_),
                  DesSeqN=paste0(stringr::str_sub(DesSeqN,beg1,end1),{BNG},
                                 stringr::str_sub(DesSeqN,mid1,mid2),{BNG},
                                 stringr::str_sub(DesSeqN,beg2,end2)),
                  
                  DesBscU=paste0(stringr::str_sub(DesBscU,beg1,end1),{BNG},
                                 stringr::str_sub(DesBscU,mid1,mid2),{BNG},
                                 stringr::str_sub(DesBscU,beg2,end2)),
                  
                  DesBscM=paste0(stringr::str_sub(DesBscM,beg1,end1),{BNG},
                                 stringr::str_sub(DesBscM,mid1,mid2),{BNG},
                                 stringr::str_sub(DesBscM,beg2,end2)),
                  
                  DesBscD=paste0(stringr::str_sub(DesBscD,beg1,end1),{BNG},
                                 stringr::str_sub(DesBscD,mid1,mid2),{BNG},
                                 stringr::str_sub(DesBscD,beg2,end2))) %>%
    dplyr::arrange(StrandFR, StrandCO) %>% split(.$StrandFR)
  
}

prbsToStr = function(tib, pr='cg',
                     verbose=0,vt=1,tc=1,tabsStr='') {
  funcTag <- 'prbsToStr'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  fr1Key <- tib$StrandFR[1]
  fr2Key <- tib$StrandFR[2]
  
  co1Key <- tib$StrandCO[1]
  co2Key <- tib$StrandCO[2]
  
  mud <- list('U'='U','M'='M','D'='D')

  desSeq <- 'DesSeqN'
  bscKey <- lapply(mud, function(x) { paste('DesBsc',x,sep='')} )
  nxbKey <- lapply(mud, function(x) { paste('NXB',x,sep='_')} )
  cpnKey <- lapply(mud, function(x) { paste('CPN',x,sep='_')} )
  tarKey <- lapply(mud, function(x) { paste('TAR',x,sep='_')} )
  secKey <- lapply(mud, function(x) { paste('SEC',x,sep='_')} )
  bodKey <- lapply(mud, function(x) { paste('BOD',x,sep='_')} )
  endKey <- lapply(mud, function(x) { paste('END',x,sep='_')} )

  # TBD:: Note on the Opposite Strand we should reverse all framents, but currently fragLen==1 are left alone for effiecntcy...
  # Sketch Output::
  #
  # F_C_N    DesSeqN[CG]DesSeqN
  #
  #                  D2 22222
  #                N M1 1111
  #                N U1 1111
  # F_C_U    DesBscU[tG]DesBscU
  # F_O_U    DesBscU[Ca]DesBscU
  #             1111 1U N
  #             1111 1M N
  #            22222 2D
  #
  if (pr=='rs'||pr=='ch') {
    if (fr1Key=='F' && fr2Key=='F') {
      bufC <- 0
      bufO <- 0
      str <- glue::glue(
        "{fr1Key}_{co1Key}_II{TAB}",paste0(rep(" ", 61+bufC), collapse=''),"{tib[[tarKey$D]][1]}{tib[[secKey$D]][1]}{BNG}{tib[[bodKey$D]][1]}{tib[[endKey$D]][1]}{RET}",
        "{fr1Key}_{co1Key}_IM{TAB}",paste0(rep(" ", 59+bufC), collapse=''),"{tib[[nxbKey$M]][1]}{BNG}{tib[[tarKey$M]][1]}{tib[[secKey$M]][1]}{BNG}{tib[[bodKey$M]][1]}{RET}",
        "{fr1Key}_{co1Key}_IU{TAB}",paste0(rep(" ", 59+bufC), collapse=''),"{tib[[nxbKey$U]][1]}{BNG}{tib[[tarKey$U]][1]}{tib[[secKey$U]][1]}{BNG}{tib[[bodKey$U]][1]}{RET}",
        "{fr1Key}_{co1Key}_U {TAB}{tib[[bscKey$U]][1]}{RET}",
        "{fr1Key}_{co1Key}_M {TAB}{tib[[bscKey$M]][1]}{RET}",
        "{fr1Key}_{co1Key}_D {TAB}{tib[[bscKey$D]][1]}{RET}",
        "{fr1Key}_{co1Key}_N {TAB}{tib[[desSeq]][1]}{RET}",
        # "FwdSeq{TAB}{tib$Forward_Sequence[1]}{RET}",
        "{RET}",
        "{fr1Key}_{co2Key}_N {TAB}{cmpl(tib[[desSeq]][2])}{RET}",
        "{fr2Key}_{co2Key}_D {TAB}{Biostrings::reverse(tib[[bscKey$D]][2])}{RET}",
        "{fr2Key}_{co2Key}_M {TAB}{Biostrings::reverse(tib[[bscKey$M]][2])}{RET}",
        "{fr2Key}_{co2Key}_U {TAB}{Biostrings::reverse(tib[[bscKey$U]][2])}{RET}",
        "{fr2Key}_{co2Key}_IU{TAB}",paste0(rep(" ", 11-bufO), collapse=''),"{Biostrings::reverse(tib[[bodKey$U]][2])}{tib[[secKey$U]][2]}{BNG}{tib[[tarKey$U]][2]}{tib[[nxbKey$U]][2]}{RET}",
        "{fr2Key}_{co2Key}_IM{TAB}",paste0(rep(" ", 11-bufO), collapse=''),"{Biostrings::reverse(tib[[bodKey$M]][2])}{tib[[secKey$M]][2]}{BNG}{tib[[tarKey$M]][2]}{tib[[nxbKey$M]][2]}{RET}",
        "{fr2Key}_{co2Key}_II{TAB}",paste0(rep(" ", 10-bufO), collapse=''),"{tib[[endKey$D]][2]}{Biostrings::reverse(tib[[bodKey$D]][2])}{tib[[secKey$D]][2]}{BNG}{tib[[tarKey$D]][2]}{RET}",
        "{RET}")
    } else if (fr1Key=='R' && fr2Key=='R') {
      bufO <- 0
      bufC <- 0
      str <- glue::glue(
        "{fr2Key}_{co2Key}_II{TAB}",paste0(rep(" ", 61+bufO), collapse=''),"{tib[[tarKey$D]][2]}{tib[[secKey$D]][2]}{BNG}{tib[[bodKey$D]][2]}{tib[[endKey$D]][2]}{RET}",
        "{fr2Key}_{co2Key}_IM{TAB}",paste0(rep(" ", 59+bufO), collapse=''),"{tib[[nxbKey$M]][2]}{BNG}{tib[[tarKey$M]][2]}{tib[[secKey$M]][2]}{BNG}{tib[[bodKey$M]][2]}{RET}",
        "{fr2Key}_{co2Key}_IU{TAB}",paste0(rep(" ", 59+bufO), collapse=''),"{tib[[nxbKey$U]][2]}{BNG}{tib[[tarKey$U]][2]}{tib[[secKey$U]][2]}{BNG}{tib[[bodKey$U]][2]}{RET}",
        "{fr2Key}_{co2Key}_D {TAB}{tib[[bscKey$D]][2]}{RET}",
        "{fr2Key}_{co2Key}_N {TAB}{revCmp(tib[[desSeq]][2])}{RET}",
        # "{fr2Key}_{co2Key}_N {TAB}{tib[[desSeq]][2]}{RET}",
        "{RET}",
        # "{fr1Key}_{co1Key}_N {TAB}{cmpl(tib[[desSeq]][1])}{RET}",
        "{fr1Key}_{co1Key}_N {TAB}{Biostrings::reverse(tib[[desSeq]][1])}{RET}",
        "{fr1Key}_{co1Key}_D {TAB}{Biostrings::reverse(tib[[bscKey$D]][1])}{RET}",
        "{fr1Key}_{co1Key}_IU{TAB}",paste0(rep(" ", 11-bufC), collapse=''),"{Biostrings::reverse(tib[[bodKey$U]][1])}{tib[[secKey$U]][1]}{BNG}{tib[[tarKey$U]][1]}{tib[[nxbKey$U]][1]}{RET}",
        "{fr1Key}_{co1Key}_IM{TAB}",paste0(rep(" ", 11-bufC), collapse=''),"{Biostrings::reverse(tib[[bodKey$M]][1])}{tib[[secKey$M]][1]}{BNG}{tib[[tarKey$M]][1]}{tib[[nxbKey$M]][1]}{RET}",
        "{fr1Key}_{co1Key}_II{TAB}",paste0(rep(" ", 10-bufC), collapse=''),"{tib[[endKey$D]][1]}{Biostrings::reverse(tib[[bodKey$D]][1])}{tib[[secKey$D]][1]}{BNG}{tib[[tarKey$D]][1]}{RET}",
        "{RET}")
    } else {
      stop(glue::glue("{RET}[{funcTag}]: ERROR: fr1Key={fr1Key}, fr2Key={fr2Key}, Allowed Values=[F,R]!{RET}{RET}"))
    }
  } else {
    if (fr1Key=='F' && fr2Key=='F') {
      buf <- 0
      str <- glue::glue(
        "{fr1Key}_{co1Key}_II{TAB}",paste0(rep(" ", 61+buf), collapse=''),"{tib[[tarKey$D]][1]}{tib[[secKey$D]][1]}{BNG}{tib[[bodKey$D]][1]}{tib[[endKey$D]][1]}{RET}",
        "{fr1Key}_{co1Key}_IM{TAB}",paste0(rep(" ", 59+buf), collapse=''),"{tib[[nxbKey$M]][1]}{BNG}{tib[[tarKey$M]][1]}{tib[[secKey$M]][1]}{BNG}{tib[[bodKey$M]][1]}{RET}",
        "{fr1Key}_{co1Key}_IU{TAB}",paste0(rep(" ", 59+buf), collapse=''),"{tib[[nxbKey$U]][1]}{BNG}{tib[[tarKey$U]][1]}{tib[[secKey$U]][1]}{BNG}{tib[[bodKey$U]][1]}{RET}",
        "{fr1Key}_{co1Key}_D {TAB}{tib[[bscKey$D]][1]}{RET}",
        "{fr1Key}_{co1Key}_N {TAB}{tib[[desSeq]][1]}{RET}",
        "{RET}",
        "{fr1Key}_{co2Key}_N {TAB}{cmpl(tib[[desSeq]][2])}{RET}",
        "{fr2Key}_{co2Key}_D {TAB}{Biostrings::reverse(tib[[bscKey$D]][2])}{RET}",
        "{fr2Key}_{co2Key}_IU{TAB}",paste0(rep(" ", 11-buf), collapse=''),"{Biostrings::reverse(tib[[bodKey$U]][2])}{tib[[secKey$U]][2]}{BNG}{tib[[tarKey$U]][2]}{tib[[nxbKey$U]][2]}{RET}",
        "{fr2Key}_{co2Key}_IM{TAB}",paste0(rep(" ", 11-buf), collapse=''),"{Biostrings::reverse(tib[[bodKey$M]][2])}{tib[[secKey$M]][2]}{BNG}{tib[[tarKey$M]][2]}{tib[[nxbKey$M]][2]}{RET}",
        "{fr2Key}_{co2Key}_II{TAB}",paste0(rep(" ", 10-buf), collapse=''),"{tib[[endKey$D]][2]}{Biostrings::reverse(tib[[bodKey$D]][2])}{tib[[secKey$D]][2]}{BNG}{tib[[tarKey$D]][2]}{RET}",
        "{RET}")
    } else if (fr1Key=='R' && fr2Key=='R') {
      buf <- 0
      buf <- 1
      str <- glue::glue(
        "{fr2Key}_{co2Key}_II{TAB}",paste0(rep(" ", 61+buf), collapse=''),"{tib[[tarKey$D]][2]}{BNG}{tib[[secKey$D]][2]}{tib[[bodKey$D]][2]}{tib[[endKey$D]][2]}{RET}",
        "{fr2Key}_{co2Key}_IM{TAB}",paste0(rep(" ", 60+buf), collapse=''),"{tib[[nxbKey$M]][2]}{tib[[tarKey$M]][2]}{BNG}{tib[[secKey$M]][2]}{tib[[bodKey$M]][2]}{RET}",
        "{fr2Key}_{co2Key}_IU{TAB}",paste0(rep(" ", 60+buf), collapse=''),"{tib[[nxbKey$U]][2]}{tib[[tarKey$U]][2]}{BNG}{tib[[secKey$U]][2]}{tib[[bodKey$U]][2]}{RET}",
        "{fr2Key}_{co2Key}_D {TAB}{tib[[bscKey$D]][2]}{RET}",
        "{fr2Key}_{co2Key}_N {TAB}{revCmp(tib[[desSeq]][2])}{RET}",
        # "{fr2Key}_{co2Key}_N {TAB}{tib[[desSeq]][2]}{RET}",
        "{RET}",
        # "{fr1Key}_{co1Key}_N {TAB}{cmpl(tib[[desSeq]][1])}{RET}",
        "{fr1Key}_{co1Key}_N {TAB}{Biostrings::reverse(tib[[desSeq]][1])}{RET}",
        "{fr1Key}_{co1Key}_D {TAB}{Biostrings::reverse(tib[[bscKey$D]][1])}{RET}",
        "{fr1Key}_{co1Key}_IU{TAB}",paste0(rep(" ", 11+buf), collapse=''),"{Biostrings::reverse(tib[[bodKey$U]][1])}{BNG}{tib[[secKey$U]][1]}{tib[[tarKey$U]][1]}{BNG}{tib[[nxbKey$U]][1]}{RET}",
        "{fr1Key}_{co1Key}_IM{TAB}",paste0(rep(" ", 11+buf), collapse=''),"{Biostrings::reverse(tib[[bodKey$M]][1])}{BNG}{tib[[secKey$M]][1]}{tib[[tarKey$M]][1]}{BNG}{tib[[nxbKey$M]][1]}{RET}",
        "{fr1Key}_{co1Key}_II{TAB}",paste0(rep(" ", 10+buf), collapse=''),"{tib[[endKey$D]][1]}{Biostrings::reverse(tib[[bodKey$D]][1])}{BNG}{tib[[secKey$D]][1]}{tib[[tarKey$D]][1]}{RET}",
        "{RET}")
    } else {
      stop(glue::glue("{RET}[{funcTag}]: ERROR: fr1Key={fr1Key}, fr2Key={fr2Key}, Allowed Values=[F,R]!{RET}{RET}"))
    }
  }
  if (verbose>=vt) cat(str)
  
  str
}

prbsToStrMUD = function(tib, mu='U',
                        verbose=0,vt=1,tc=1,tabsStr='') {
  funcTag <- 'prbsToStr'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  fr1Key <- tib$StrandFR[1]
  fr2Key <- tib$StrandFR[2]
  
  co1Key <- tib$StrandCO[1]
  co2Key <- tib$StrandCO[2]
  
  bscKey <- paste('DesBsc',mu, sep='')
  
  nxbKey <- paste('NXB',mu, sep='_')
  cpnKey <- paste('CPN',mu, sep='_')
  tarKey <- paste('TAR',mu, sep='_')
  secKey <- paste('SEC',mu, sep='_')
  bodKey <- paste('BOD',mu, sep='_')
  endKey <- paste('END',mu, sep='_')
  
  str <- glue::glue(
    paste0(rep(" ", 59+10), collapse=''),"{tib[[cpnKey]][1]}{tib[[secKey]][1]}{BNG}{tib[[bodKey]][1]}{tib[[endKey]][1]}{RET}",
    paste0(rep(" ", 59+8),  collapse=''),"{tib[[nxbKey]][1]}{BNG}{tib[[cpnKey]][1]}{tib[[secKey]][1]}{BNG}{tib[[bodKey]][1]}{RET}",
    "{fr1Key}_{co1Key}_{mu}{TAB}{tib[[bscKey]][1]}{RET}",
    "{fr2Key}_{co2Key}_{mu}{TAB}{Biostrings::reverse(tib[[bscKey]][2])}{RET}",
    paste0(rep(" ", 9+10), collapse=''),"{Biostrings::reverse(tib[[bodKey]][2])}{tib[[secKey]][2]}{BNG}{tib[[cpnKey]][2]}{tib[[nxbKey]][2]}{RET}",
    paste0(rep(" ", 9+9),  collapse=''),"{tib[[endKey]][2]}{Biostrings::reverse(tib[[bodKey]][2])}{tib[[secKey]][2]}{BNG}{tib[[cpnKey]][2]}{RET}",
    "{RET}")
  
  if (verbose>=vt) cat(str)
  
  str
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                Infinium Methylation Probe Design Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

desAllPrbs = function(tib) {
  fr <- tib %>% dplyr::distinct(FR) %>% base::as.logical()
  co <- tib %>% dplyr::distinct(CO) %>% base::as.logical()
  pr <- tib %>% dplyr::distinct(PRB_DES) %>% base::as.character()
  
  dplyr::bind_rows(tib %>% 
                     des2prbs(fwd=fr, con=co, pr=pr, mu='U', strand='DesBscU') %>%
                     des2prbs(fwd=fr, con=co, pr=pr, mu='M', strand='DesBscM') %>%
                     des2prbs(fwd=fr, con=co, pr=pr, mu='D', strand='DesBscD') )
}

des2prbs = function(tib, fwd, con, pr, mu, strand='DesSeqN', len=48, del='_',QC_CPN=TRUE,
                    verbose=0,vt=1,tc=1,tabsStr='') {
  funcTag <- 'des2prbs'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  stopifnot(is.logical(fwd))
  stopifnot(is.logical(con))
  
  if (mu!='N' && mu!='U' && mu!='M' && mu!='D')
    stop(glue::glue("{RET}[{funcTag}]: ERROR: mu={mu} Only Supported=[N,U,M,D]!{RET}{RET}"))
  if (pr!='cg' && pr!='ch' & pr!='rs' & pr!='rp')
    stop(glue::glue("{RET}[{funcTag}]: ERROR: pr={pr} Only Supported=[cg,ch,rp,rs]!{RET}{RET}"))
  
  strand <- rlang::sym(strand)
  if (pr=='rs'|| pr=='ch') {
    # $prb_F_C_I   = revCmp(substr($des_F_C, 60, 50));
    # $prb_R_C_I   = revCmp(substr($des_R_C, 61, 50));
    # $prb_F_C_II  = revCmp(substr($des_F_C, 60, 50));
    # $prb_R_C_II  = revCmp(substr($des_R_C, 61, 50));
    #
    # $prb_F_O_I   = revCmp(substr($des_F_O, 61, 50));
    # $prb_R_O_I   = revCmp(substr($des_R_O, 62, 50));
    # $prb_F_O_II  = revCmp(substr($des_F_O, 61, 50));
    # $prb_R_O_II  = revCmp(substr($des_R_O, 62, 50));
    if      ( fwd &&  con) nxb_pos <- 60
    else if (!fwd &&  con) nxb_pos <- 61
    else if ( fwd && !con) nxb_pos <- 61
    else if (!fwd && !con) nxb_pos <- 60
    # else if (!fwd && !con) nxb_pos <- 62  # Original 
    else {
      stop(glue::glue("{RET}[{funcTag}]: ERROR: unsupported combination fwd={fwd}, con={con}!!!{RET}{RET}"))
    }
    
  } else if (pr=='cg' || pr=='rp') {
    # $prb_F_C_I  = revCmp(substr($des_F_C, 60, 50));
    # $prb_F_O_I  = revCmp(substr($des_F_O, 61, 50));
    # $prb_R_C_I  = revCmp(substr($des_R_C, 60, 50));
    # $prb_R_O_I  = revCmp(substr($des_R_O, 61, 50));
    # 
    # $prb_F_C_II  = revCmp(substr($des_F_C, 61, 50));
    # $prb_F_O_II  = revCmp(substr($des_F_O, 62, 50));
    # $prb_R_C_II  = revCmp(substr($des_R_C, 61, 50));
    # $prb_R_O_II  = revCmp(substr($des_R_O, 62, 50));
    nxb_pos <- 60
    if (!con) nxb_pos <- 61

  } else {
    stop(glue::glue("[{funcTag}]: ERROR: NOT READY!!!{RET}{RET}"))
  }
  cpg_pos <- nxb_pos + 1
  sec_pos <- cpg_pos + 1
  bod_pos <- sec_pos + 1
  end_pos <- bod_pos + len
  
  # Special consideration is needed for U/M strands at the query site. 
  #  For CN (i.e. cg or ch) this is actually done naturally in U/M conversion
  #  However, for  non-CN probes (i.e. rs) this needs to be forced to U/M
  #
  # This is handled by the TAR (Target/Query Nucleotide). This should only change
  #  for U/M (QMAP_U/QMAP_M) for D its just itself.
  #
  tib <- tib %>% dplyr::mutate(
    NXB=stringr::str_sub(!!strand, nxb_pos, nxb_pos),
    CPN=stringr::str_sub(!!strand, cpg_pos, cpg_pos),
    TAR=qmaps(CPN, mu=mu),
    SEC=stringr::str_sub(!!strand, sec_pos, sec_pos),
    BOD=stringr::str_sub(!!strand, bod_pos, end_pos-1),
    END=stringr::str_sub(!!strand, end_pos, end_pos)
  )

  #  QC TEST:: for CpN (cg or ch) verify that the probes are equal. Well call this
  #   PRB0 (CGN) and PRB1 (TAR). After testing remove PRB0
  #
  if (QC_CPN && (pr=='cg')) {
    tib <- tib %>%
      tidyr::unite(PRB0, CPN,SEC,BOD, sep='', remove=FALSE) %>%
      tidyr::unite(PRB1, TAR,SEC,BOD, sep='', remove=FALSE) %>%
      tidyr::unite(PRB2, SEC,BOD,END, sep='', remove=FALSE) %>%
      dplyr::mutate(PRB0=revCmp(PRB0), PRB1=revCmp(PRB1), PRB2=revCmp(PRB2))
    
    qc_tib <- tib %>% filter(PRB0!=PRB1)
    qc_len <- qc_tib %>% base::nrow()
    if (qc_len != 0) {
      qc_tib %>% dplyr::select(1,PRB0,PRB1) %>% print()
      stop(glue::glue("{RET}[{funcTag}]: ERROR: pr={pr}, qc_len={qc_len} != 0!!!{RET}{RET}"))
    }
  } else {
    tib <- tib %>%
      tidyr::unite(PRB1, TAR,SEC,BOD, sep='', remove=FALSE) %>%
      tidyr::unite(PRB2, SEC,BOD,END, sep='', remove=FALSE) %>%
      dplyr::mutate(PRB1=revCmp(PRB1), PRB2=revCmp(PRB2))
  }

  # Add suffix to sequences for merging later
  tib <- tib %>%
    dplyr::select(PRB1,PRB2, NXB,CPN,TAR,BOD,END, everything()) %>%
    dplyr::rename(!!paste('PRB1',mu, sep=del):=PRB1,
                  !!paste('PRB2',mu, sep=del):=PRB2,
                  !!paste('NXB', mu, sep=del):=NXB,
                  !!paste('CPN', mu, sep=del):=CPN,
                  !!paste('TAR', mu, sep=del):=TAR,
                  !!paste('SEC', mu, sep=del):=SEC,
                  !!paste('BOD', mu, sep=del):=BOD,
                  !!paste('END', mu, sep=del):=END)
  
  tib
}

des2prbsNOTES = function(srd, desSeq,
                         vt=1, tc=0, tabsStr='') {
  funcTag <- 'des2prbs'
  for (i in c(1:tc)) tabsStr <- paste0(tabsStr, TAB)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  # bscNewU <- bscNewU(desSeq)
  # bscNewM <- bscNewM(desSeq)
  # bscNewD <- bscNewD(desSeq)
  
  # my @desSetU = unpack("A11"."A48"."AAAA"."A47"."A12", $bscNewU);
  # my @desSetM = unpack("A11"."A48"."AAAA"."A47"."A12", $bscNewM);
  # my @desSetD = unpack("A11"."A48"."AAAA"."A47"."A12", $bscNewD);
  
  #    my @desSetU=unpack("A11"."A48"."AAAA"."A47"."A12",bscNewU($desSeq,$retUC));                                                                                                             
  #    my @desSetM=unpack("A11"."A48"."AAAA"."A47"."A12",bscNewM($desSeq,$retUC));                                                                                                             
  
  # my ($desNxbU, $desCpgU, $desSeqU, $desEndU);
  # my ($desNxbM, $desCpgM, $desSeqM, $desEndM);
  # my ($desNxbD, $desCpgD, $desSeqD, $desEndD);
  
  #         (  $desNxbU,    $desCpgU,               $desSeqU,                $desEndU) =                                                                                                       
  #   return( $desSetU[2], $desSetU[3], $desSetU[4].$desSetU[5].$desSetU[6], $desSetU[7],                                                                                                      
  #           $desSetM[2], $desSetM[3], $desSetM[4].$desSetM[5].$desSetM[6], $desSetM[7]) if ($desCO eq $C);                                                                                   
  #          ( $desNxbM,    $desCpgM,               $desSeqM,                $desEndM) =                                                                                                       
  
  #                  ( $desNxbU,            $desCpgU,                    $desSeqU,                  $desEndU) =                                                                                
  #   return( revCmpl($desSetU[4]), revCmpl($desSetU[3]), revCmpl($desSetU[1].$desSetU[2]), revCmpl($desSetU[0]),                                                                              
  #           revCmpl($desSetM[4]), revCmpl($desSetM[3]), revCmpl($desSetM[1].$desSetM[2]), revCmpl($desSetM[0])) if ($desCO eq $O);                                                           
  #                  ( $desNxbM,            $desCpgM,                    $desSeqM,                  $desEndM) =                                                                                
  
  
  # $$prbRef[$srd][$iU]=[ $desSetU[2], $desSetU[3], $desSetU[4].$desSetU[5].$desSetU[6], $desSetU[7] ];
  # $$prbRef[$srd][$iM]=[ $desSetM[2], $desSetM[3], $desSetM[4].$desSetM[5].$desSetM[6], $desSetM[7] ];
  # $$prbRef[$srd][$iD]=[ $desSetD[2], $desSetD[3], $desSetD[4].$desSetD[5].$desSetD[6], $desSetD[7] ];
  # 
  # $srd++;
  # $$prbRef[$srd][$iU]=[ revCmpl($desSetU[4]), revCmpl($desSetU[3]), revCmpl($desSetU[1].$desSetU[2]), revCmpl($desSetU[0]) ];
  # $$prbRef[$srd][$iM]=[ revCmpl($desSetM[4]), revCmpl($desSetM[3]), revCmpl($desSetM[1].$desSetM[2]), revCmpl($desSetM[0]) ];
  # $$prbRef[$srd][$iD]=[ revCmpl($desSetD[4]), revCmpl($desSetD[3]), revCmpl($desSetD[1].$desSetD[2]), revCmpl($desSetD[0]) ];
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}{RET}"))
  
  NULL
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                         Output improbe Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
writeImprobeInput = function(tib, name, dir, run=FALSE, 
                             exe=NULL, impTango=NULL, imp13mer=NULL, 
                             tbVar='BOTH', coVar='BOTH',
                             verbose=0,vt=3,tc=1,tabsStr='') {
  funcTag <- 'writeImprobeInput'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} name={name}, dir={dir}.{RET}"))
  
  if (!dir.exists(dir)) dir.create(dir, recursive=TRUE)
  
  imp_run_sh  <- file.path(dir, paste0(name,'.improbe-input.sh'))
  imp_inp_tsv <- file.path(dir, paste0(name,'.improbe-input.tsv'))
  imp_out_tsv <- file.path(dir, paste0(name,'.improbe-output.tsv'))
  imp_log_txt <- file.path(dir, paste0(name,'.improbe-output.log'))
  imp_fin_txt <- file.path(dir, paste0(name,'.improbe-output.fin.txt'))
  imp_cmd_str <- ''

  readr::write_tsv(tib, imp_inp_tsv)
  
  if (!is.null(exe) && ! is.null(impTango) && !is.null(imp13mer)) {
    # stopifnot(file.exists(exe))
    # stopifnot(file.exists(impTango))
    # stopifnot(file.exists(imp13mer))
    
    # ${CMD} -oASPE -tBOTH -cBoth -n${fn_13mer} -a${fn_tango} -V ${fn_in} >$fn_out 2>${fn_log}
    imp_cmd_str <- glue::glue("{exe} -oASPE -t{tbVar} -c{coVar} -n{imp13mer} -a{impTango} -V ",
                              "{imp_inp_tsv} >{imp_out_tsv} 2>{imp_log_txt}{RET}",
                              "touch {imp_fin_txt}{RET}")
    
    readr::write_file(imp_cmd_str, imp_run_sh)
    Sys.chmod(imp_run_sh, mode='0777')
    
    if (run) system(imp_run_sh, wait = TRUE)
  }
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}{RET}"))
  
  imp_out_tsv
}
# 
# tibToImprobeInput = function(tib, verbose=0,vt=5,tc=1,tabsStr='') {
#   funcTag <- 'tibToImprobeInput'
#   tabsStr <- paste0(rep(TAB, tc), collapse='')
#   
#   if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
#   
#   imp_tib <- tib %>% dplyr::mutate(Sequence=stringr::str_replace(Forward_Sequence, '\\[[A-Za-z][A-Za-z]\\]', '[CG]'), 
#                                    CpG_Island='FALSE' ) %>% 
#     dplyr::select(IlmnID, Sequence, Genome_Build, CHR, MAPINFO, CpG_Island) %>% 
#     dplyr::rename(Seq_ID=IlmnID, Chromosome=CHR, Coordinate=MAPINFO)
#   
#   imp_tib
# }

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                     Load/Format Input File Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
loadManifestRS = function(file, swap, revAllele=FALSE, verbose=0,vt=1,tc=1,tabsStr='') {
  funcTag <- 'loadManifestRS'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Loading manifest(RS)={file}.{RET}"))
  
  snp_man_tib <- suppressMessages(suppressWarnings(readr::read_csv(file) ))
  snp_swp_tib <- suppressMessages(suppressWarnings(readr::read_tsv(swap) ))
  
  # NOTE: Rediculous Swapping Needed for Infinium I M/U -> U/M from old manufacturing mistake::
  #  - This method should not be used...
  if (revAllele) {
    snp_man_tib <- snp_man_tib %>% 
      dplyr::mutate(
        TMP_ADD_A=AddressA_ID,
        TMP_ADD_B=AddressB_ID,
        TMP_SEQ_A=AlleleA_ProbeSeq,
        TMP_SEQ_B=AlleleB_ProbeSeq,
        IS_SWAP=case_when(
          Infinium_Design_Type=='I' & stringr::str_ends(AlleleA_ProbeSeq,'C') ~ TRUE,
          TRUE ~ FALSE
        ),
        AddressA_ID=case_when( IS_SWAP ~ TMP_ADD_B, TRUE ~ AddressA_ID),
        AddressB_ID=case_when( IS_SWAP ~ TMP_ADD_A, TRUE ~ AddressB_ID),
        AlleleA_ProbeSeq=case_when( IS_SWAP ~ TMP_SEQ_B, TRUE ~ AlleleA_ProbeSeq ),
        AlleleB_ProbeSeq=case_when( IS_SWAP ~ TMP_SEQ_A, TRUE ~ AlleleB_ProbeSeq )
      )
  }
  
  snp_man_tib <-  snp_man_tib %>% 
    dplyr::mutate(CGN=stringr::str_remove(Name, '_\\.*\\$') ) %>%
    dplyr::select(IlmnID, Name, AddressA_ID, AlleleA_ProbeSeq, AddressB_ID, AlleleB_ProbeSeq, 
                  Infinium_Design_Type, Next_Base, Color_Channel, CGN)
  
  snp_swp_tib <- snp_swp_tib %>% tidyr::separate(Seq_ID, into=c('IlmnID', 'diNUC'), sep='_') %>% 
    dplyr::rename(Forward_Sequence=Sequence, CHR=Chromosome, MAPINFO=Coordinate) %>% 
    dplyr::mutate(Forward_Sequence=stringr::str_replace(Forward_Sequence,'\\[CG\\]', paste0('[',diNUC,']') ) ) %>% 
    dplyr::select(-CpG_Island)
  
  snp_man_tib <- snp_man_tib %>% dplyr::inner_join(snp_swp_tib, by="IlmnID") %>%
    dplyr::select(IlmnID, Name, AddressA_ID, AlleleA_ProbeSeq, AddressB_ID, AlleleB_ProbeSeq,
                  Infinium_Design_Type, Next_Base, Color_Channel, Forward_Sequence, 
                  Genome_Build, CHR, MAPINFO, CGN, diNUC)
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}{RET}"))
  
  snp_man_tib
}                          

loadManifestCG = function(file, pr='cg', verbose=0,vt=1,tc=1,tabsStr='') {
  funcTag <- 'loadManifestCG'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Loading manifest(CG)={file}.{RET}"))
  
  cpg_man_tib <- suppressMessages(suppressWarnings(readr::read_csv(file, skip=7) )) %>%
    dplyr::filter(stringr::str_starts(IlmnID,pr)) %>%
    dplyr::mutate(CGN=stringr::str_remove(Name, '_\\.*\\$'),
                  diNUC=stringr::str_remove(
                    stringr::str_replace(Forward_Sequence, '^.*\\[([A-Za-z][A-Za-z])].*$', "\\$1"),'^\\\\') ) %>%
    dplyr::select(IlmnID, Name, AddressA_ID, AlleleA_ProbeSeq, AddressB_ID, AlleleB_ProbeSeq,
                  Infinium_Design_Type, Next_Base, Color_Channel, Forward_Sequence, 
                  Genome_Build, CHR, MAPINFO, CGN, diNUC)
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}{RET}"))
  
  cpg_man_tib
}                          

loadManifestCH = function(file, pr='ch', ry='R', verbose=0,vt=1,tc=1,tabsStr='') {
  funcTag <- 'loadManifestCH'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Loading manifest(CH)={file}.{RET}"))
  
  cph_man_tib <- suppressMessages(suppressWarnings(readr::read_csv(cph_man_csv ))) %>% 
    dplyr::filter(stringr::str_starts(IlmnID,pr)) %>%
    dplyr::mutate(Strand=stringr::str_sub(IlmnID, -1), FR=Strand, CO='O',
                  TP=case_when( is.na(AlleleB_ProbeSeq) ~ 'II', TRUE ~ 'I' )) %>%
    dplyr::mutate(CGN=stringr::str_remove(Name, '_\\.*\\$'),
                  diNUC=paste0(ry,stringr::str_remove(
                    stringr::str_replace(Forward_Sequence, '^.*\\[[A-Za-z]([A-Za-z])].*$', "\\$1"),'^\\\\')),
                  Forward_Sequence=stringr::str_replace(Forward_Sequence, '\\[[A-Za-z][A-Za-z]\\]', paste0('[',diNUC,']')) ) %>%
    dplyr::select(IlmnID, Name, AddressA_ID, AlleleA_ProbeSeq, AddressB_ID, AlleleB_ProbeSeq,
                  Infinium_Design_Type, Next_Base, Color_Channel, Forward_Sequence, 
                  Genome_Build, CHR, MAPINFO, CGN, diNUC)
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}{RET}"))
  
  cph_man_tib
}

loadManifestYH = function(file, pr='ch', verbose=0,vt=1,tc=1,tabsStr='') {
  funcTag <- 'loadManifestCH'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Loading manifest(CH)={file}.{RET}"))
  
  cph_man_tib <- suppressMessages(suppressWarnings(readr::read_csv(cph_man_csv ))) %>% 
    dplyr::filter(stringr::str_starts(IlmnID,pr)) %>%
    dplyr::mutate(Strand=stringr::str_sub(IlmnID, -1), FR=Strand, CO='O',
                  TP=case_when( is.na(AlleleB_ProbeSeq) ~ 'II', TRUE ~ 'I' )) %>%
    dplyr::mutate(CGN=stringr::str_remove(Name, '_\\.*\\$'),
                  diNUC=paste0('Y',stringr::str_remove(
                    stringr::str_replace(Forward_Sequence, '^.*\\[[A-Za-z]([A-Za-z])].*$', "\\$1"),'^\\\\')),
                  Forward_Sequence=stringr::str_replace(Forward_Sequence, '\\[[A-Za-z][A-Za-z]\\]', paste0('[',diNUC,']')) ) %>%
    dplyr::select(IlmnID, Name, AddressA_ID, AlleleA_ProbeSeq, AddressB_ID, AlleleB_ProbeSeq,
                  Infinium_Design_Type, Next_Base, Color_Channel, Forward_Sequence, 
                  Genome_Build, CHR, MAPINFO, CGN, diNUC)
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}{RET}"))
  
  cph_man_tib
}                          

loadManifestCH_OLD = function(file, pr='ch', verbose=0,vt=1,tc=1,tabsStr='') {
  funcTag <- 'loadManifestCH_old'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Loading manifest(CH)={file}.{RET}"))
  
  cph_man_tib <- suppressMessages(suppressWarnings(readr::read_csv(cph_man_csv ))) %>% 
    dplyr::filter(stringr::str_starts(IlmnID,pr)) %>%
    dplyr::mutate(Strand=stringr::str_sub(IlmnID, -1), FR=Strand, CO='O',
                  TP=case_when( is.na(AlleleB_ProbeSeq) ~ 'II', TRUE ~ 'I' )) %>%
    dplyr::mutate(CGN=stringr::str_remove(Name, '_\\.*\\$'),
                  diNUC=stringr::str_remove(
                    stringr::str_replace(Forward_Sequence, '^.*\\[([A-Za-z][A-Za-z])].*$', "\\$1"),'^\\\\') ) %>%
    dplyr::select(IlmnID, Name, AddressA_ID, AlleleA_ProbeSeq, AddressB_ID, AlleleB_ProbeSeq,
                  Infinium_Design_Type, Next_Base, Color_Channel, Forward_Sequence, 
                  Genome_Build, CHR, MAPINFO, CGN, diNUC)
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}{RET}"))
  
  cph_man_tib
}                          

loadImprobeDesign = function(file=NULL, src_des_tib=NULL, verbose=0,vt=1,tc=1) {
  funcTag <- 'loadImprobeDesign'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Loading improbe={file}.{RET}"))

  if (!is.null(file)) src_des_tib <- suppressMessages(suppressWarnings(readr::read_tsv(file)))
  stopifnot(!is.null(src_des_tib))
  
  src_des_tib <- src_des_tib %>%
    dplyr::rename(PRB1_U=UnMethyl_Probe_Sequence,
                  PRB1_M=Methyl_Probe_Sequence,
                  NXB_U=UnMethyl_Next_Base,
                  NXB_M=Methyl_Next_Base) %>%
    dplyr::mutate(Probe_ID=paste(Seq_ID, Methyl_Allele_FR_Strand, stringr::str_sub(Methyl_Allele_TB_Strand,1,1), 
                                 Methyl_Allele_CO_Strand, sep='_'),
                  FR=case_when(Methyl_Allele_FR_Strand=='F'   ~ TRUE, Methyl_Allele_FR_Strand=='R'   ~ FALSE, TRUE ~ NA),
                  TB=case_when(Methyl_Allele_TB_Strand=='TOP' ~ TRUE, Methyl_Allele_TB_Strand=='BOT' ~ FALSE, TRUE ~ NA),
                  CO=case_when(Methyl_Allele_CO_Strand=='C'   ~ TRUE, Methyl_Allele_CO_Strand=='O'   ~ FALSE, TRUE ~ NA),
                  diNUC='CG',
                  NXB_IMP=case_when(NXB_U==NXB_M ~ NXB_U,
                                    TRUE ~ NA_character_),
                  COL_U=case_when(NXB_U=='A'|NXB_U=='T' ~ FALSE, # Red == FALSE
                                  NXB_U=='C'|NXB_U=='G' ~ TRUE,  # Grn == TRUE
                                  TRUE ~ NA),
                  COL_M=case_when(NXB_M=='A'|NXB_M=='T' ~ FALSE, # Red == FALSE
                                  NXB_M=='C'|NXB_M=='G' ~ TRUE,  # Grn == TRUE
                                  TRUE ~ NA),
                  COL_IMP=case_when(COL_U==COL_M ~ COL_U,
                                    TRUE ~ NA)
                  
    ) %>%
    
    # Design Score Parameters::
    dplyr::rename(
      PRB_SCR_U=UnMethyl_Final_Score,
      PRB_SCR_M=Methyl_Final_Score,
      PRB_SCR_S=Probeset_Score,
      
      TM_RAW_M=Methyl_Tm,
      TM_SCR_M=Methyl_Tm_Score,
      GC_RAW_M=Methyl_GC_Percent,
      GC_SCR_M=Methyl_GC_Score,
      KM_RAW_M=Methyl_13mer_Count,
      KM_SCR_M=Methyl_13mer_Score,
      AD_RAW_M=Methyl_Address_Count,
      AD_SCR_M=Methyl_Address_Score,
      CM_RAW_M=Methyl_Self_Complementarity,
      CM_SCR_M=Methyl_Self_Complementarity_Score,
      MO_RAW_M=Methyl_Mono_Run,
      MO_SCR_M=Methyl_Mono_Run_Score,
      EP_RAW_M=Methyl_Ectopic_Count,
      EP_SCR_M=Methyl_Ectopic_Score,
      CG_RAW_M=Methyl_Underlying_CpG_Count,
      MD_RAW_M=Methyl_Underlying_CpG_Min_Dist,
      CG_SCR_M=Methyl_Underlying_CpG_Score,
      NB_SCR_M=Methyl_Next_Base_Score,
      
      TM_RAW_U=UnMethyl_Tm,
      TM_SCR_U=UnMethyl_Tm_Score,
      GC_RAW_U=UnMethyl_GC_Percent,
      GC_SCR_U=UnMethyl_GC_Score,
      KM_RAW_U=UnMethyl_13mer_Count,
      KM_SCR_U=UnMethyl_13mer_Score,
      AD_RAW_U=UnMethyl_Address_Count,
      AD_SCR_U=UnMethyl_Address_Score,
      CM_RAW_U=UnMethyl_Self_Complementarity,
      CM_SCR_U=UnMethyl_Self_Complementarity_Score,
      MO_RAW_U=UnMethyl_Mono_Run,
      MO_SCR_U=UnMethyl_Mono_Run_Score,
      EP_RAW_U=UnMethyl_Ectopic_Count,
      EP_SCR_U=UnMethyl_Ectopic_Score,
      CG_RAW_U=UnMethyl_Underlying_CpG_Count,
      MD_RAW_U=UnMethyl_Underlying_CpG_Min_Dist,
      CG_SCR_U=UnMethyl_Underlying_CpG_Score,
      NB_SCR_U=UnMethyl_Next_Base_Score
    ) %>%
    dplyr::select(Probe_ID, Seq_ID, FR, TB, CO, diNUC, NXB_IMP, COL_IMP, PRB1_U, PRB1_M, 
                  Forward_Sequence, 
                  Genome_Build, Chromosome, Coordinate,
                  dplyr::contains("_RAW_"), dplyr::contains("_SCR_"))
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}{RET}"))
  
  src_des_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                     Basic Bisulfite Conversion Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
shearBrac = function(x) {
  x %>% stringr::str_remove('\\[') %>% stringr::str_remove('\\]')
}

bscU = function(x, uc=FALSE) {
  x <- stringr::str_to_upper(x)
  if (uc) x <- tr(x, 'CYSMBHV', 'TTKWKWD')
  else    x <- tr(x, 'CYSMBHV', 'ttkwkwd')
  x
}
bscUs = function(x, uc=FALSE) { bscU(x) }

MAPM = function(x) {
  if (length(MAP_M[[x]])==0) return(x)
  MAP_M[[x]]
}
bscM = function(x) { stringr::str_replace_all(x, '([CYSMBHV][GRSKBDV])', MAPM) }
bscMs = function(x, uc=FALSE) {
  x <- stringr::str_to_upper(x)
  x <- lapply(x, bscM) %>% BiocGenerics::unlist()
  x <- tr(x, 'CYSMBHV', 'ttkwkwd')
  if (uc) x <- stringr::str_to_upper(x)
  x
}

MAPD = function(x) {
  if (length(MAP_D[[x]])==0) return(x)
  MAP_D[[x]]
}
bscD = function(x) { stringr::str_replace_all(x, '([CYSMBHV][GRSKBDV])', MAPD) }
bscDs = function(x, uc=FALSE) {
  x <- stringr::str_to_upper(x)
  x <- lapply(x, bscD) %>% BiocGenerics::unlist()
  x <- tr(x, 'CYSMBHV', 'ttkwkwd')
  if (uc) x <- stringr::str_to_upper(x)
  x
}

QMAP = function(x, mu) {
  if (mu=='U') {
    return(QMAP_U[[x]])
  } else if (mu=='M') {
    return(QMAP_M[[x]])
  }
  x
}
qmaps = function(x, mu) {
  x <- lapply(x, QMAP, mu) %>% BiocGenerics::unlist()
}

cmpIUPAC = function(x) {
  if (base::is.element(x, names(IUPAC_EQ) )) return(IUPAC_EQ[[x]])
  # if (is.null(IUPAC_EQ[[x]])) return(FALSE)
  # if (length(IUPAC_EQ[[x]])==0) return(FALSE)
  FALSE
}

cmpIUPACs = function(x) {
  x <- lapply(x, cmpIUPAC) %>% BiocGenerics::unlist()
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                      MisMatch Probe Comparison Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

cmpInfIMU_MisMatch = function(tib, fieldA, fieldB, mu, del='_',
                              verbose=0,vt=4,tc=1,tabsStr='') {
  funcTag <- 'cmpInfIMU_MisMatch'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]: fieldA={fieldA}, fieldB={fieldB} mu={mu}{RET}"))
  
  fieldA <- rlang::sym(fieldA)
  fieldB <- rlang::sym(fieldB)
  
  tib <- tib %>% 
    dplyr::mutate(
      BOD_NumMM=mapply(
        adist,
        stringr::str_sub(stringr::str_to_upper(!!fieldA),1,stringr::str_length(!!fieldA)-1),
        stringr::str_sub(stringr::str_to_upper(!!fieldB),1,stringr::str_length(!!fieldB)-1) ),
      
      DI_NUC_AB=paste0(
        stringr::str_to_upper(stringr::str_sub(!!fieldA,stringr::str_length(!!fieldA),stringr::str_length(!!fieldA)) ),
        stringr::str_to_upper(stringr::str_sub(!!fieldB,stringr::str_length(!!fieldB),stringr::str_length(!!fieldB)) )
      ),
      TAR_EQU=cmpIUPACs(DI_NUC_AB)
    ) %>%
    dplyr::rename(!!paste('BOD_NumMM',mu, sep=del):=BOD_NumMM,
                  !!paste('TAR_EQU',  mu, sep=del):=TAR_EQU)
  
  tib
}

cmpInfI_MisMatch = function(tib, fieldAU, fieldBU, fieldAM, fieldBM, del='_',
                            verbose=0,vt=4,tc=1,tabsStr='') {
  funcTag <- 'cmpInfI_MisMatch'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]: fieldAU={fieldAU}, fieldBU={fieldBU}{RET}"))
  if (verbose>=vt) cat(glue::glue("[{funcTag}]: fieldAM={fieldAM}, fieldBM={fieldBM}{RET}"))

  tib <- tib %>% cmpInfIMU_MisMatch(fieldAU, fieldBU, mu='U', del=del,verbose=verbose, vt=vt+1,tc=tc+1)
  tib <- tib %>% cmpInfIMU_MisMatch(fieldAM, fieldBM, mu='M', del=del,verbose=verbose, vt=vt+1,tc=tc+1)
  tib <- tib %>% dplyr::mutate(
    Man_MisMatch=(BOD_NumMM_U+BOD_NumMM_M)/2, #, na.rm=TRUE),
    Man_TarMatch=case_when(TAR_EQU_U & TAR_EQU_M ~ TRUE, TRUE ~ FALSE) ) %>%
    dplyr::select(-c(BOD_NumMM_U,BOD_NumMM_M,TAR_EQU_U,TAR_EQU_M))
  
  tib
}

cmpInfII_MisMatch = function(tib, fieldA, fieldB, mu='D', del='_',
                    verbose=0,vt=4,tc=1,tabsStr='') {
  funcTag <- 'cmpInfI'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  fieldA <- rlang::sym(fieldA)
  fieldB <- rlang::sym(fieldB)
  
  tib <- tib %>% cmpInfIMU_MisMatch(fieldA, fieldB, mu='D', del=del,verbose=verbose, vt=vt+1,tc=tc+1) %>%
    dplyr::rename(
      Man_MisMatch=BOD_NumMM_D,
      Man_TarMatch=TAR_EQU_D)
    # dplyr::select(-c(BOD_NumMM_U,BOD_NumMM_M,TAR_EQU_U,TAR_EQU_M))

  tib
}


# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                       Exact Probe Comparison Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

cmpInfIMU= function(tib, fieldA, fieldB, mu, del='_',
                    verbose=0,vt=4,tc=1,tabsStr='') {
  funcTag <- 'cmpInfMU'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]: fieldA={fieldA}, fieldB={fieldB} mu={mu}{RET}"))
  
  fieldA <- rlang::sym(fieldA)
  fieldB <- rlang::sym(fieldB)
  
  tib <- tib %>% 
    dplyr::mutate(
      SUB_SEQ_A=stringr::str_to_upper(stringr::str_sub(!!fieldA,1,stringr::str_length(!!fieldA)-1)),
      SUB_SEQ_B=stringr::str_to_upper(stringr::str_sub(!!fieldB,1,stringr::str_length(!!fieldB)-1)),
      DI_NUC_AB=paste0(
        stringr::str_to_upper(stringr::str_sub(!!fieldA,stringr::str_length(!!fieldA),stringr::str_length(!!fieldA)) ),
        stringr::str_to_upper(stringr::str_sub(!!fieldB,stringr::str_length(!!fieldB),stringr::str_length(!!fieldB)) )
      ),
      
      BOD_EQU=case_when(SUB_SEQ_A==SUB_SEQ_B ~ TRUE, TRUE ~ FALSE),
      TAR_EQU=cmpIUPACs(DI_NUC_AB),
      Inf1_Match=case_when(BOD_EQU & BOD_EQU==TAR_EQU ~ TRUE, TRUE ~ FALSE)
    ) %>%
    dplyr::select(-c(SUB_SEQ_A,SUB_SEQ_B,BOD_EQU,DI_NUC_AB,TAR_EQU)) %>%
    dplyr::rename(!!paste('Inf1_Match',mu, sep=del):=Inf1_Match)
  
  if (verbose>=vt+1) print(tib)
  
  tib
}

cmpInfI = function(tib, fieldAU, fieldBU, fieldAM, fieldBM, del='_',
                   verbose=0,vt=4,tc=1,tabsStr='') {
  funcTag <- 'cmpInfI'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]: fieldAU={fieldAU}, fieldBU={fieldBU}{RET}"))
  if (verbose>=vt) cat(glue::glue("[{funcTag}]: fieldAM={fieldAM}, fieldBM={fieldBM}{RET}"))
  
  tib <- tib %>% cmpInfIMU(fieldAU, fieldBU, mu='U', del=del,verbose=verbose, vt=vt+1,tc=tc+1)
  tib <- tib %>% cmpInfIMU(fieldAM, fieldBM, mu='M', del=del,verbose=verbose, vt=vt+1,tc=tc+1)
  tib <- tib %>%
    dplyr::mutate(Man_Match=case_when(Inf1_Match_U & Inf1_Match_M ~ TRUE, TRUE ~ FALSE) ) %>%
    dplyr::select(-c(Inf1_Match_U,Inf1_Match_M))
  
  tib
}

cmpInfII = function(tib, fieldA, fieldB, mu='D', del='_',
                    verbose=0,vt=4,tc=1,tabsStr='') {
  funcTag <- 'cmpInfI'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  fieldA <- rlang::sym(fieldA)
  fieldB <- rlang::sym(fieldB)
  
  tib <- tib %>% 
    dplyr::mutate(Man_Match=stringr::str_to_upper(!!fieldA)==stringr::str_to_upper(!!fieldB) )
  
  tib
}


# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                       Basic Reverse/Complement Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
revCmp = function(x) {
  Biostrings::reverse(x) %>% cmpl()
}

cmpl = function(x) {
  tr(x, 'ACTGRYSWKMBDHVactgryswkmbdhv','TGACYRSWMKVHDBtgacyrswmkvhdb')
  # x <- tr(x, 'ACTGRYSWKMBDHV','TGACYRSWMKVHDB')
  # x <- tr(x, 'actgryswkmbdhv','tgacyrswmkvhdb')
  # x
}

tr = function(x, old, new) {
  Biostrings::chartr(old, new, x)
}

# mapA = function(x) {
#   if (length(MAP_A[[x]])==0) return(x)
#   MAP_A[[x]]
# }
# INIT_MAP_A = function() {
#   MAP_A <- NULL
#   MAP_A[['AA']] <- 'aa'
#   MAP_A[['Aa']] <- 'aa'
#   MAP_A[['aA']] <- 'aa'
#   MAP_A[['aa']] <- 'aa'
#   
#   MAP_A
# }
# MAP_A <- INIT_MAP_A()

# Data generated with:: /Users/bbarnes/Documents/Projects/scripts/mapMD.pl
INIT_MAP_M = function() {
  MAP <- NULL
  MAP[['CG']] <- 'cG'
  MAP[['CR']] <- 'cR'
  MAP[['CS']] <- 'cS'
  MAP[['CK']] <- 'cK'
  MAP[['CB']] <- 'cB'
  MAP[['CD']] <- 'cD'
  MAP[['CV']] <- 'cV'
  MAP[['YG']] <- 'yG'
  MAP[['YR']] <- 'yR'
  MAP[['YS']] <- 'yS'
  MAP[['YK']] <- 'yK'
  MAP[['YB']] <- 'yB'
  MAP[['YD']] <- 'yD'
  MAP[['YV']] <- 'yV'
  MAP[['SG']] <- 'sG'
  MAP[['SR']] <- 'sR'
  MAP[['SS']] <- 'sS'
  MAP[['SK']] <- 'sK'
  MAP[['SB']] <- 'sB'
  MAP[['SD']] <- 'sD'
  MAP[['SV']] <- 'sV'
  MAP[['MG']] <- 'mG'
  MAP[['MR']] <- 'mR'
  MAP[['MS']] <- 'mS'
  MAP[['MK']] <- 'mK'
  MAP[['MB']] <- 'mB'
  MAP[['MD']] <- 'mD'
  MAP[['MV']] <- 'mV'
  MAP[['BG']] <- 'bG'
  MAP[['BR']] <- 'bR'
  MAP[['BS']] <- 'bS'
  MAP[['BK']] <- 'bK'
  MAP[['BB']] <- 'bB'
  MAP[['BD']] <- 'bD'
  MAP[['BV']] <- 'bV'
  MAP[['HG']] <- 'hG'
  MAP[['HR']] <- 'hR'
  MAP[['HS']] <- 'hS'
  MAP[['HK']] <- 'hK'
  MAP[['HB']] <- 'hB'
  MAP[['HD']] <- 'hD'
  MAP[['HV']] <- 'hV'
  MAP[['VG']] <- 'vG'
  MAP[['VR']] <- 'vR'
  MAP[['VS']] <- 'vS'
  MAP[['VK']] <- 'vK'
  MAP[['VB']] <- 'vB'
  MAP[['VD']] <- 'vD'
  MAP[['VV']] <- 'vV'
  
  MAP
}
MAP_M <- INIT_MAP_M()

INIT_MAP_D = function() {
  MAP <- NULL
  
  MAP[['CG']] <- 'yG'
  MAP[['CR']] <- 'yR'
  MAP[['CS']] <- 'yS'
  MAP[['CK']] <- 'yK'
  MAP[['CB']] <- 'yB'
  MAP[['CD']] <- 'yD'
  MAP[['CV']] <- 'yV'
  MAP[['YG']] <- 'yG'
  MAP[['YR']] <- 'yR'
  MAP[['YS']] <- 'yS'
  MAP[['YK']] <- 'yK'
  MAP[['YB']] <- 'yB'
  MAP[['YD']] <- 'yD'
  MAP[['YV']] <- 'yV'
  MAP[['SG']] <- 'bG'
  MAP[['SR']] <- 'bR'
  MAP[['SS']] <- 'bS'
  MAP[['SK']] <- 'bK'
  MAP[['SB']] <- 'bB'
  MAP[['SD']] <- 'bD'
  MAP[['SV']] <- 'bV'
  MAP[['MG']] <- 'hG'
  MAP[['MR']] <- 'hR'
  MAP[['MS']] <- 'hS'
  MAP[['MK']] <- 'hK'
  MAP[['MB']] <- 'hB'
  MAP[['MD']] <- 'hD'
  MAP[['MV']] <- 'hV'
  MAP[['BG']] <- 'bG'
  MAP[['BR']] <- 'bR'
  MAP[['BS']] <- 'bS'
  MAP[['BK']] <- 'bK'
  MAP[['BB']] <- 'bB'
  MAP[['BD']] <- 'bD'
  MAP[['BV']] <- 'bV'
  MAP[['HG']] <- 'hG'
  MAP[['HR']] <- 'hR'
  MAP[['HS']] <- 'hS'
  MAP[['HK']] <- 'hK'
  MAP[['HB']] <- 'hB'
  MAP[['HD']] <- 'hD'
  MAP[['HV']] <- 'hV'
  MAP[['VG']] <- 'nG'
  MAP[['VR']] <- 'nR'
  MAP[['VS']] <- 'nS'
  MAP[['VK']] <- 'nK'
  MAP[['VB']] <- 'nB'
  MAP[['VD']] <- 'nD'
  MAP[['VV']] <- 'nV'
  
  MAP
}
MAP_D <- INIT_MAP_D()

INIT_QMAP_U = function() {
  MAP <- NULL
  
  # Upper Case::
  MAP[['A']] <- 'A'
  MAP[['C']] <- 'C'  # '-'
  MAP[['G']] <- 'G'  # '-'
  MAP[['T']] <- 'T'
  
  MAP[['R']] <- 'A'
  MAP[['Y']] <- 'T'
  MAP[['S']] <- 'S'  # '-'
  MAP[['W']] <- 'W'  # '-'
  
  MAP[['K']] <- 'T'
  MAP[['M']] <- 'A'
  MAP[['B']] <- 'T'
  MAP[['D']] <- 'W'
  MAP[['H']] <- 'W'
  MAP[['V']] <- 'A'
  
  MAP[['N']] <- 'W'
  
  # Lower Case::
  MAP[['a']] <- 'a'
  MAP[['c']] <- 'c'  # '-'
  MAP[['g']] <- 'g'  # '-'
  MAP[['t']] <- 't'
  
  MAP[['r']] <- 'a'
  MAP[['y']] <- 't'
  MAP[['s']] <- 's'  # '-'
  MAP[['w']] <- 'w'  # '-'
  
  MAP[['k']] <- 't'
  MAP[['m']] <- 'a'
  MAP[['b']] <- 't'
  MAP[['d']] <- 'w'
  MAP[['h']] <- 'w'
  MAP[['v']] <- 'a'
  
  MAP[['n']] <- 'w'
  
  MAP
}
QMAP_U <- INIT_QMAP_U()

INIT_QMAP_M = function() {
  MAP <- NULL
  
  # Upper Case::
  MAP[['A']] <- 'A'  # '-'
  MAP[['C']] <- 'C'
  MAP[['G']] <- 'G'
  MAP[['T']] <- 'T'  # '-'
  
  MAP[['R']] <- 'G'
  MAP[['Y']] <- 'C'
  MAP[['S']] <- 'S'  # '-'
  MAP[['W']] <- 'W'  # '-'
  
  MAP[['K']] <- 'G'
  MAP[['M']] <- 'C'
  MAP[['B']] <- 'S'
  MAP[['D']] <- 'G'
  MAP[['H']] <- 'C'
  MAP[['V']] <- 'S'
  
  MAP[['N']] <- 'S'
  
  # Lower Case::
  MAP[['a']] <- 'a'  # '-'
  MAP[['c']] <- 'c'
  MAP[['g']] <- 'g'
  MAP[['t']] <- 't'  # '-'
  
  MAP[['r']] <- 'g'
  MAP[['y']] <- 'c'
  MAP[['s']] <- 's'  # '-'
  MAP[['w']] <- 'w'  # '-'
  
  MAP[['k']] <- 'g'
  MAP[['m']] <- 'c'
  MAP[['b']] <- 's'
  MAP[['d']] <- 'g'
  MAP[['h']] <- 'c'
  MAP[['v']] <- 's'
  
  MAP[['n']] <- 's'
  
  MAP
}
QMAP_M <- INIT_QMAP_M()

INIT_IUPAC_EQ = function() {
  MAP <- NULL
  
  MAP[['AA']] <- TRUE
  MAP[['AD']] <- TRUE
  MAP[['AH']] <- TRUE
  MAP[['AM']] <- TRUE
  MAP[['AN']] <- TRUE
  MAP[['AR']] <- TRUE
  MAP[['AV']] <- TRUE
  MAP[['AW']] <- TRUE
  MAP[['BB']] <- TRUE
  MAP[['BC']] <- TRUE
  MAP[['BG']] <- TRUE
  MAP[['BK']] <- TRUE
  MAP[['BN']] <- TRUE
  MAP[['BS']] <- TRUE
  MAP[['BT']] <- TRUE
  MAP[['BY']] <- TRUE
  MAP[['CB']] <- TRUE
  MAP[['CC']] <- TRUE
  MAP[['CH']] <- TRUE
  MAP[['CM']] <- TRUE
  MAP[['CN']] <- TRUE
  MAP[['CS']] <- TRUE
  MAP[['CV']] <- TRUE
  MAP[['CY']] <- TRUE
  MAP[['DA']] <- TRUE
  MAP[['DD']] <- TRUE
  MAP[['DG']] <- TRUE
  MAP[['DK']] <- TRUE
  MAP[['DN']] <- TRUE
  MAP[['DR']] <- TRUE
  MAP[['DT']] <- TRUE
  MAP[['DW']] <- TRUE
  MAP[['GB']] <- TRUE
  MAP[['GD']] <- TRUE
  MAP[['GG']] <- TRUE
  MAP[['GK']] <- TRUE
  MAP[['GN']] <- TRUE
  MAP[['GR']] <- TRUE
  MAP[['GS']] <- TRUE
  MAP[['GV']] <- TRUE
  MAP[['HA']] <- TRUE
  MAP[['HC']] <- TRUE
  MAP[['HH']] <- TRUE
  MAP[['HM']] <- TRUE
  MAP[['HN']] <- TRUE
  MAP[['HT']] <- TRUE
  MAP[['HW']] <- TRUE
  MAP[['HY']] <- TRUE
  MAP[['KB']] <- TRUE
  MAP[['KD']] <- TRUE
  MAP[['KG']] <- TRUE
  MAP[['KK']] <- TRUE
  MAP[['KN']] <- TRUE
  MAP[['KT']] <- TRUE
  MAP[['MA']] <- TRUE
  MAP[['MC']] <- TRUE
  MAP[['MH']] <- TRUE
  MAP[['MM']] <- TRUE
  MAP[['MN']] <- TRUE
  MAP[['MV']] <- TRUE
  MAP[['NA']] <- TRUE
  MAP[['NB']] <- TRUE
  MAP[['NC']] <- TRUE
  MAP[['ND']] <- TRUE
  MAP[['NG']] <- TRUE
  MAP[['NH']] <- TRUE
  MAP[['NK']] <- TRUE
  MAP[['NM']] <- TRUE
  MAP[['NN']] <- TRUE
  MAP[['NR']] <- TRUE
  MAP[['NS']] <- TRUE
  MAP[['NT']] <- TRUE
  MAP[['NV']] <- TRUE
  MAP[['NW']] <- TRUE
  MAP[['NY']] <- TRUE
  MAP[['RA']] <- TRUE
  MAP[['RD']] <- TRUE
  MAP[['RG']] <- TRUE
  MAP[['RN']] <- TRUE
  MAP[['RR']] <- TRUE
  MAP[['RV']] <- TRUE
  MAP[['SB']] <- TRUE
  MAP[['SC']] <- TRUE
  MAP[['SG']] <- TRUE
  MAP[['SN']] <- TRUE
  MAP[['SS']] <- TRUE
  MAP[['SV']] <- TRUE
  MAP[['TB']] <- TRUE
  MAP[['TD']] <- TRUE
  MAP[['TH']] <- TRUE
  MAP[['TK']] <- TRUE
  MAP[['TN']] <- TRUE
  MAP[['TT']] <- TRUE
  MAP[['TW']] <- TRUE
  MAP[['TY']] <- TRUE
  MAP[['VA']] <- TRUE
  MAP[['VC']] <- TRUE
  MAP[['VG']] <- TRUE
  MAP[['VM']] <- TRUE
  MAP[['VN']] <- TRUE
  MAP[['VR']] <- TRUE
  MAP[['VS']] <- TRUE
  MAP[['VV']] <- TRUE
  MAP[['WA']] <- TRUE
  MAP[['WD']] <- TRUE
  MAP[['WH']] <- TRUE
  MAP[['WN']] <- TRUE
  MAP[['WT']] <- TRUE
  MAP[['WW']] <- TRUE
  MAP[['YB']] <- TRUE
  MAP[['YC']] <- TRUE
  MAP[['YH']] <- TRUE
  MAP[['YN']] <- TRUE
  MAP[['YT']] <- TRUE
  MAP[['YY']] <- TRUE
  
  MAP
}
IUPAC_EQ <- INIT_IUPAC_EQ()



# End of file
