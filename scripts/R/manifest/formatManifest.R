
rm(list=ls(all=TRUE))
while (dev.cur()>1) dev.off()

suppressWarnings(suppressPackageStartupMessages(require("optparse",quietly=TRUE)))

suppressWarnings(suppressPackageStartupMessages(require("tidyverse")) )
suppressWarnings(suppressPackageStartupMessages(require("stringr")) )
suppressWarnings(suppressPackageStartupMessages(require("glue")) )
suppressWarnings(suppressPackageStartupMessages(require("scales")) )
suppressWarnings(suppressPackageStartupMessages(require("matrixStats")) )

# Sesame Packages:: NOT SURE WE NEEDS THIS CALL
# suppressWarnings(suppressPackageStartupMessages(require("sesame")) )

# Parallel Computing Packages
suppressWarnings(suppressPackageStartupMessages(require("BiocParallel")) )
suppressWarnings(suppressPackageStartupMessages(require("foreach")) )
suppressWarnings(suppressPackageStartupMessages(require("parallel")) )
suppressWarnings(suppressPackageStartupMessages(require("doParallel")) )
suppressWarnings(suppressPackageStartupMessages(require("Rcpp")) )

suppressWarnings(suppressPackageStartupMessages(require("GGally")) )
suppressWarnings(suppressPackageStartupMessages(require("hexbin")) )

opt <- NULL
opt$topDir <- '/Users/bbarnes/Documents/Projects/darkmatter/Projects/workhorse'
opt$manDir <- '/Users/bbarnes/Documents/Projects/manifests/methylation'

epic.B2.csv <- file.path(opt$manDir, 'MethylationEPIC_v-1-0_B2.csv.gz')
epic.B2.tib <- readr::read_csv(epic.B2.csv, skip=7)
ctl.B2.idx <- grep("Control", epic.B2.tib$IlmnID)
epic.B2.tib <- epic.B2.tib %>% head(n=ctl.B2.idx-1)

epic.hg19.man <- sesameData::sesameDataGet('EPIC.hg19.manifest')
epic.hg19.tib <- epic.hg19.man %>% as.data.frame() %>% rownames_to_column('Probe_ID') %>% tibble::as_tibble()


b0.man <- pre.man.tib.add %>% filter(!str_starts(Probe_ID, 'ctl')) %>% arrange(Probe_ID)
b2.man <- epic.B2.tib %>% dplyr::rename(Probe_ID=IlmnID) %>% filter(!str_starts(Probe_ID, 'ctl')) %>% arrange(Probe_ID)
s4.man <- epic.hg19.tib %>% filter(!str_starts(Probe_ID, 'ctl')) %>% arrange(Probe_ID)

bo.cpgs <- b0.man %>% dplyr::select(Probe_ID) %>% dplyr::mutate(SRC='B0')
b2.cpgs <- b2.man %>% dplyr::select(Probe_ID) %>% dplyr::mutate(SRC='B2')
s4.cpgs <- s4.man %>% dplyr::select(Probe_ID) %>% dplyr::mutate(SRC='S4')

cpgs <- bo.cpgs %>% 
  dplyr::full_join(b2.cpgs, by='Probe_ID', suffix=c('.b0', '.b2')) %>%
  dplyr::full_join(s4.cpgs, by='Probe_ID') %>%
  dplyr::rename(SRC.s4=SRC)

cpgs %>% summarise(B0=count(is.na(SRC.b0)), B2=count(is.na(SRC.b2)), S4=count(is.na(SRC.s4)))

sesameManifest = function(platform='EPIC', build='hg19', verbose=0, vt=1) {
  funcTag <- 'sesameManifest'
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:\tStarting..."),"\n", sep='')
  name <- paste(platform,build,'manifest', sep='.')
  tib <- sesameData::sesameDataGet(name) %>% 
    as.data.frame() %>% 
    rownames_to_column('Probe_ID') %>% 
    tibble::as_tibble() %>% 
    dplyr::mutate(col=NA, Masked=NA) %>%
    dplyr::mutate(col=str_sub(col, 1,1),
                  Masked=case_when(
                    MASK_mapping==TRUE ~ 'Mapping',
                    MASK_general==TRUE ~ 'General'
                  )) %>%
    dplyr::mutate(M=address_B,
                  U=address_A,
                  COLOR_CHANNEL=channel,
                  Probe_Type=probeType,
                  Probe_Source='EPIC-B4',
                  Next_Base=nextBaseRef,
                  Probe_CpG_Cnt=probeCpGcnt) %>%
    dplyr::select(Probe_ID, M, U, COLOR_CHANNEL, col, Probe_Type, Probe_Source, Next_Base, 
                  seqnames, start, end, width, strand, Probe_CpG_Cnt, Masked) %>%
    arrange(Probe_ID)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:\tDone."),"\n", sep='')

  tib
}
ses_man.tib <- sesameManifest(verbose=2)

# Probe_ID M U DESIGN COLOR_CHANNEL col Probe_Type Probe_Source Next_Base

tmp %>% dplyr::select(Probe_ID, address_B, address_A, channel, designType, nextBaseRef, probeType, orientation, everything()) %>% head

tmp %>% head %>%
  dplyr::mutate(col=NA, Masked=NA) %>%
  dplyr::mutate(col=str_sub(col, 1,1),
                Masked=case_when(
                  MASK_mapping==TRUE ~ 'Mapping',
                  MASK_general==TRUE ~ 'General'
                )) %>%
  dplyr::mutate(M=address_B,
                U=address_A,
                COLOR_CHANNEL=channel,
                Probe_Type=probeType,
                Probe_Source='EPIC-B4',
                Next_Base=nextBaseRef,
                Probe_CpG_Cnt=probeCpGcnt) %>%
  dplyr::select(seqnames, start, end, width, strand,
                Probe_ID, M, U, COLOR_CHANNEL, col, Probe_Type, Probe_Source, Next_Base, Probe_CpG_Cnt, Masked)


# End of file
