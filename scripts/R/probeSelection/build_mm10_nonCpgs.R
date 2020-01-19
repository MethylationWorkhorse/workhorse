
rm(list=ls(all=TRUE))

suppressWarnings(suppressPackageStartupMessages(require("optparse",quietly=TRUE)))

suppressWarnings(suppressPackageStartupMessages(require("tidyverse")) )
suppressWarnings(suppressPackageStartupMessages(require("stringr")) )
suppressWarnings(suppressPackageStartupMessages(require("glue")) )
suppressWarnings(suppressPackageStartupMessages(require("scales")) )
suppressWarnings(suppressPackageStartupMessages(require("matrixStats")) )
suppressWarnings(suppressPackageStartupMessages(require("grid")) )

suppressWarnings(suppressPackageStartupMessages(require("Biostrings")) )

# Parallel Computing Packages
suppressWarnings(suppressPackageStartupMessages(require("doParallel")) )

doParallel::registerDoParallel()
num_cores   <- detectCores()
num_workers <- getDoParWorkers()

COM <- ","
TAB <- "\t"
RET <- "\n"
BNG <- "|"

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                             Default Params::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

opt <- NULL
opt$prgmTag <- 'improbe_design'

opt$single   <- FALSE
opt$parallel <- FALSE

opt$desDir <- NULL
opt$manDir <- NULL
opt$outDir <- NULL

opt$name   <- NULL

opt$prbSrc <- NULL
opt$format <- NULL
opt$subCnt <- NULL
opt$blockSize <- 10000

opt$exactMatch <- FALSE
opt$impExecute <- FALSE

opt$verbosity <- 3

improbeHeader <- c('Seq_ID', 'Sequence', 'Genome_Build', 'Chromosome', 'Coordinate', 'CpG_Island')

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                       System (Linux/Mac) Directories::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
opt$mac_dir <- '/Users/bbarnes/Documents/Projects'
opt$lix_dir <- '/illumina/scratch/darkmatter/data'
opt$sys_key <- NULL

if (dir.exists(opt$mac_dir)) {
  opt$sysName <- 'mac'
  opt$prgmTag  <- paste(opt$prgmTag, 'RStudio', sep='_')
  opt$prgDir   <- opt$mac_dir
  opt$manDir   <- '/Users/bbarnes/Documents/Projects/manifests/methylation'
  opt$desDir   <- '/Users/bbarnes/Documents/Projects/methylation/NZT_Limitless/data/imDesignOutput'
  opt$topDir   <- file.path(opt$mac_dir, 'workhorse')
  opt$gitDir   <- file.path(opt$topDir, 'git/workhorse')
  
  # NOT REAL ON MAC::
  opt$improbe  <- file.path(opt$lix_dir, 'bin/improbe')
  opt$imp13mer <- file.path(opt$lix_dir, 'dat/improbe/human-36.1-methyl-and-unmethyl-13mer-s3-for-infinium-methylation.dat')
  opt$impTango <- file.path(opt$lix_dir, 'dat/improbe/Tango_A_or_B_11mer_s1.dat')
} else if (dir.exists(opt$lix_dir)) {
  opt$sysName <- 'linux'
  opt$prgDir   <- opt$lix_dir
  opt$manDir   <- '/illumina/scratch/darkmatter/manifests/data/methylation'
  opt$desDir   <- NULL
  opt$topDir   <- file.path(opt$lix_dir)
  opt$gitDir   <- file.path(opt$topDir, 'git/workhorse')
  
  opt$improbe  <- file.path(opt$lix_dir, 'bin/improbe')
  opt$imp13mer <- file.path(opt$lix_dir, 'dat/improbe/human-36.1-methyl-and-unmethyl-13mer-s3-for-infinium-methylation.dat')
  opt$impTango <- file.path(opt$lix_dir, 'dat/improbe/Tango_A_or_B_11mer_s1.dat')
} else {
  stop(glue::glue("[ERROR]: Unrecognized prgDir=({opt$mac_dir}, {opt$lix_dir})"),"\n", sep='')
  q()
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                           Command Line Options::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
args.dat <- commandArgs(trailingOnly = FALSE)
opt_parser <- NULL
if (args.dat[1]=='RStudio') {
  opt$srcDir   <- file.path(opt$gitDir, 'scripts')
  opt$prgmPath <- file.path(opt$srcDir, 'R/plotsByExperiment.R')
  
  opt$format <- 'des'
  opt$format <- 'man'
  
  opt$prbSrc <- 'cg'
  opt$prbSrc <- 'rs'
  opt$prbSrc <- 'ch'
  
  opt$subCnt <- 20000
  
  opt$name   <- 'EPIC-B2'
  
  opt$exactMatch <- FALSE
  opt$impExecute <- FALSE
  
  opt$outDir <- file.path(opt$topDir, 'workspace/improbe-non-cpg')
  
} else {
  opt$prgmPath <- substring(args.dat[grep("--file=", args.dat)], 8)
  opt$runMode  <- 'CommandLine'
  opt$prgmTag  <- sub('\\.R$', '', basename(substring(args.dat[grep("--file=", args.dat)], 8)))
  # opt$srcDir   <- dirname(normalizePath(dirname(substring(args.dat[grep("--file=", args.dat)], 8)) ))
  opt$srcDir   <- dirname(dirname(normalizePath(dirname(substring(args.dat[grep("--file=", args.dat)], 8)) )))
  opt$gitDir   <- file.path(opt$srcDir, '../../..')
  
  args.dat <- commandArgs(trailingOnly = TRUE)
  option_list = list(
    make_option(c("--Rscript"), type="character", default=opt$Rscript, 
                help="Rscript path [default= %default]", metavar="character"),
    make_option(c("--prgmTag"), type="character", default=opt$prgmTag, 
                help="Program name [default= %default]", metavar="character"),
    make_option(c("--sysName"), type="character", default=opt$sysName, 
                help="System Name. Should be auto-detected [linux, mac] [default= %default]", metavar="character"),
    
    # Directories/Files::
    make_option(c("--topDir"), type="character", default=opt$topDir, 
                help="Top data directory location [default= %default]", metavar="character"),
    make_option(c("--srcDir"), type="character", default=opt$srcDir, 
                help="Source data directory location [default= %default]", metavar="character"),
    make_option(c("--manDir"), type="character", default=opt$manDir, 
                help="Manifest directory location [default= %default]", metavar="character"),
    make_option(c("--desDir"), type="character", default=opt$desDir, 
                help="Improbe Design data directory location [default= %default]", metavar="character"),
    make_option(c("--outDir"), type="character", default=opt$outDir, 
                help="Output directory location [default= %default]", metavar="character"),
    make_option(c("--name"), type="character", default=opt$name, 
                help="Output name [default= %default]", metavar="character"),
    
    # Processing Parameters::
    make_option(c("--single"), action="store_true", default=opt$single, 
                help="Boolean variable to run a single sample on a single-core [default= %default]", metavar="boolean"),
    make_option(c("--parallel"), action="store_true", default=opt$parallel, 
                help="Boolean variable to run parallel on multi-core [default= %default]", metavar="boolean"),
    
    # Run Parameters::
    make_option(c("--prbSrc"), type="character", default=opt$prbSrc,
                help="Probe Type (cg,ch,rs) [default= %default]", metavar="character"),
    make_option(c("--format"), type="character", default=opt$format,
                help="Input file format  [default= %default]", metavar="character"),
    make_option(c("--subCnt"), type="double", default=opt$subCnt, 
                help="Sub Sample Input Count [default= %default]", metavar="double"),
    make_option(c("--blockSize"), type="double", default=opt$blockSize, 
                help="Loci (x4) block size (mostly for cluster and memory partitioning) [default= %default]", metavar="double"),
    make_option(c("--exactMatch"), action="store_true", default=opt$exactMatch, 
                help="Boolean variable for exact probe matching vs. mismatch counting [default= %default]", metavar="boolean"),
    
    # improbe Parameters::
    make_option(c("--impExecute"), action="store_true", default=opt$impExecute, 
                help="Boolean variable for executing improbe [default= %default]", metavar="boolean"),
    
    make_option(c("--improbe"), type="character", default=opt$improbe,
                help="improbe executable [default= %default]", metavar="character"),
    make_option(c("--impTango"), type="character", default=opt$impTango,
                help="improbe tangos file [default= %default]", metavar="character"),
    make_option(c("--imp13mer"), type="character", default=opt$imp13mer,
                help="improbe 13mer file [default= %default]", metavar="character"),
    
    # Verbosity::
    make_option(c("-v", "--verbosity"), type="double", default=opt$verbosity, 
                help="0-5 (5 is very verbosity) [default= %default]", metavar="double")
  )
  opt_parser = OptionParser(option_list=option_list)
  opt = parse_args(opt_parser)
}

if (is.null(opt$prgmTag)  || 
    is.null(opt$prbSrc)   ||
    is.null(opt$format)   ||
    is.null(opt$topDir)   ||
    is.null(opt$manDir)   ||
    is.null(opt$outDir)   ||
    is.null(opt$name)     ||
    is.null(opt$srcDir) ) {
  
  if (!is.null(opt_parser)) print_help(opt_parser)
  print(as.data.frame(tidyr::gather(dplyr::bind_rows(opt),"Options", "Value")))
  stop("Null arguments!\n\n")
  q()
}
if (is.null(opt$format) && is.null(opt$desDir)) {
  if (!is.null(opt_parser)) print_help(opt_parser)
  print(as.data.frame(tidyr::gather(dplyr::bind_rows(opt),"Options", "Value")))
  stop("If using improbe design files for analysis you must provide improbe design directory (--desDir)!\n\n")
  q()
}

# Automatic Params
opt$MinMisMatch <- 0
if (opt$prbSrc=='ch') opt$MinMisMatch <- 1

# R Source Code Files
opt$fIOSrc   <- file.path(opt$srcDir, 'R/workhorse_IO_functions.R')
opt$wrkSrc   <- file.path(opt$srcDir, 'R/workhorse_functions.R')
opt$plotSrc  <- file.path(opt$srcDir, 'R/workhorse_plotting_functions.R')
opt$impSrc   <- file.path(opt$srcDir, 'R/probeSelection/improbe_functions.R')

if (opt$verbosity>0) print(as.data.frame(tidyr::gather(dplyr::bind_rows(opt),"Options", "Value")))

stopifnot(file.exists(opt$fIOSrc))
stopifnot(file.exists(opt$plotSrc))
stopifnot(file.exists(opt$impSrc))
# stopifnot(file.exists(opt$wrkSrc))

source(opt$fIOSrc)
source(opt$plotSrc)
source(opt$impSrc)
# source(opt$wrkSrc)

# Output Directory::
if (!dir.exists(opt$outDir)) dir.create(opt$outDir, recursive=TRUE)
opt$out_man_csv <- file.path(opt$outDir, paste(opt$name,'_',opt$prbSrc,'.manifest.csv.gz', sep=''))

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                             Default Files::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
mm10_des_tsv <- file.path(opt$topDir, 'designs/merged_with_raw_ordered.with-header.tsv.gz')
mm10_non_csv <- file.path(opt$topDir, 'designs/merged_with_raw_ordered.with-header.non-cpg.csv.gz')

clean_mm10 <- FALSE
if (clean_mm10) {
  mm10_des_tib <- suppressMessages(suppressWarnings(readr::read_tsv(mm10_des_tsv) ))
  mm10_des_tib <- mm10_des_tib %>% 
    dplyr::mutate(
      Seq_ID=stringr::str_remove_all(Seq_ID, ' '),
      Probe_Type=case_when(
        stringr::str_starts(Seq_ID, 'cg')            ~ 'cg',
        stringr::str_starts(User_Transcript, 'SNP')  ~ 'rs',
        stringr::str_starts(User_Transcript, 'CpH')  ~ 'ch',
        stringr::str_starts(User_Transcript, 'RMSK') ~ 'rp',
        
        TRUE ~ NA_character_),
      diNUC=case_when(
        Probe_Type=='cg' ~ 'CG',
        Probe_Type=='rp' ~ 'CG',
        Probe_Type=='ch' ~ stringr::str_remove(User_Transcript, '^CpH;'),
        Probe_Type=='rs' ~ stringr::str_remove_all(stringr::str_replace(User_Transcript, '^.*;([A-Z]+)>([A-Z]+);.*$', '\\$1\\$2'), '\\\\'),
        TRUE ~ NA_character_
      )
    ) %>% dplyr::select(Probe_Type, diNUC, everything())
  
  mm10_non_tib <- mm10_des_tib %>% dplyr::filter(Probe_Type!='cg')
  readr::write_csv(mm10_non_tib, mm10_non_csv)
}
mm10_non_tib <- suppressMessages(suppressWarnings(readr::read_csv(mm10_non_csv) ))

# Extract improbe design format:: imp_des_tib
imp_mm10_tib <- mm10_non_tib %>% dplyr::select(Seq_ID:UnMethyl_Next_Base_Score)
imp_des_tib  <- loadImprobeDesign(src_des_tib=imp_mm10_tib) %>% dplyr::mutate(
  Chromosome=as.character(Chromosome),
  Genome_Build=case_when(
    is.na(Genome_Build) ~ 'mm10',
    Genome_Build=='.'   ~ 'mm10',
    TRUE ~ Genome_Build
  )
)

# Format all IUPAC Design Fields:: src_man_tib
src_man_tib <- mm10_non_tib %>% 
  dplyr::select(Seq_ID, Forward_Sequence, Genome_Build, Chromosome, Coordinate, Probe_Type, diNUC,User_Infinium) %>%
  dplyr::mutate(
    Genome_Build=case_when(
    is.na(Genome_Build) ~ 'mm10',
    Genome_Build=='.'   ~ 'mm10',
    TRUE ~ Genome_Build),
    CpG_Island='FALSE',
    PRB_CGN=Seq_ID,
    IUPAC_Forward_Sequence=Forward_Sequence,
    Strand_FR=case_when(
      Probe_Type=='rs' ~ stringr::str_remove_all(stringr::str_replace(Seq_ID, '^.*_([RF])$', '\\$1'), '\\\\'),
      TRUE ~ 'F'
    ),
    diSNP=case_when(
      Probe_Type=='rs' ~ stringr::str_remove_all(stringr::str_replace(Seq_ID, '^.*_([A-Z]+)_[RF]$', '\\$1'), '\\\\'),
      Probe_Type=='ch' ~ paste0('R',stringr::str_sub(diNUC,2,2)),
      TRUE ~ diNUC
    ),
    diSNP=case_when(
      Probe_Type=='rs' & Strand_FR=='R' ~ revCmp(diSNP),
      # Probe_Type=='rs' & Strand_FR=='F' ~ revCmp(diSNP),
      TRUE ~ diSNP
    ),
    Infinium_Design_Type=User_Infinium,
    IUPAC_Forward_Sequence=stringr::str_replace(IUPAC_Forward_Sequence, '\\[[A-Za-z][A-Za-z]\\]', paste0('[',diSNP,']'))
  ) %>% dplyr::rename(Sequence=Forward_Sequence, PRB_DES=Probe_Type) %>%
  dplyr::select(Seq_ID,Sequence,Genome_Build,Chromosome,Coordinate,CpG_Island,
                PRB_CGN,PRB_DES,diNUC,IUPAC_Forward_Sequence,Strand_FR,diSNP,Infinium_Design_Type)

# src_man_tib %>% dplyr::filter(PRB_DES=='ch') %>% dplyr::select(Seq_ID, diNUC, diSNP, Strand_FR, IUPAC_Forward_Sequence)
# src_man_tib %>% dplyr::filter(PRB_DES=='rp') %>% dplyr::select(Seq_ID, diNUC, diSNP, Strand_FR, IUPAC_Forward_Sequence)
# src_man_tib %>% dplyr::filter(PRB_DES=='rs') %>% dplyr::select(Seq_ID, diNUC, diSNP, Strand_FR, IUPAC_Forward_Sequence)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                    Traditional EPIC Code for Validation::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
if (FALSE) {
  # improbe Design Files::
  cpg_des_tsv <- file.path(opt$desDir, 'hs_ref_both_both_RefSeq37_probeDesignOutput.tsv.gz')
  cpg_des_tsv <- file.path(opt$desDir, 'hs_ref_both_both_RefSeq37_probeDesignOutput.1k.tsv.gz')
  cph_des_tsv <- file.path(opt$desDir, 'ch.improbe-designOutput.tsv.gz')
  
  # 450k Formats::
  # cph_man_csv <- file.path(opt$manDir, 'ch-repair/HumanMethylation450_15017482_v.1.2.ch-noAnnotation.csv.gz')
  # snp_man_csv <- file.path(opt$manDir, 'rs-repair/HumanMethylation450_15017482_v.1.2.rs.no-FwdSeq.csv.gz')
  # snp_swp_tsv <- file.path(opt$manDir, 'rs-repair/hs_ref_RefSeq37_probeDesign_SWAP.tsv.gz')
  
  # EPIC Formats::
  # cpg_man_csv <- file.path(opt$manDir, 'MethylationEPIC_v-1-0_B4.core.cpg-only.table.1k.csv.gz')
  # cpg_man_csv <- file.path(opt$manDir, 'MethylationEPIC_v-1-0_B4.core.cpg-only.table.csv.gz')
  cpg_man_csv <- file.path(opt$manDir, 'MethylationEPIC_v-1-0_B2.csv.gz')
  cpg_man_csv <- file.path(opt$manDir, 'MethylationEPIC_v-1-0_B2.1k.csv.gz')
  
  cph_man_csv <- file.path(opt$manDir, 'ch-repair/MethylationEPIC_v-1-0_B4.ch-noAnnotation.csv.gz')
  snp_man_csv <- file.path(opt$manDir, 'rs-repair/MethylationEPIC_v-1-0_B4.rs.no-header.csv.gz')
  snp_swp_tsv <- file.path(opt$manDir, 'rs-repair/rs_swap_improbe_input.tsv.gz')
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                         Load/Properly-Build SNP Manifest::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  # if (opt$format=='des' && opt$prbSrc=='cg') src_man_tib <- loadImprobeDesign(cpg_des_tsv)
  # if (opt$format=='man' && opt$prbSrc=='cg') cpg_man_tib <- loadManifestCG(cpg_man_csv)
  # if (opt$format=='man' && opt$prbSrc=='ch') cph_man_tib <- loadManifestCH(cph_man_csv, ry='R')
  # if (opt$format=='man' && opt$prbSrc=='rs') snp_man_tib <- loadManifestRS(snp_man_csv, snp_swp_tsv)
  
  # MAN - Manifest
  # DES - Template Sequence
  # PRB - Probes
  # IMP - Improbe
  
  cpg_man_tib <- loadManifestCG(cpg_man_csv) %>% dplyr::mutate(PRB_DES='cg', CHR=as.character(CHR))
  cph_man_tib <- loadManifestCH(cph_man_csv, ry='R') %>% dplyr::mutate(PRB_DES='ch', CHR=as.character(CHR))
  snp_man_tib <- loadManifestRS(snp_man_csv, snp_swp_tsv) %>% dplyr::mutate(PRB_DES='rs', CHR=as.character(CHR))
  src_man_tib <- dplyr::bind_rows(cpg_man_tib, cph_man_tib, snp_man_tib)
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                                 Sub Sample::
  if (!is.null(opt$subCnt) && opt$subCnt!=0)
    src_man_tib <- src_man_tib %>% head(n=opt$subCnt)
  src_man_len <- src_man_tib %>% base::nrow()
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                            Format Manifest::
  src_man_tib <- src_man_tib %>% select(IlmnID,Forward_Sequence,Genome_Build,CHR,MAPINFO,PRB_DES,diNUC) %>%
    dplyr::rename(Seq_ID=IlmnID,Chromosome=CHR,Coordinate=MAPINFO,IUPAC_Forward_Sequence=Forward_Sequence) %>%
    dplyr::mutate(Sequence=stringr::str_replace(IUPAC_Forward_Sequence, '\\[[A-Za-z][A-Za-z]\\]', '[CG]'),
                  CpG_Island='FALSE',PRB_CGN=Seq_ID) %>%
    dplyr::select(improbeHeader, PRB_CGN,PRB_DES,diNUC,IUPAC_Forward_Sequence)
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                    Include Full Design Data from improbe::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  imp_inp_tib <- src_man_tib %>% dplyr::select(improbeHeader)
  imp_out_key <- paste(opt$name,opt$prbSrc,paste0(1,'-',src_man_len), sep='_')
  imp_des_tsv <- writeImprobeInput(tib=imp_inp_tib,name=imp_out_key,dir=opt$outDir,run=opt$impExecute,
                                   exe=opt$improbe,impTango=opt$impTango,imp13mer=opt$imp13mer,
                                   verbose=opt$verbosity)
  
  if (opt$sysName=='mac') {
    imp_cpg_tsv <- file.path(opt$topDir, 'improbe/fromCluster/EPIC-B2_cg_1-1000.improbe-output.tsv.gz')
    imp_cph_tsv <- file.path(opt$topDir, 'improbe/fromCluster/EPIC-B2_ch_1-2932.improbe-output.tsv.gz')
    imp_snp_tsv <- file.path(opt$topDir, 'improbe/fromCluster/EPIC-B2_rs_1-58.improbe-output.tsv.gz')
    
    imp_cpg_tib <- loadImprobeDesign(file=imp_cpg_tsv) %>% dplyr::mutate(Chromosome=as.character(Chromosome))
    imp_cph_tib <- loadImprobeDesign(file=imp_cph_tsv) %>% dplyr::mutate(Chromosome=as.character(Chromosome))
    imp_snp_tib <- loadImprobeDesign(file=imp_snp_tsv) %>% dplyr::mutate(Chromosome=as.character(Chromosome))
    imp_des_tib <- dplyr::bind_rows(imp_cpg_tib,imp_cph_tib,imp_snp_tib) # %>%
    #  dplyr::select(-c(Genome_Build,Chromosome,Coordinate))
    
  }
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                      Calcluate Probes on All Strands::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
des_seq_F_C <- src_man_tib %>% 
  dplyr::mutate(FR=TRUE, CO=TRUE, DesSeqN=shearBrac(IUPAC_Forward_Sequence))

des_seq_R_C <- des_seq_F_C %>% dplyr::mutate(
  FR=!FR,CO=CO, DesSeqN=revCmp(DesSeqN) )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                        Build All Design Strands::
bsc_tibs <- NULL

# BSC-Forward-Converted::
bsc_tibs$FC <- des_seq_F_C %>% dplyr::mutate(
  DesBscU = bscUs(DesSeqN),
  DesBscM = bscMs(DesSeqN),
  DesBscD = bscDs(DesSeqN) ) %>%
  dplyr::mutate(Seq_ID)
# BSC-Foward-Opposite::
bsc_tibs$FO <- bsc_tibs$FC %>% dplyr::mutate(
  FR=FR,CO=!CO, 
  # DesSeqN=revCmp(DesSeqN),
  DesBscU=revCmp(DesBscU),
  DesBscM=revCmp(DesBscM),
  DesBscD=revCmp(DesBscD) ) %>%
  dplyr::mutate(Seq_ID)

# BSC-Reverse-Converted::
bsc_tibs$RC <- des_seq_R_C %>% dplyr::mutate(
  DesBscU = bscUs(DesSeqN),
  DesBscM = bscMs(DesSeqN),
  DesBscD = bscDs(DesSeqN) ) %>%
  dplyr::mutate(Seq_ID)
# BSC-Reverse-Opposite::
bsc_tibs$RO <- bsc_tibs$RC %>% dplyr::mutate(
  FR=FR,CO=!CO, 
  # DesSeqN=revCmp(DesSeqN),
  DesBscU=revCmp(DesBscU),
  DesBscM=revCmp(DesBscM),
  DesBscD=revCmp(DesBscD) ) %>%
  dplyr::mutate(Seq_ID)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                  Compare Manifest with Calculated Probes::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
prb_tibs <- NULL

# TEST CODE::
# for (srd in names(bsc_tibs)) { lapply(split(bsc_tibs[[srd]], bsc_tibs[[srd]]$PRB_DES), desAllPrbs) %>% dplyr::bind_rows() }

# Original Code::
prb_tibs <- foreach (srd=names(bsc_tibs), .combine=rbind) %dopar% {
  lapply(split(bsc_tibs[[srd]], bsc_tibs[[srd]]$PRB_DES), desAllPrbs) %>% dplyr::bind_rows()
}

# TBD:: Remove probes that just can't work::
#  - prb_tibs %>% dplyr::filter(PRB1_U!=PRB1_M)
prb_tibs <- prb_tibs %>% dplyr::mutate(
  Bad_Design=case_when(
    stringr::str_to_upper(PRB1_U)==stringr::str_to_upper(PRB1_M) ~ TRUE, 
    TRUE ~ FALSE),
  Valid_Extension=case_when(
    # stringr::str_to_upper(stringr::str_sub(PRB1_U_IUP,50,50))==stringr::str_to_upper(stringr::str_sub(PRB1_M_IUP,50,50)) ~ FALSE,
    stringr::str_to_upper(stringr::str_sub(PRB1_U,50,50))==stringr::str_to_upper(stringr::str_sub(PRB1_M,50,50)) ~ FALSE,
    TRUE ~ TRUE
  )
)

# Original inner join::
# prb_des_tibs <- prb_tibs %>% 
#   dplyr::inner_join(imp_des_tib, 
#                     by=c('Seq_ID','FR','CO','Genome_Build','Chromosome','Coordinate'),
#                     suffix=c('_IUP','_IMP'))
prb_des_tibs <- prb_tibs %>% 
  dplyr::left_join(imp_des_tib, 
                    by=c('Seq_ID','FR','CO','Genome_Build','Chromosome','Coordinate'),
                    suffix=c('_IUP','_IMP'))

# Now we have comparisons::
prb_mat_tibs <- prb_des_tibs %>% 
  cmpInfI_MisMatch('PRB1_U_IUP', 'PRB1_U_IMP', 'PRB1_M_IUP', 'PRB1_M_IMP',
                   verbose=opt$verbosity) # %>% dplyr::select(Man_MisMatch,Man_TarMatch)

select_plotting <- FALSE
if (select_plotting) {
  plot_ord_tib <- prb_mat_tibs %>% dplyr::arrange(-Man_MisMatch) %>% 
    dplyr::group_by(FR,TB,CO,PRB_DES,Man_MisMatch,Man_TarMatch) %>% 
    dplyr::top_n(1,wt=Seq_ID) %>% 
    dplyr::select(Seq_ID,FR,TB,CO,PRB_DES,Man_MisMatch,Man_TarMatch)
  plot_nrows <- plot_ord_tib %>% base::nrow()
} else {
  plot_ord_tib <- prb_mat_tibs %>% dplyr::arrange(-Man_MisMatch) %>% dplyr::distinct(Seq_ID, .keep_all=TRUE)
  plot_nrows <- plot_ord_tib %>% base::nrow()
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                            Generate Plots::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
plot_txt <- file.path(opt$outDir, 'mm10_LEGX_nonCpG_probes.plot.txt')
unlink(plot_txt)

for (ii in seq(1,plot_nrows)) {
  tag_tib <- plot_ord_tib[ii,]
  print(tag_tib)
  
  cur_tib <- prb_mat_tibs %>% dplyr::filter(Seq_ID==tag_tib$Seq_ID[1])
  tibs <- cur_tib %>% srdsToBrac()
  
  # Build Printable Strings for reach strand::
  IlmnId <- cur_tib$Seq_ID
  cat(glue::glue("ii={ii}; IlmnID={IlmnId}:{RET}"))
  strs <- NULL
  strs$F <- prbsToStr(tibs$F, pr=tag_tib$PRB_DES[1], verbose=opt$verbosity)
  strs$R <- prbsToStr(tibs$R, pr=tag_tib$PRB_DES[1], verbose=opt$verbosity)
  
  readr::write_file(paste(unique(cur_tib$Seq_ID), collapse='\t'), plot_txt, append=TRUE)
  readr::write_file("\n", plot_txt, append=TRUE)
  readr::write_file(str_c(strs$F), plot_txt, append=TRUE)
  readr::write_file(str_c(strs$R), plot_txt, append=TRUE)
  
  # cur_tib <- mat_tib %>% dplyr::filter(IlmnID==IlmnId)
  if (FALSE && cur_tib$Infinium_Design_Type[1]=='II') {
    cur_tib %>% dplyr::select(FR,CO, PRB2_D_IUP, PRB2_D_IMP,Man_MisMatch,Man_TarMatch,Bad_Design) %>% print()
  } else {
    prb1U_iup <- cur_tib %>% dplyr::select(FR,CO, PRB1_U_IUP, PRB1_U_IMP,Man_MisMatch,Man_TarMatch,Bad_Design)
    prb1M_iup <- cur_tib %>% dplyr::select(FR,CO, PRB1_M_IUP, PRB1_M_IMP,Man_MisMatch,Man_TarMatch,Bad_Design)
    
    # readr::write_file(str_c(prb1U_iup), plot_txt, append=TRUE)
    # readr::write_file(str_c(prb1M_iup), plot_txt, append=TRUE)
  }
  
  break
}

# Selection::
#  1. Pick all !Bad_Design with !is.na(PRB1_U_IMP)
top_pick_tib  <- prb_mat_tibs %>% dplyr::filter(Bad_Design==FALSE) %>% dplyr::filter(!is.na(PRB1_U_IMP))
top_pick_cgns <- top_pick_tib %>% dplyr::select(Seq_ID) %>% distinct(Seq_ID)
#  2. Pick all !Bad_Design with is.na(PRB1_U_IMP)
mis_pick_cgns <- prb_mat_tibs %>% dplyr::anti_join(top_pick_cgns, by=c('Seq_ID') ) %>% distinct(Seq_ID)

sec_pick_tib <- prb_mat_tibs %>% dplyr::filter(CO==TRUE & Bad_Design==FALSE & Valid_Extension==TRUE) %>%
  dplyr::anti_join(top_pick_cgns, by=c('Seq_ID') ) # %>% 
  # group_by(Seq_ID) %>% summarise(CNT=n()) %>% dplyr::arrange(-CNT)

fin_pick_tib <- top_pick_tib %>% dplyr::bind_rows(sec_pick_tib) %>% dplyr::arrange(Seq_ID)

# Assay_Design_Id, AlleleA_Probe_Id, AlleleA_Probe_Sequence, AlleleB_Probe_Id, AlleleB_Probe_Sequence, Normalization_Bin
fin_ord_tib <- fin_pick_tib %>% dplyr::mutate(
  FR_Char=case_when(FR ~ 'F', TRUE ~ 'R'),
  TB_Char=case_when(TB ~ 'T', TRUE ~ 'B'),
  CO_Char=case_when(CO ~ 'C', TRUE ~ 'O'),
  ID_Char=case_when(
    PRB_DES=='rs' & !stringr::str_starts(Seq_ID,'rs') ~ paste0('rs',Seq_ID),
    TRUE ~ Seq_ID
  ),
  Assay_Design_Id=paste(ID_Char,FR_Char,TB_Char,CO_Char,Infinium_Design_Type, sep='_'),
  AlleleA_Probe_Id=paste(Assay_Design_Id,'A',sep='_'),
  AlleleA_Probe_Sequence=case_when(
    Infinium_Design_Type=='II' ~ PRB2_D,
    Infinium_Design_Type=='I'  ~ PRB1_U_IUP,
    TRUE ~ NA_character_
  ),

  AlleleB_Probe_Id=case_when(
    Infinium_Design_Type=='I' ~ paste(Assay_Design_Id,'B',sep='_'),
    TRUE ~ ''
  ),
  AlleleB_Probe_Sequence=case_when(
    Infinium_Design_Type=='II' ~ '',
    Infinium_Design_Type=='I'  ~ PRB1_M_IUP,
    TRUE ~ NA_character_
  ),
  Normalization_Bin=case_when(
    Infinium_Design_Type=='II' ~ 'C',
    Infinium_Design_Type=='I' & 
      (stringr::str_to_upper(NXB_U)=='A' | stringr::str_to_upper(NXB_M)=='A' | 
         stringr::str_to_upper(NXB_U)=='T' | stringr::str_to_upper(NXB_M)=='T') ~ 'A',
    Infinium_Design_Type=='I' & 
      (stringr::str_to_upper(NXB_U)=='C' | stringr::str_to_upper(NXB_M)=='C' | 
         stringr::str_to_upper(NXB_U)=='G' | stringr::str_to_upper(NXB_M)=='G') ~ 'B',
    TRUE ~ NA_character_
  )
) %>% dplyr::select(Assay_Design_Id, 
                    AlleleA_Probe_Id, AlleleA_Probe_Sequence, 
                    AlleleB_Probe_Id, AlleleB_Probe_Sequence, 
                    Normalization_Bin, Seq_ID) %>%
  dplyr::arrange(Assay_Design_Id)

fin_csv <- file.path(opt$outDir, 'mm10_LEGX_nonCpG_probes.order.csv.gz')
readr::write_csv(fin_ord_tib %>% dplyr::select(-Seq_ID), fin_csv)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                      Solid Cross Method Validate rp's::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
validate_rps <- FALSE
if (validate_rps) {
  rp_csv <- '/Users/bbarnes/Documents/Projects/methylation/LifeEpigentics/to-order/designW.Final_Jan-15-2020.rp-only/Mus_musculus.order_BP1.csv'
  rp_tib <- readr::read_csv(rp_csv, skip=8)
  
  # TBD:: Validate Bin as well!!!
  rp_tib %>% tidyr::separate(Assay_Design_Id, into=c('Seq_ID', 'Strand_FR', 'Strand_TB', 'Strand_CO', 'Infinium'), sep='_', remove=FALSE) %>%
    dplyr::left_join(fin_pick_tib, by='Seq_ID')
  
  rp_tib %>% tidyr::separate(Assay_Design_Id, into=c('Seq_ID', 'Strand_FR', 'Strand_TB', 'Strand_CO', 'Infinium'), sep='_', remove=FALSE) %>%
    dplyr::left_join(fin_pick_tib, by='Seq_ID') %>% 
    dplyr::mutate(
      SeqA_Match=case_when(
        Infinium=='II' & stringr::str_to_upper(AlleleA_Probe_Sequence)==stringr::str_to_upper(PRB2_D) ~ TRUE,
        Infinium=='I' & 
          stringr::str_to_upper(AlleleA_Probe_Sequence)==stringr::str_to_upper(PRB1_U_IUP) &
          stringr::str_to_upper(AlleleB_Probe_Sequence)==stringr::str_to_upper(PRB1_M_IUP) ~ TRUE,
        
        TRUE ~ NA)
    ) %>% dplyr::group_by(SeqA_Match) %>% dplyr::summarise(n())
  
  rp_tib %>% tidyr::separate(Assay_Design_Id, into=c('Seq_ID', 'Strand_FR', 'Strand_TB', 'Strand_CO', 'Infinium'), sep='_', remove=FALSE) %>%
    dplyr::left_join(fin_ord_tib, by='Seq_ID', suffix=c('_1', '_2')) %>% 
    dplyr::mutate(
      Seq_Match=case_when(
        Infinium=='II' & stringr::str_to_upper(AlleleA_Probe_Sequence_1)==stringr::str_to_upper(AlleleA_Probe_Sequence_2) ~ TRUE,
        Infinium=='I' & 
          stringr::str_to_upper(AlleleA_Probe_Sequence_1)==stringr::str_to_upper(AlleleA_Probe_Sequence_2) &
          stringr::str_to_upper(AlleleB_Probe_Sequence_1)==stringr::str_to_upper(AlleleB_Probe_Sequence_2) ~ TRUE,
        
        TRUE ~ NA),
      Bin_Match=case_when(
        Normalization_Bin_1==Normalization_Bin_2 ~ TRUE,
        TRUE ~ FALSE
      )
    ) %>% dplyr::group_by(Seq_Match,Bin_Match) %>% dplyr::summarise(n())
}


# End of file
