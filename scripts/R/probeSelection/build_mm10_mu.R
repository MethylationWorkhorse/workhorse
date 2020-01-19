
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
  
  opt$outDir <- file.path(opt$topDir, 'workspace/improbe-mu-mm10')
  
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
#                      Load target multi-unique probes::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
imp_des_tsv <- '/Users/bbarnes/Documents/Projects/workhorse/designs/target-mm10.mu.improbe-output.tsv.gz'
imp_raw_tib <- suppressMessages(suppressWarnings( readr::read_tsv(imp_des_tsv) )) %>% dplyr::mutate(Chromosome=as.character(paste0('chr',Chromosome)))

# imp_des_tib <- loadImprobeDesign(file=imp_des_tsv) %>% dplyr::mutate(Chromosome=as.character(paste0('chr',Chromosome)))
mup_des_tib <- loadImprobeDesign(src_des_tib=imp_raw_tib)

imp_des_tib <- imp_raw_tib %>% dplyr::arrange(Seq_ID, Methyl_Allele_FR_Strand, Methyl_Allele_TB_Strand, Methyl_Allele_CO_Strand) %>% 
  dplyr::mutate(Probe_ID=paste(Seq_ID, Methyl_Allele_FR_Strand, stringr::str_sub(Methyl_Allele_TB_Strand,1,1), Methyl_Allele_CO_Strand, sep='_')) %>% 
  dplyr::distinct(Probe_ID, .keep_all=TRUE) 

imp_cnt_tib <- mup_des_tib %>% dplyr::select(Probe_ID, Chromosome, Coordinate) %>% dplyr::distinct() %>% dplyr::group_by(Probe_ID) %>% dplyr::summarise(IMP_Counts=n())
imp_des_cnt <- mup_des_tib %>% base::nrow()
cat(glue::glue("[{opt$prgmTag}]: imp_des_cnt={imp_des_cnt}.{RET}"))

# Don't think we need this anymore with the full original improbe designs above...
#
# mu_des_tsv <- '/Users/bbarnes/Documents/Projects/workhorse/designs/mm10-2020.data.cgn-sorted.bed.gz'
# mu_des_tib <- suppressMessages(suppressWarnings( readr::read_tsv(mu_des_tsv, col_names=FALSE) )) %>% 
#   purrr::set_names('Chromosome', 'start', 'end', 'Probe_ID', 'scrMin', 'Strand_FR', 'cpgCnt', 'scrU', 'scrM', 'tmU', 'tmM', 'prbU', 'prbM') %>%
#   dplyr::mutate(scrMin=as.integer(scrMin*10))
# mu_des_cnt  <- mu_des_tib %>% base::nrow()
# cat(glue::glue("[{opt$prgmTag}]:  mu_des_cnt={mu_des_cnt}.{RET}"))

mu_aln_tsv <- '/Users/bbarnes/Documents/Projects/workhorse/designs/cpg-hg19.s12.v4.p16.n1.r2.41.less.noOF-sorted.CGI-intersect.bed'
mu_aln_tib <- suppressMessages(suppressWarnings(readr::read_tsv(mu_aln_tsv, col_names=FALSE) )) %>% 
  purrr::set_names('Chromosome', 'start', 'end', 'Probe_ID', 'scrStr', 'Strand_BSP', 
                   'qSeq', 'aSeq', 'alnTag', 'misCntA', 'misCntB', 'cgiChr', 'cgiBeg', 'cgiEnd', 'cgiKey') %>% 
  tidyr::separate(scrStr, into=c('hitCnt0', 'hitCnt1', 'hitCnt2', 'hitCnt3', 'hitCnt4'), sep=':', remove=FALSE) %>% dplyr::distinct()
mu_aln_cnt  <- mu_aln_tib %>% base::nrow()
cat(glue::glue("[{opt$prgmTag}]:  mu_aln_cnt={mu_aln_cnt}.{RET}"))

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                      Data Merging and Filtering::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
qc_too_stringent <- FALSE
if (qc_too_stringent) {
  # 1. Filter Probes based on Alignment Count vs. Expected Design Count::
  # 2. Add back UNIQUE design metrics (probes, scores, non-positional stuff)
  mu_filt_tib <- mu_aln_tib %>% dplyr::inner_join(imp_cnt_tib, by="Probe_ID") %>% dplyr::filter(hitCnt0==IMP_Counts) %>%
    dplyr::left_join(dplyr::distinct(dplyr::select(mup_des_tib, -Chromosome, -Coordinate)), by="Probe_ID")
  
  # 3. Sanity Check that all probe sequences are unique after removing redundant fields::
  qc_seq_unique <- FALSE
  if (qc_seq_unique) {
    mu_filt_tib %>% dplyr::select(-Chromosome, -start, -end, -Strand_BSP, -qSeq, -aSeq, -cgiChr, -cgiBeg, -cgiEnd, -cgiKey) %>% dplyr::distinct() %>% 
      group_by(PRB1_U) %>% dplyr::summarise(Prb1U_Count=n()) %>% dplyr::arrange(-Prb1U_Count)
    mu_filt_tib %>% dplyr::select(-Chromosome, -start, -end, -Strand_BSP, -qSeq, -aSeq, -cgiChr, -cgiBeg, -cgiEnd, -cgiKey) %>% dplyr::distinct() %>% 
      group_by(PRB1_M) %>% dplyr::summarise(Prb1M_Count=n()) %>% dplyr::arrange(-Prb1M_Count)
  }
}

# This is doing #1 above, but with out the filtering. #2 (joining of design metrics) is still performed
mu_filt_tib <- mu_aln_tib %>% dplyr::inner_join(imp_cnt_tib, by="Probe_ID") %>%
  dplyr::left_join(dplyr::distinct(dplyr::select(mup_des_tib, -Chromosome, -Coordinate)), by="Probe_ID") %>%
  dplyr::mutate(PRB_SCR_MIN=pmin(PRB_SCR_U,PRB_SCR_M))

# 4. Look at coverage of Islands/Probes::
#    CONCLUSION: We loose too much stuff above. Let's keep everything with
#    A. Only 2 hits and amazing scores or CGI selection
#    B. Everything else should be 3 or greater
mu_cgi_sorted_tib <- mu_filt_tib %>% dplyr::group_by(cgiKey) %>% summarise(CGI_Cov_Count=n()) %>% dplyr::arrange(-CGI_Cov_Count)
mu_prb_sorted_tib <- mu_filt_tib %>% dplyr::group_by(Probe_ID) %>% summarise(prbRep_Count=n()) %>% dplyr::arrange(-prbRep_Count)
mu_prbCGI_Count_tib <- mu_filt_tib %>% dplyr::group_by(cgiKey, Probe_ID) %>% summarise(cgiPrb_Count=n()) %>% group_by(Probe_ID) %>% 
  summarise(prbCGI_Count=n()) %>% dplyr::arrange(-prbCGI_Count)

# mu_prbCGI_Count_tib %>% dplyr::left_join(mu_prb_sorted_tib, by="Probe_ID") %>% dplyr::left_join(mu_filt_tib, by="Probe_ID") %>%
#   dplyr::arrange(-prbCGI_Count, -prbRep_Count) %>% 
#   dplyr::select(Probe_ID, prbCGI_Count, prbRep_Count, Chromosome, start, PRB_SCR_MIN, PRB_SCR_U, PRB_SCR_M, TM_RAW_U, TM_RAW_M) %>%
#   as.data.frame()

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                             Group Selection::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# Group1::top prb/cgi that looks really good
#   - 4 unique Probe_ID's
#   - size variants Infinium IM: 4+2
grp1_tib <- mu_prbCGI_Count_tib %>% dplyr::left_join(mu_prb_sorted_tib, by="Probe_ID") %>% dplyr::left_join(mu_filt_tib, by="Probe_ID") %>%
  dplyr::arrange(-prbCGI_Count, -prbRep_Count) %>% 
  dplyr::select(Probe_ID, prbCGI_Count, prbRep_Count, Chromosome, start, end, PRB_SCR_MIN, PRB_SCR_U, PRB_SCR_M, TM_RAW_U, TM_SCR_U, TM_RAW_M, TM_SCR_M) %>%
  head(n=12)

# Group2:: MHC (HLA for mice)
#  - 6 unique Probe_ID's
#   - size variants Infinium IM: 6+3
grp2_tib <- mu_prbCGI_Count_tib %>% dplyr::left_join(mu_prb_sorted_tib, by="Probe_ID") %>% dplyr::left_join(mu_filt_tib, by="Probe_ID") %>%
  dplyr::arrange(-prbCGI_Count, -prbRep_Count) %>% 
  dplyr::select(Probe_ID, prbCGI_Count, prbRep_Count, Chromosome, start, end, PRB_SCR_MIN, PRB_SCR_U, PRB_SCR_M, TM_RAW_U, TM_SCR_U, TM_RAW_M, TM_SCR_M) %>%
  dplyr::filter(Chromosome=='chr17' & start >= 33681276 & end <= 38548659)

# Group3:: top prb rankes (repetitive in serial position)
#  - 20 unique Probe_ID's
#  - Need to make most of these Infinium II
grp3_tib <- mu_prbCGI_Count_tib %>% dplyr::left_join(mu_prb_sorted_tib, by="Probe_ID") %>% dplyr::left_join(mu_filt_tib, by="Probe_ID") %>%
  dplyr::arrange(-prbRep_Count,-prbCGI_Count) %>% 
  dplyr::select(Probe_ID, prbCGI_Count, prbRep_Count, Chromosome, start, end, PRB_SCR_MIN, PRB_SCR_U, PRB_SCR_M, TM_RAW_U, TM_SCR_U, TM_RAW_M, TM_SCR_M) %>%
  head(n=80)
# head(n=74)

mup_sel_tib <- dplyr::bind_rows(grp1_tib, grp2_tib, grp3_tib)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                       Further Selection for Trimming::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
imp_des_tib <- dplyr::distinct(mup_sel_tib, Probe_ID) %>% 
  dplyr::left_join(imp_des_tib, by="Probe_ID") %>% dplyr::select(Seq_ID:UnMethyl_Next_Base_Score)

src_man_tib <- imp_des_tib %>% 
  dplyr::mutate(
    Sequence=Forward_Sequence,
    CpG_Island='FALSE',
    PRB_CGN=Seq_ID,
    PRB_DES='cg',
    diNUC='CG',
    IUPAC_Forward_Sequence=Sequence,
    Strand_FR=Methyl_Allele_FR_Strand,
    diSNP=diNUC,
    Infinium_Design_Type='II') %>%
  dplyr::select(Seq_ID,Sequence,Genome_Build,Chromosome,Coordinate,CpG_Island,
                PRB_CGN,PRB_DES,diNUC,IUPAC_Forward_Sequence,Strand_FR,diSNP,Infinium_Design_Type) # %>%
  # dplyr::filter(Seq_ID!='cg29215192')
  # dplyr::filter(Probe_ID != 'cg29215192_F_T_C')
  # 29215192

# This is the annotaiton file need for silly LIMS::
imp_ann_tib <- imp_des_tib
# We can now switch to the abbreviated version::
imp_des_tib <- mup_des_tib

# 26 Infinium II (ALL)
#    mup_sel_tib %>% dplyr::distinct(Probe_ID) %>% base::nrow()
#  8 Infinium I  (grp1)
# 12 Infinium I  (grp2)
# 
# 11 For size selection::
#
#  3 grp1 (-5)
#  3 grp1 (-10)
#  2 grp2 (-5)
#  2 grp2 (-10)
#  9 grp3 (-5)

# Criteria for smaller Infinium IM Probes::
mup_tp2_tib <- mup_sel_tib %>% dplyr::filter( !(TM_RAW_M>60 & PRB_SCR_U-PRB_SCR_M > 0.1)) %>% distinct(Probe_ID)
mup_tp1_tib <- mup_sel_tib %>% dplyr::filter(TM_RAW_M>60 & PRB_SCR_U-PRB_SCR_M > 0.1) %>% distinct(Probe_ID)

grp1_inf1_tib <- grp1_tib %>% dplyr::filter(TM_RAW_M>60 & PRB_SCR_U-PRB_SCR_M > 0.1) %>% distinct(Probe_ID)
# grp2_tib %>% dplyr::filter(TM_RAW_M>60 & PRB_SCR_U-PRB_SCR_M > 0.1) %>% distinct(Probe_ID)
grp2_inf1_tib <- grp2_tib %>% distinct(Probe_ID)
grp3_inf1_tib <- grp3_tib %>% dplyr::filter(TM_RAW_M>60 & PRB_SCR_U-PRB_SCR_M > 0.1) %>% distinct(Probe_ID)

mu_inf1_tib <- dplyr::bind_rows(grp1_inf1_tib,grp2_inf1_tib,grp3_inf1_tib) %>% dplyr::distinct(Probe_ID)
mu_cut5_tib <- grp2_tib %>% dplyr::filter(TM_RAW_M>60 & PRB_SCR_U-PRB_SCR_M > 0.1) %>% 
  distinct(Probe_ID) %>% dplyr::bind_rows(grp1_inf1_tib)
mu_cut9_tib <- mu_cut5_tib %>% 
  dplyr::filter(Probe_ID!='cg43772155_R_B_C') # %>%
  # dplyr::filter(Probe_ID!='cg43772154_F_T_C') %>%
  # dplyr::filter(Probe_ID!='cg34872457_F_T_C')

if (FALSE) {
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
for (ii in seq(1,plot_nrows)) {
  tag_tib <- plot_ord_tib[ii,]
  # if (tag_tib$Seq_ID != 'rs33779096_GY_R') next
  # if (tag_tib$Seq_ID != 'rs220086268_GR_R') next
  print(tag_tib)
  
  cur_tib <- prb_mat_tibs %>% dplyr::filter(Seq_ID==tag_tib$Seq_ID[1])
  tibs <- cur_tib %>% srdsToBrac()
  
  # Build Printable Strings for reach strand::
  IlmnId <- cur_tib$Seq_ID
  cat(glue::glue("ii={ii}; IlmnID={IlmnId}:{RET}"))
  strs <- NULL
  strs$F <- prbsToStr(tibs$F, pr=tag_tib$PRB_DES[1], verbose=opt$verbosity)
  strs$R <- prbsToStr(tibs$R, pr=tag_tib$PRB_DES[1], verbose=opt$verbosity)
  
  # cur_tib <- mat_tib %>% dplyr::filter(IlmnID==IlmnId)
  if (FALSE && cur_tib$Infinium_Design_Type[1]=='II') {
    cur_tib %>% dplyr::select(FR,CO, PRB2_D_IUP, PRB2_D_IMP,Man_MisMatch,Man_TarMatch,Bad_Design) %>% print()
  } else {
    cur_tib %>% dplyr::select(FR,CO, PRB1_U_IUP, PRB1_U_IMP,Man_MisMatch,Man_TarMatch,Bad_Design) %>% print()
    cur_tib %>% dplyr::select(FR,CO, PRB1_M_IUP, PRB1_M_IMP,Man_MisMatch,Man_TarMatch,Bad_Design) %>% print()
  }
  
  break
}

# Selection::
#  1. Pick all !Bad_Design with !is.na(PRB1_U_IMP)
# top_pick_tib  <- prb_mat_tibs %>% dplyr::filter(Bad_Design==FALSE) %>% dplyr::filter(!is.na(PRB1_U_IMP))
top_pick_tib  <- prb_mat_tibs %>% dplyr::filter(Bad_Design==FALSE) %>% dplyr::filter(!is.na(PRB1_U_IMP)) %>% 
  dplyr::inner_join(dplyr::distinct(dplyr::select(mup_sel_tib,Probe_ID)), by="Probe_ID") %>% dplyr::select(-Strand_FR) %>% dplyr::distinct()

top_pick_cgns <- top_pick_tib %>% dplyr::select(Seq_ID) %>% distinct(Seq_ID)
#  2. Pick all !Bad_Design with is.na(PRB1_U_IMP)
mis_pick_cgns <- prb_mat_tibs %>% dplyr::anti_join(top_pick_cgns, by=c('Seq_ID') ) %>% distinct(Seq_ID)

sec_pick_tib <- prb_mat_tibs %>% dplyr::filter(CO==TRUE & Bad_Design==FALSE & Valid_Extension==TRUE) %>%
  dplyr::anti_join(top_pick_cgns, by=c('Seq_ID') ) # %>% 
# group_by(Seq_ID) %>% summarise(CNT=n()) %>% dplyr::arrange(-CNT)

fin_inf2_tib <- top_pick_tib %>% dplyr::bind_rows(sec_pick_tib) %>% dplyr::arrange(Seq_ID)
fin_inf1_tib <- mu_inf1_tib %>% dplyr::left_join(fin_inf2_tib, by="Probe_ID") %>% dplyr::mutate(Infinium_Design_Type='I')

fin_cut5_tib <- mu_cut5_tib %>% dplyr::left_join(fin_inf2_tib, by="Probe_ID") %>% 
  dplyr::mutate(Infinium_Design_Type='I_45M',
                # Probe_ID=paste(Probe_ID,'45M', sep='_'),
                PRB1_U_IUP=stringr::str_sub(PRB1_U_IUP,5),
                PRB1_M_IUP=stringr::str_sub(PRB1_M_IUP,5) )

fin_cut9_tib <- mu_cut9_tib %>% dplyr::left_join(fin_inf2_tib, by="Probe_ID") %>% 
  dplyr::mutate(Infinium_Design_Type='I_41M',
                # Probe_ID=paste(Probe_ID,'41M', sep='_'),
                PRB1_U_IUP=stringr::str_sub(PRB1_U_IUP,9),
                PRB1_M_IUP=stringr::str_sub(PRB1_M_IUP,9) )

fin_pick_tib <- bind_rows(fin_inf2_tib, fin_inf1_tib, fin_cut5_tib, fin_cut9_tib) %>% dplyr::arrange(Probe_ID)

# Assay_Design_Id, AlleleA_Probe_Id, AlleleA_Probe_Sequence, AlleleB_Probe_Id, AlleleB_Probe_Sequence, Normalization_Bin
fin_ord_tib <- fin_pick_tib %>% dplyr::mutate(
  FR_Char=case_when(FR ~ 'F', TRUE ~ 'R'),
  TB_Char=case_when(TB ~ 'T', TRUE ~ 'B'),
  CO_Char=case_when(CO ~ 'C', TRUE ~ 'O'),
  # ID_Char=case_when(
  #   PRB_DES=='rs' & !stringr::str_starts(Seq_ID,'rs') ~ paste0('rs',Seq_ID),
  #   TRUE ~ Seq_ID
  # ),
  ID_Char=stringr::str_replace(Probe_ID,'^cg','mu'),
  # Assay_Design_Id=paste(ID_Char,FR_Char,TB_Char,CO_Char,Infinium_Design_Type, sep='_'),
  Assay_Design_Id=paste(ID_Char,Infinium_Design_Type, sep='_'),
  AlleleA_Probe_Id=paste(Assay_Design_Id,'A',sep='_'),
  AlleleA_Probe_Sequence=case_when(
    Infinium_Design_Type=='II'     ~ PRB2_D,
    Infinium_Design_Type=='I_41M'  ~ PRB1_M_IUP,
    Infinium_Design_Type=='I_45M'  ~ PRB1_M_IUP,
    Infinium_Design_Type=='I'      ~ PRB1_U_IUP,
    TRUE ~ NA_character_
  ),
  
  AlleleB_Probe_Id=case_when(
    Infinium_Design_Type=='I' ~ paste(Assay_Design_Id,'B',sep='_'),
    TRUE ~ ''
  ),
  AlleleB_Probe_Sequence=case_when(
    Infinium_Design_Type=='II'     ~ '',
    Infinium_Design_Type=='I_41M'  ~ '',
    Infinium_Design_Type=='I_45M'  ~ '',
    Infinium_Design_Type=='I'      ~ PRB1_M_IUP,
    TRUE ~ NA_character_
  ),
  Normalization_Bin=case_when(
    Infinium_Design_Type=='II' ~ 'C',
    Infinium_Design_Type=='I_41M' ~ 'C',
    Infinium_Design_Type=='I_45M' ~ 'C',
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

opt$writeOrder <- FALSE
if (opt$writeOrder) {
  fin_csv <- file.path(opt$outDir, 'mm10_LEGX_mu_probes.order.csv.gz')
  readr::write_csv(fin_ord_tib %>% dplyr::select(-Seq_ID), fin_csv)
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                Validation/Error Checking of the whole order::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
fin_all_csv <- file.path(opt$topDir, 'workspace/improbe-mu-mm10/mm10_LEGX_nonCpG_probes.Jan16-2020.order.csv')
fin_all_tib <- readr::read_csv(fin_all_csv, skip=8)

mm10_non_csv <- file.path(opt$topDir, 'designs/merged_with_raw_ordered.with-header.non-cpg.csv.gz')
mm10_raw_tib <- readr::read_csv(mm10_non_csv) %>% dplyr::select(Seq_ID:UnMethyl_Next_Base_Score) %>%
  mutate(Probe_Type='ASPE') %>% 
  dplyr::select(Seq_ID, Forward_Sequence, Genome_Build, Chromosome, Coordinate, Design_State, Seq_Length, 
                Forward_CpG_Coord, TB_Strand, Top_Sequence, Top_CpG_Coord, Probe_Type, Probeset_ID, Probeset_Score, everything())
mm10_non_tib <- suppressMessages(suppressWarnings(readr::read_csv(mm10_non_csv) ))

full_cnt <- fin_all_tib %>% base::nrow()
col1_cnt <- fin_all_tib %>% dplyr::distinct(Assay_Design_Id) %>% base::nrow()
col2_cnt <- fin_all_tib %>% dplyr::distinct(AlleleA_Probe_Id) %>% base::nrow()
col3_cnt <- fin_all_tib %>% dplyr::distinct(AlleleA_Probe_Sequence) %>% base::nrow()
col4_cnt <- fin_all_tib %>% dplyr::distinct(AlleleB_Probe_Id) %>% base::nrow()
col5_cnt <- fin_all_tib %>% dplyr::distinct(AlleleB_Probe_Sequence) %>% base::nrow()
col6_cnt <- fin_all_tib %>% dplyr::distinct(Normalization_Bin) %>% base::nrow()

cat(glue::glue("full_cnt={full_cnt}, col1_cnt={col1_cnt}, col2_cnt={col2_cnt}, col3_cnt={col3_cnt}, ",
               "col4_cnt={col4_cnt}, col5_cnt={col5_cnt}, col6_cnt={col6_cnt}.{RET}"))

imp_cpg_full_tsv <- file.path('/Users/bbarnes/Documents/Projects/methylation/LifeEpigentics/to-order/designW.July-9-2019.cpg-only/Mus_musculus.annotation.tsv.gz')
imp_cpg_full_tib <- readr::read_tsv(imp_cpg_full_tsv)
imp_cpg_full_tib <- imp_cpg_full_tib %>% dplyr::mutate(Chromosome=as.character(Chromosome))

mu_ann_tib <- imp_raw_tib %>% dplyr::inner_join(dplyr::select(fin_ord_tib, 'Assay_Design_Id', 'Seq_ID'), by="Seq_ID")

# We need to combine these three::
# 1. CpG Sites
# imp_cpg_full_tib
# 2. Non-CpG Sites
# mm10_raw_tib <- mm10_raw_tib %>% dplyr::mutate(Seq_ID=case_when(
#  stringr::str_replace(Seq_ID, pattern='^([0-9])')
# )
# 3. Mup-CpG Sites
imp_raw_tib <- imp_raw_tib %>% dplyr::mutate(Seq_ID=stringr::str_replace(Seq_ID,'cg','mu'))

inf1_rows <- dplyr::bind_rows(mm10_raw_tib, imp_raw_tib) %>%
  dplyr::mutate(UnMethyl_assay_id=paste(Seq_ID,Methyl_Allele_FR_Strand, stringr::str_sub(Methyl_Allele_TB_Strand, 1,1),Methyl_Allele_CO_Strand,'I', sep='_'))
inf2_rows <- dplyr::bind_rows(mm10_raw_tib, imp_raw_tib) %>%
  dplyr::mutate(UnMethyl_assay_id=paste(Seq_ID,Methyl_Allele_FR_Strand, stringr::str_sub(Methyl_Allele_TB_Strand, 1,1),Methyl_Allele_CO_Strand,'II', sep='_'))

# new_data <- dplyr::bind_rows(imp_cpg_full_tib,inf1_rows,inf2_rows)
cur_data <- dplyr::bind_rows(inf1_rows,inf2_rows) %>% dplyr::mutate(Seq_ID=UnMethyl_assay_id)
org_data <- dplyr::bind_rows(inf1_rows,inf2_rows)
new_data <- dplyr::bind_rows(inf1_rows,inf2_rows)
new_data <- new_data %>% dplyr::mutate(
  Seq_ID=case_when(
    str_starts(Seq_ID, '[0-9]') ~ paste0('rs',Seq_ID), 
    TRUE ~ Seq_ID), 
  UnMethyl_assay_id=Seq_ID)

# new_data_tsv <- file.path(opt$topDir, 'workspace/improbe-mu-mm10/mm10_LEGX_nonCpG_probes.Jan16-2020.annotation.tsv.gz')
# readr::write_csv(new_data, new_data_tsv)

# imp_check_csv <- file.path(opt$topDir, 'workspace/improbe-mu-mm10/mm10_LEGX_nonCpG_probes.Jan16-2020.annotation.tsv.gz')
# imp_check_tib <- readr::read_csv(imp_check_csv)

# imp_ann_full_tsv <- file.path(opt$topDir, 'workspace/improbe-mu-mm10/mm10_LEGX_nonCpG_probes.Jan16-2020.annotation.tsv.gz')
# imp_ann_full_tib <- mm10_des_tib %>% dplyr::bind_rows(imp_des_tib) %>% dplyr::select(Seq_ID:UnMethyl_Next_Base_Score)
# readr::write_tsv(imp_ann_full_tib, imp_ann_full_tsv)

#
# Latest Check::
#
dat_csv <- '/Users/bbarnes/Documents/Projects/workhorse/workspace/improbe-mu-mm10/mm10_LEGX_nonCpG_probes.Jan16-2020.annotation.csv.gz'
dat_tsv <- '/Users/bbarnes/Documents/Projects/workhorse/workspace/improbe-mu-mm10/mm10_LEGX_nonCpG_probes.Jan16-2020.annotation.tsv.gz'
ord_csv <- '/Users/bbarnes/Documents/Projects/workhorse/workspace/improbe-mu-mm10/mm10_LEGX_nonCpG_probes.Jan16-2020.order.csv'
dat_tib <- readr::read_csv(dat_csv)
ord_tib <- readr::read_csv(ord_csv, skip=8)

# 301 Total Missing Data
mis_ord_tib <- ord_tib %>% dplyr::anti_join( (ord_tib %>% dplyr::inner_join(dat_tib, by=c('Assay_Design_Id'='UnMethyl_assay_id')) ), by='Assay_Design_Id')

snp1_ord_tib <- mis_ord_tib %>% dplyr::mutate(
  New_ID=case_when(stringr::str_starts(Assay_Design_Id, 'rs([A-Z0-9]+)_[0-9]+_[A-Z]') ~ stringr::str_remove(Assay_Design_Id, 'rs'),
                   TRUE ~ Assay_Design_Id),
  New_ID=stringr::str_remove(New_ID, '_[FR]_[TB]_[CO]_[I]+')
) %>% 
  dplyr::inner_join(new_data, by=c("New_ID"="UnMethyl_assay_id") ) %>%
  dplyr::distinct(New_ID, .keep_all=TRUE) %>%
  dplyr::mutate(Seq_ID=Assay_Design_Id, UnMethyl_assay_id=Seq_ID) # %>%
  # dplyr::select(Seq_ID:UnMethyl_assay_id)

# Finds 36 of the missing data
snp2_ord_tib <- mis_ord_tib %>% dplyr::mutate(New_ID=stringr::str_remove(Assay_Design_Id,'rs') ) %>% 
  dplyr::inner_join(org_data, by=c('New_ID'='UnMethyl_assay_id')) %>%
  dplyr::mutate(Seq_ID=Assay_Design_Id,UnMethyl_assay_id=Seq_ID)

found_ids <- dplyr::bind_rows(snp2_ord_tib %>% 
                                dplyr::inner_join(ord_tib, by='Assay_Design_Id') %>% 
                                dplyr::distinct(Assay_Design_Id), snp1_ord_tib %>% 
                                dplyr::inner_join(ord_tib, by='Assay_Design_Id') %>% 
                                dplyr::distinct(Assay_Design_Id))

# Remaining missing 26
snp3_ord_tib <- mis_ord_tib %>% dplyr::anti_join(found_ids, by="Assay_Design_Id") %>% 
  dplyr::mutate(New_ID=stringr::str_remove(Assay_Design_Id,'^rs'), 
                New_ID=stringr::str_remove(New_ID,'_41M'), 
                New_ID=stringr::str_remove(New_ID,'_45M')) %>% 
  dplyr::left_join(new_data, by=c('New_ID'='Seq_ID') )

# Data to add to::
add_dat_tib <- dplyr::bind_rows(snp1_ord_tib, snp2_ord_tib, snp3_ord_tib) %>% 
  # dplyr::mutate(Seq_ID=Assay_Design_Id, UnMethyl_assay_id=Assay_Design_Id) %>% 
  dplyr::select(Seq_ID:UnMethyl_assay_id) %>% dplyr::distinct()

# imp_ann_full_tib <- dplyr::bind_rows(add_dat_tib, new_data)
imp_ann_full_tib <- dplyr::bind_rows(add_dat_tib, cur_data)

imp_val_full_tib <- ord_tib %>% dplyr::inner_join(imp_ann_full_tib, by=c("Assay_Design_Id"='UnMethyl_assay_id')) %>% 
  dplyr::mutate(UnMethyl_assay_id=Assay_Design_Id) %>% dplyr::distinct(UnMethyl_assay_id, .keep_all=TRUE) %>% 
  dplyr::select(Seq_ID:UnMethyl_assay_id) %>% dplyr::mutate(UnMethyl_Allele_Type=NA)

imp_ann_full_tsv <- file.path(opt$topDir, 'workspace/improbe-mu-mm10/mm10_LEGX_nonCpG_probes.Jan17-2020.annotation.tsv.gz')
readr::write_tsv(imp_val_full_tib, imp_ann_full_tsv)

# imp_cpg_full_tsv <- file.path('/Users/bbarnes/Documents/Projects/methylation/LifeEpigentics/to-order/designW.July-9-2019.cpg-only/Mus_musculus.annotation.tsv.gz')
imp_cpg_full_tsv <- '/Users/bbarnes/Documents/Projects/workhorse/designs/Mus_musculus.annotation.tsv.gz'
imp_cpg_full_tib <- readr::read_tsv(imp_cpg_full_tsv)
imp_cpg_full_tib <- imp_cpg_full_tib %>% dplyr::mutate(Chromosome=as.character(Chromosome))

imp_combine_tib <- dplyr::bind_rows(imp_val_full_tib, imp_cpg_full_tib) %>% 
  dplyr::mutate(Genome_Build='mm10')
imp_combine_tsv <- '/Users/bbarnes/Documents/Projects/workhorse/designs/Mus_musculus_Jan17-2020.annotation.tsv.gz'
readr::write_tsv(imp_combine_tib, imp_combine_tsv)


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
