
suppressWarnings(suppressPackageStartupMessages(require("tidyverse")) )
suppressWarnings(suppressPackageStartupMessages(require("stringr")) )
suppressWarnings(suppressPackageStartupMessages(require("glue")) )
suppressWarnings(suppressPackageStartupMessages(require("scales")) )
suppressWarnings(suppressPackageStartupMessages(require("matrixStats")) )

COM <- ","
TAB <- "\t"
RET <- "\n"

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                             idat Date Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

getIdatTimeStampTib = function(idat, method='Extract', order='latest', 
                               verbose=0, vt=6, tc=1, tabsStr='') {
  funcTag <- 'getIdatTimeStampTib'
  for (i in c(1:tc)) tabsStr <- paste0(tabsStr, TAB)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting..."),"\n", sep='')

  time_tib  <- idat$RunInfo %>% tibble::as_tibble()
  method_idxs <- grep(method, time_tib$BlockType)
  stopifnot(length(method_idxs)>0)
  if (verbose>=vt+5) print(method_idxs)

  stopifnot(order=='latest')
  if (order=='latest') 
    metod_idx <- method_idxs %>% tail(n=1)
  else
    metod_idx <- method_idxs %>% head(n=1)

  date_str <- time_tib$RunTime[metod_idx] %>% as.POSIXct(format="%m/%d/%Y %H:%M:%S")
  name_str <- paste('iscan',method, sep='_')
  
  time_tib <- date_str %>% stringr::str_split('[\\-\\/ \\:]') %>% BiocGenerics::unlist() %>% 
    purrr::set_names('Year','Mon','Day','Hour','Min','Sec') %>% 
    tibble::enframe() %>% spread(name, value) %>% dplyr::mutate_all(.funs = (as.integer)) %>%
    tibble::add_column('Date' = date_str) %>%
    dplyr::select(Date, Year, Mon, Day, Hour, Min, Sec) %>%
    purrr::set_names(paste(name_str, names(.), sep='_'))
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Date({method})={date_str}"),"\n\n", sep='')

  time_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                          idat ChipFormat Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
getIdatFormatTib = function(idat, verbose=0, vt=6, tc=1, tabsStr='') {
  funcTag <- 'getIdatFormatTib'
  for (i in c(1:tc)) tabsStr <- paste0(tabsStr, TAB)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting..."),"\n", sep='')
 
  chipType <- idat$ChipType
  if (chipType=='BeadChip 8x5') chipFormat         <- '8x1'
  else if (chipType=='BeadChip 12x8') chipFormat   <- '12x1'
  else if (chipType=='BeadChip 24x1x4') chipFormat <- '24x1'
  else if (chipType=='BeadChip 24x1x2') chipFormat <- '24x1'
  else if (chipType=='Beadchip 24x1x2') chipFormat <- '24x1'
  else {
    stop("\n",glue::glue("[{funcTag}]: Fatal ERROR: Unrecognized ChipType={chipType}"),"\n\n", sep='')
    q()
  }
  tib <- tibble::tibble('ChipType'=chipType,
                        'ChipFormat'=chipFormat)
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} ChipType={chipType}, ChipFormat={chipFormat}"),"\n\n", sep='')

  tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                          idat ChipBarcode Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

getIdatBarcodeTib = function(idat, verbose=0, vt=6, tc=1, tabsStr='') {
  funcTag <- 'getIdatBarcodeTib'
  for (i in c(1:tc)) tabsStr <- paste0(tabsStr, TAB)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting..."),"\n", sep='')
  
  sentrixName <- paste(idat$Barcode,idat$Unknowns$MostlyA, sep='_')
  rowcol_df <- idat$Unknowns$MostlyA %>% stringr::str_remove("^R") %>% stringr::str_split('C', simplify=TRUE, n=2) %>% as.numeric()

  tib <- tibble::tibble(Sentrix_Name=sentrixName,
                        Sentrix_Barcode=idat$Barcode,
                        Sentrix_Poscode=idat$Unknowns$MostlyA,
                        Sentrix_Row=as.integer(rowcol_df[1]),
                        Sentrix_Col=as.integer(rowcol_df[2]) )
  
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Sentrix_Name(R{tib$Sentrix_Row}:C{tib$Sentrix_Col})={tib$Sentrix_Name}"),"\n\n", sep='')
  
  tib
}

getIdatSignalTib = function(idat, channel, verbose=0, vt=6, tc=1, tabsStr='') {
  funcTag <- 'getIdatBarcodeTib'
  for (i in c(1:tc)) tabsStr <- paste0(tabsStr, TAB)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting..."),"\n", sep='')
  
  meanName  <- paste0(channel,'_Mean')
  sdName    <- paste0(channel,'_SD')
  nbeadName <- paste0(channel,'_NBeads')
  
  tib <- idat$Quants %>% tibble::as_tibble(rownames="Address") %>%
    dplyr::mutate(Address=as.integer(Address)) %>%
    purrr::set_names('Address',meanName,sdName,nbeadName)
  tib_nrow <- base::nrow(tib)
  
  if (verbose>=vt+2) print(tib)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Idat nrows={tib_nrow}"),"\n\n", sep='')

  tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                            Parse idat Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

parseIdat = function(file, channel, verbose=0, vt=5, tc=1, tabsStr='') {
  funcTag <- 'parseIdat'
  for (i in c(1:tc)) tabsStr <- paste0(tabsStr, TAB)

  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{TAB} Starting Parse, channel={channel}, file={file}"),"\n", sep='')
  stopifnot(file.exists(file))
  
  dat  <- NULL
  idat <- illuminaio::readIDAT(file)
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                          Extract Signal Data::
  signal_tib <- getIdatSignalTib(idat, channel=channel, verbose=verbose,vt=vt+1,tc=tc+1)
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                            Add Chip Barcodes::
  barcode_tib <- getIdatBarcodeTib(idat, verbose=verbose,vt=vt+1,tc=tc+1)
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                           Add Valid Chip Types::
  format_tib <- getIdatFormatTib(idat, verbose=verbose,vt=vt+1,tc=tc+1)
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                        Add Decoding/Extraction Dates::
  decodeDate_tib  <- getIdatTimeStampTib(idat, method='Decoding', verbose=verbose,vt=vt+1,tc=tc+1)
  extractDate_tib <- getIdatTimeStampTib(idat, method='Extract', verbose=verbose,vt=vt+1,tc=tc+1)

  dat$sig <- signal_tib
  dat$sss <- dplyr::bind_cols(barcode_tib, format_tib, decodeDate_tib, extractDate_tib)

  dat
}

prefixToTib = function(prefix, verbose=0,vt=4,tc=1,tabsStr='') {
  funcTag <- 'prefixToTib'
  for (i in c(1:tc)) tabsStr <- paste0(tabsStr, TAB)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting prefix={prefix}"),"\n", sep='')
  
  dat <- NULL
  for (channel in c('Grn','Red')) {
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                        Ensure Idats are Gzipped::
    channel.file <- paste0(prefix,"_",channel,".idat")
    if (file.exists(channel.file)) system(paste0("gzip ",channel.file))
    channel.file <- paste0(prefix,"_",channel,".idat.gz")
    stopifnot(file.exists(channel.file))
    
    if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Channel={channel}, file={channel.file}"),"\n", sep='')
    cur <- parseIdat(channel.file, channel, verbose=verbose, vt=vt+1, tc=tc+1)
    
    if (is.null(dat)) {
      dat <- cur
    } else {
      if (dat$sss$ChipType != cur$sss$ChipType) 
        stop("\n",glue::glue("[{funcTag}]: ERROR ChipType's don't match {dat$sss$ChipType} != {cur$sss$ChipType}"),"\n", sep='')
      if (dat$sss$ChipFormat != cur$sss$ChipFormat) 
        stop("\n",glue::glue("[{funcTag}]: ERROR ChipFormat's don't match {dat$sss$ChipFormat} != {cur$sss$ChipFormat}"),"\n", sep='')
      if (dat$sss$Sentrix_Name != cur$sss$Sentrix_Name)
        stop("\n",glue::glue("[{funcTag}]: ERROR Barcode's don't match {dat$sss$Sentrix_Name} != {cur$sss$Sentrix_Name}"),"\n", sep='')
      
      dat$sig <- dat$sig %>% dplyr::full_join(cur$sig, by="Address")
    }
  }
  if (verbose>=vt) cat("\t",glue::glue("[{funcTag}]:{tabsStr} Done"),"\n\n", sep='')
  
  dat
}

fetchPrefix = function(prefix, sig_csv, sss_csv, loadIdat=FALSE, saveIdat=FALSE,
                       verbose=0,vt=4,tc=1,tabsStr='') {
  funcTag <- 'fetchPrefix'
  for (i in c(1:tc)) tabsStr <- paste0(tabsStr, TAB)
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting prefix={prefix}"),"\n", sep='')
  
  if (!file.exists(sig_csv) || !file.exists(sss_csv)) loadIdat <- FALSE
  
  dat <- NULL
  if (loadIdat && file.mtime(sig_csv) < file.mtime(sss_csv)) {
    dat$sig <- suppressMessages(suppressWarnings(readr::read_csv(sig_csv) ))
    dat$sss <- suppressMessages(suppressWarnings(readr::read_csv(sss_csv) ))
  } else {
    dat <- prefixToTib(prefix=prefix, verbose=0,vt=4,tc=1,tabsStr='')
    # print(dat)
    
    if (saveIdat) {
      if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Saving(CSV)={sig_csv}"),"\n", sep='')
      readr::write_csv(dat$sig, sig_csv)
      if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Saving(CSV)={sss_csv}"),"\n", sep='')
      readr::write_csv(dat$sss, sss_csv)
    }
  }
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done."),"\n\n", sep='')

  dat
}

# prefix_path <- "/Users/bbarnes/Documents/Projects/workhorse/idats/ReferenceBETA/201502830033/201502830033_R01C01"
# prefixToTib(prefix_path, verbose=5)

# End of file
