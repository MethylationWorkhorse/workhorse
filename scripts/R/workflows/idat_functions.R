
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                          Basic IDAT Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
 
suppressWarnings(suppressPackageStartupMessages( base::require("illuminaio") ))

suppressWarnings(suppressPackageStartupMessages( base::require("tidyverse") ))
suppressWarnings(suppressPackageStartupMessages( base::require("stringr") ))
suppressWarnings(suppressPackageStartupMessages( base::require("glue") ))

COM <- ","
TAB <- "\t"
RET <- "\n"

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                             Idat File I/O::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

prefixesToChipTib = function(prefixes) {
  basecodes <- names(prefixes)
  
  tib <- NULL
  for (basecode in basecodes) {
    sentrix_codes <- basecode %>% stringr::str_split(pattern='_', n=2, simplify=TRUE)
    tib <- tib %>% dplyr::bind_rows(
      tibble::tibble(barcode=sentrix_codes[1],
                     poscode=sentrix_codes[2],
                     basecode=basecode,
                     path=prefixes[[basecode]])
    )
  }
  
  tib
}

prefixToIdat = function(prefix, load=FALSE, save=FALSE, rds=NULL, gzip=TRUE, validate=TRUE, 
                        verbose=0,vt=4,tc=1,tt=NULL) {
  funcTag <- 'prefixToIdat'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Loading prefix={prefix}.{RET}"))
  
  stime <- system.time({
    if (load && !is.null(rds) && file.exists(rds)) {
      if (verbose>=vt+1) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Loading RDS={rds}.{RET}"))
      dat <- readr::read_rds(rds)
    } else {
      dat <- NULL
      
      grn_idat <- loadIdat(prefix, 'Grn', gzip=gzip, verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
      red_idat <- loadIdat(prefix, 'Red', gzip=gzip, verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
      
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      #                          Extract Signal Data::
      grn_sig <- getIdatSignalTib(grn_idat, channel='Grn', verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
      red_sig <- getIdatSignalTib(red_idat, channel='Red', verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
      
      sigs <- dplyr::full_join(grn_sig, red_sig, by="Address")
      
      grn_ann <- dplyr::bind_cols(
        getIdatBarcodeTib(grn_idat, verbose=verbose,vt=vt+1,tc=tc+1,tt=tt),
        getIdatFormatTib(grn_idat, verbose=verbose,vt=vt+1,tc=tc+1,tt=tt),
        getIdatTimeStampTib(grn_idat, method='Decoding', verbose=verbose,vt=vt+1,tc=tc+1,tt=tt),
        getIdatTimeStampTib(grn_idat, method='Extract', verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
      )
      
      if (validate) {
        red_ann <- dplyr::bind_cols(
          getIdatBarcodeTib(red_idat, verbose=verbose,vt=vt+1,tc=tc+1,tt=tt),
          getIdatFormatTib(red_idat, verbose=verbose,vt=vt+1,tc=tc+1,tt=tt),
          getIdatTimeStampTib(red_idat, method='Decoding', verbose=verbose,vt=vt+1,tc=tc+1,tt=tt),
          getIdatTimeStampTib(red_idat, method='Extract', verbose=verbose,vt=vt+1,tc=tc+1,tt=tt)
        )
        
        unique_cnt <- dplyr::bind_rows(grn_ann, red_ann) %>% dplyr::distinct() %>% base::nrow()
        stopifnot(unique_cnt==1)
        if (verbose>=vt+1) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Passed idat pair validation!{RET}"))
      }
      
      dat$sig <- sigs
      dat$ann <- grn_ann
      
      if (save && !is.null(rds)) {
        if (verbose>=vt+1) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Writing RDS={rds}.{RET}"))
        readr::write_rds(dat, rds, compress="gz")
      }
    }
  })
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}{RET}"))
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  
  dat
}

loadIdat = function(prefix, col, gzip=TRUE, verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'loadIdat'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Loading prefix={prefix}.{RET}"))
  
  stime <- system.time({
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                        Ensure Idats are Gzipped::
    idat_file <- paste0(prefix,"_",col,".idat")
    if (file.exists(idat_file)) {
      if (gzip) {
        if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr}{TAB} Compressing file={idat_file}.{RET}"))
        system(paste0("gzip ",idat_file))
        idat_file <- paste0(prefix,"_",col,".idat.gz")
      }
    } else {
      idat_file <- paste0(prefix,"_",col,".idat.gz")
    }
    if (!file.exists(idat_file)) stop(glue::glue("{RET}[{funcTag}]: ERROR: idat_file={idat_file} does NOT exist!!!{RET}{RET}"))
    stopifnot(file.exists(idat_file))
    
    idat <- illuminaio::readIDAT(idat_file)
  })
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}{RET}"))
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  
  idat
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                             Idat Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

getIdatSignalTib = function(idat, channel, del='_', verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'getIdatSignalTib'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  stime <- system.time({
    datTag   <- 'sig'
    meanName <- paste('Raw', channel,datTag, sep=del)
    sdName   <- paste('SD',  channel,datTag, sep=del)
    beadName <- paste('Bead',channel, sep=del)
    
    tib <- idat$Quants %>% tibble::as_tibble(rownames="Address") %>%
      dplyr::mutate(Address=as.integer(Address)) %>%
      purrr::set_names('Address',meanName,sdName,beadName)
    tib_nrow <- base::nrow(tib)
  })
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done. Idat nrows={tib_nrow}.{RET}{RET}"))
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  
  tib
}

getIdatBarcodeTib = function(idat, verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'getIdatBarcodeTib'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  stime <- system.time({
    sentrixName <- paste(idat$Barcode,idat$Unknowns$MostlyA, sep='_')
    rowcol_df <- idat$Unknowns$MostlyA %>% stringr::str_remove("^R") %>% stringr::str_split('C', simplify=TRUE, n=2) %>% as.numeric()
    tib <- tibble::tibble(Sentrix_Name=sentrixName,
                          Sentrix_Barcode=idat$Barcode,
                          Sentrix_Poscode=idat$Unknowns$MostlyA,
                          Sentrix_Row=as.integer(rowcol_df[1]),
                          Sentrix_Col=as.integer(rowcol_df[2]) )
  })
  if (verbose>=vt) 
    cat(glue::glue("[{funcTag}]:{tabsStr} Done. Sentrix_Name(R{tib$Sentrix_Row}:C{tib$Sentrix_Col})={tib$Sentrix_Name}{RET}{RET}"))
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  
  tib # %>% gather()
}

getIdatFormatTib = function(idat, verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'getIdatFormatTib'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  stime <- system.time({
    chipType <- idat$ChipType
    if (chipType=='BeadChip 8x5') chipFormat         <- '8x1'
    else if (chipType=='BeadChip 12x8') chipFormat   <- '12x1'
    else if (chipType=='BeadChip 24x1x4') chipFormat <- '24x1'
    else if (chipType=='BeadChip 24x1x2') chipFormat <- '24x1'
    else if (chipType=='Beadchip 24x1x2') chipFormat <- '24x1'
    else stop(glue::glue("{RET}[{funcTag}]: ERROR: Unrecognized ChipType={chipType}!{RET}{RET}"))
    tib <- tibble::tibble('Chip_Type'=chipType,'Chip_Format'=chipFormat)
  })
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Chip_Type={chipType}, Chip_Format={chipFormat}.{RET}{RET}"))
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  
  tib # %>% gather()
}

getIdatTimeStampTib = function(idat, method='Extract', sherlockID='sherlockID', order='latest', 
                               verbose=0,vt=3,tc=1,tt=NULL) {
  funcTag <- 'getIdatTimeStampTib'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting...{RET}"))
  
  stime <- system.time({
    time_tib  <- idat$RunInfo %>% tibble::as_tibble()
    method_idxs <- grep(method, time_tib$BlockType)
    machine_idxs <- grep(sherlockID, time_tib$BlockPars)
    stopifnot(length(method_idxs)>0)
    if (verbose>=vt+5) print(method_idxs)
    
    stopifnot(order=='latest')
    if (order=='latest')  {
      metod_idx <- method_idxs %>% tail(n=1)
      machine_idx <- machine_idxs %>% tail(n=1)
    } else {
      metod_idx <- method_idxs %>% head(n=1)
      machine_idx <- machine_idxs %>% head(n=1)
    }
    mach_vec <- time_tib$BlockPars[machine_idx] %>% str_split('\\|', simplify=TRUE)
    mach_var <- mach_vec[2] %>% stringr::str_split('=', simplify=TRUE)
    date_str <- time_tib$RunTime[metod_idx] %>% as.POSIXct(format="%m/%d/%Y %H:%M:%S")
    name_str <- paste('Iscan',method, sep='_')
    
    mach_key <- mach_var[1]
    mach_val <- mach_var[2]
    
    time_tib <- date_str %>% stringr::str_split('[\\-\\/ \\:]') %>% BiocGenerics::unlist() %>% 
      purrr::set_names('Year','Mon','Day','Hour','Min','Sec') %>% 
      tibble::enframe() %>% spread(name, value) %>% dplyr::mutate_all(.funs = (as.integer)) %>%
      tibble::add_column('Date' = date_str) %>%
      tibble::add_column(!!mach_key := !!mach_val) %>%
      dplyr::select(!!mach_key, Date, Year, Mon, Day, Hour, Min, Sec) %>%
      purrr::set_names(paste(name_str, names(.), sep='_'))
  })
  if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done. Date({method})={date_str}.{RET}{RET}"))
  if (!is.null(tt)) tt$addTime(stime,funcTag)
  
  time_tib # %>% gather()
}


# End of file
