
suppressWarnings(suppressPackageStartupMessages( base::require("tidyverse") ))
suppressWarnings(suppressPackageStartupMessages( base::require("stringr") ))
suppressWarnings(suppressPackageStartupMessages( base::require("glue") ))

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                         Common Booster Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

joinTibbles = function(a, b, by, side="left", verbose=0,vt=3,tc=1) {
  funcTag <- 'joinTibbles'
  tabsStr <- paste0(rep(TAB, tc), collapse='')
  
  # if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Starting; platform={platform}.{RET}"))
  
  if (is.null(a)) return(b)
  if (is.null(b)) return(a)
  
  if (side=='left')  return(dplyr::left_join(a,  b, by=by))
  if (side=='right') return(dplyr::right_join(a, b, by=by))
  if (side=='inner') return(dplyr::inner_join(a, b, by=by))
  if (side=='full')  return(dplyr::full_join(a,  b, by=by))
  stop(glue::glue("[$func]: ERROR: Unsupported join={side}!!!{RET}{RET}"))

  # if (verbose>=vt) cat(glue::glue("[{funcTag}]:{tabsStr} Done.{RET}{RET}"))
  
  NULL
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                     Summarizing Counting Methods::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
cntPer_lte = function(x, min, prc=3) {
  round(100*length(x[which(x<=min)])/length(x),prc)
}

cntPer_gt = function(x, max, prc=3) {
  round(100*length(x[which(x>max)])/length(x),prc)
}

bool2int = function(x) {
  y <- NULL
  for (ii in seq(length(x))) {
    if (is.na(x[ii]) || x[ii]==FALSE) {
      y[ii] <- 0
    } else {
      y[ii] <- 1
    }
  }
  
  y
}

# End of file
