# Shared MSP parser for getms2pmd / getmspmd.
# `precursor = TRUE` returns the precursor m/z; set FALSE for EI-MS.
.parse_msp <- function(file, digits = 2, icf = 10, precursor = TRUE) {
        # adapted from compMS2Miner:
        # https://github.com/WMBEdmands/compMS2Miner/blob/ee20d3d632b11729d6bbb5b5b93cd468b097251d/R/metID.matchSpectralDB.R
        msp <- readLines(file)
        msp <- msp[msp != '']
        ncomp <- grep('^NAME:', msp, ignore.case = TRUE)
        splitFactorTmp <- rep(seq_along(ncomp),
                              diff(c(ncomp, length(msp) + 1)))
        li <- split(msp, f = splitFactorTmp)

        prec_re <- '^PRECURSORMZ: |^PRECURSOR M/Z: |^PRECURSOR MZ: |^PEPMASS: '

        getmsp <- function(x) {
                name <- gsub('^NAME: ', '',
                             x[grep('^NAME:', x, ignore.case = TRUE)],
                             ignore.case = TRUE)
                prec <- if (precursor) {
                        prect <- x[grep(prec_re, x, ignore.case = TRUE)]
                        as.numeric(gsub(prec_re, '', prect, ignore.case = TRUE))
                } else NA_real_
                np <- as.numeric(gsub('^Num Peaks: ', '',
                                      x[grep('^Num Peaks: ', x, ignore.case = TRUE)],
                                      ignore.case = TRUE))
                if (!length(np) || is.na(np) || np <= 0) {
                        return(list(name = name, prec = prec,
                                    msms = NULL, pmd = NULL))
                }
                # matrix of masses and intensities
                massIntIndx <- which(grepl('^[0-9]', x) & !grepl(': ', x))
                massesInts <- unlist(strsplit(x[massIntIndx], '\t| '))
                massesInts <- as.numeric(
                        massesInts[grep('^[0-9].*[0-9]$|^[0-9]$', massesInts)])
                mz  <- massesInts[seq(1L, length(massesInts), 2L)]
                ins <- massesInts[seq(2L, length(massesInts), 2L)]
                ins <- ins / max(ins) * 100
                msms <- cbind.data.frame(mz = mz, ins = ins)
                msms <- msms[msms$ins > icf, ]
                diff <- round(as.numeric(stats::dist(msms$mz, method = "manhattan")),
                              digits = digits)
                diff <- diff[order(diff)]
                list(name = name, prec = prec, msms = msms, pmd = diff)
        }

        li <- lapply(li, getmsp)
        name    <- vapply(li, function(x) x$name, character(1))
        msms    <- lapply(li, function(x) x$pmd)
        msmsraw <- lapply(li, function(x) x$msms)
        out <- list(name    = unname(name),
                    msms    = unname(msms),
                    msmsraw = unname(msmsraw))
        if (precursor) {
                out <- c(list(name = out$name,
                              mz   = unname(vapply(li, function(x) x$prec, numeric(1)))),
                         out[c("msms", "msmsraw")])
        }
        out
}

#' read in MSP file as list for ms/ms annotation
#' @param file the path to your MSP file
#' @param digits mass or mass to charge ratio accuracy for pmd, default 2
#' @param icf intensity cutoff, default 10 percentage
#' @return list a list with MSP information for MS/MS annotation
#' @export
getms2pmd <- function(file, digits = 2, icf = 10) {
        .parse_msp(file, digits = digits, icf = icf, precursor = TRUE)
}

#' read in MSP file as list for EI-MS annotation
#' @param file the path to your MSP file
#' @param digits mass or mass to charge ratio accuracy for pmd, default 0
#' @param icf intensity cutoff, default 10 percentage
#' @return list a list with MSP information for EI-MS annotation
#' @export
getmspmd <- function(file, digits = 2, icf = 10) {
        .parse_msp(file, digits = digits, icf = icf, precursor = FALSE)
}
