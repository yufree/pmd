#' read in MSP file as list for ms/ms annotation
#' @param file the path to your MSP file
#' @param digits mass or mass to charge ratio accuracy for pmd, default 2
#' @param icf intensity cutoff, default 10 percentage
#' @return list a list with MSP information for MS/MS annotation
#' @export
getms2pmd <- function(file, digits = 2, icf = 10) {
        # this part is modified from compMS2Miner's code: https://github.com/WMBEdmands/compMS2Miner/blob/ee20d3d632b11729d6bbb5b5b93cd468b097251d/R/metID.matchSpectralDB.R
        msp <- readLines(file)
        # remove empty lines
        msp <- msp[msp != '']
        ncomp <- grep('^NAME:', msp, ignore.case = TRUE)
        splitFactorTmp <-
                rep(seq_along(ncomp), diff(c(ncomp, length(msp) + 1)))

        li <- split(msp, f = splitFactorTmp)
        getmsp <- function(x) {
                namet <- x[grep('^NAME:', x, ignore.case = TRUE)]
                name <-
                        gsub('^NAME: ', '', namet, ignore.case = TRUE)
                prect <-
                        x[grep(
                                '^PRECURSORMZ: |^PRECURSOR M/Z: |^PRECURSOR MZ: |^PEPMASS: ',
                                x,
                                ignore.case = TRUE
                        )]
                prec <-
                        as.numeric(
                                gsub(
                                        '^PRECURSORMZ: |^PRECURSOR M/Z: |^PRECURSOR MZ: |^PEPMASS: ',
                                        '',
                                        prect,
                                        ignore.case = TRUE
                                )
                        )
                npt <-
                        x[grep('^Num Peaks: ', x, ignore.case = TRUE)]
                np <-
                        gsub('^Num Peaks: ', '', npt, ignore.case = TRUE)
                if (as.numeric(np) > 0) {
                        # matrix of masses and intensities
                        massIntIndx <-
                                which(grepl('^[0-9]', x) &
                                              !grepl(': ', x))
                        massesInts <-
                                unlist(strsplit(x[massIntIndx], '\t| '))
                        massesInts <-
                                as.numeric(massesInts[grep('^[0-9].*[0-9]$|^[0-9]$',
                                                           massesInts)])
                        # if any NAs remove from indx
                        mz <-
                                massesInts[seq(1, length(massesInts), 2)]
                        ins <-
                                massesInts[seq(2, length(massesInts), 2)]
                        ins <- ins / max(ins) * 100
                        msms <- cbind.data.frame(mz = mz, ins = ins)
                        msms <- msms[msms$ins > icf,]
                        dis <-
                                stats::dist(msms$mz, method = "manhattan")
                        diff <-
                                round(as.numeric(dis), digits = digits)
                        diff <- diff[order(diff)]
                        return(list(
                                name = name,
                                prec = prec,
                                msms = msms,
                                pmd = diff
                        ))
                } else{
                        return(list(
                                name = name,
                                prec = prec,
                                msms = NULL,
                                pmd = NULL
                        ))
                }

        }
        li <- lapply(li, getmsp)
        name <- vapply(li, function(x)
                x$name,'c')
        mz <- vapply(li, function(x)
                x$prec,1)
        msms <- lapply(li, function(x)
                x$pmd)
        msmsraw <- lapply(li, function(x)
                x$msms)
        return(list(
                name = unname(name),
                mz = unname(mz),
                msms = unname(msms),
                msmsraw = unname(msmsraw)
        ))
}

#' read in MSP file as list for EI-MS annotation
#' @param file the path to your MSP file
#' @param digits mass or mass to charge ratio accuracy for pmd, default 0
#' @param icf intensity cutoff, default 10 percentage
#' @return list a list with MSP information for EI-MS annotation
#' @export
getmspmd <- function(file, digits = 2, icf = 10) {
        # this part is modified from compMS2Miner's code: https://github.com/WMBEdmands/compMS2Miner/blob/ee20d3d632b11729d6bbb5b5b93cd468b097251d/R/metID.matchSpectralDB.R
        msp <- readLines(file)
        # remove empty lines
        msp <- msp[msp != '']
        ncomp <- grep('^NAME:', msp, ignore.case = TRUE)
        splitFactorTmp <-
                rep(seq_along(ncomp), diff(c(ncomp, length(msp) + 1)))

        li <- split(msp, f = splitFactorTmp)
        getmsp <- function(x) {
                namet <- x[grep('^NAME:', x, ignore.case = TRUE)]
                name <-
                        gsub('^NAME: ', '', namet, ignore.case = TRUE)
                npt <-
                        x[grep('^Num Peaks: ', x, ignore.case = TRUE)]
                np <-
                        gsub('^Num Peaks: ', '', npt, ignore.case = TRUE)
                if (as.numeric(np) > 0) {
                        # matrix of masses and intensities
                        massIntIndx <-
                                which(grepl('^[0-9]', x) &
                                              !grepl(': ', x))
                        massesInts <-
                                unlist(strsplit(x[massIntIndx], '\t| '))
                        massesInts <-
                                as.numeric(massesInts[grep('^[0-9].*[0-9]$|^[0-9]$',
                                                           massesInts)])
                        # if any NAs remove from indx
                        mz <-
                                massesInts[seq(1, length(massesInts), 2)]
                        ins <-
                                massesInts[seq(2, length(massesInts), 2)]
                        ins <- ins / max(ins) * 100
                        msms <- cbind.data.frame(mz = mz, ins = ins)
                        msms <- msms[msms$ins > icf,]
                        dis <-
                                stats::dist(msms$mz, method = "manhattan")
                        diff <-
                                round(as.numeric(dis), digits = digits)
                        diff <- diff[order(diff)]
                        return(list(
                                name = name,
                                msms = msms,
                                pmd = diff
                        ))
                } else{
                        return(list(
                                name = name,
                                msms = NULL,
                                pmd = NULL
                        ))
                }
        }
        li <- lapply(li, getmsp)
        name <- vapply(li, function(x)
                x$name,'v')
        msms <- lapply(li, function(x)
                x$pmd)
        msmsraw <- lapply(li, function(x)
                x$msms)
        return(list(
                name = unname(name),
                msms = unname(msms),
                msmsraw = unname(msmsraw)
        ))
}

#' Perform MS/MS pmd annotation for mgf file
#' @param file mgf file generated from MS/MS data
#' @param db database could be list object from `getms2pmd`
#' @param ppm mass accuracy, default 10
#' @param prems precursor mass range, default 1.1 to include M+H or M-H
#' @param pmdc pmd length percentage cutoff for annotation. 0.6(default) means 60 percentage of the pmds in your sample could be found in certain compound pmd database
#' @param scutoff relative intensity cutoff for input spectra for pmd analysis, default 0.1
#' @return list with MSMS annotation results
#' @export
pmdanno <- function(file,
                    db = NULL,
                    ppm = 10,
                    prems = 1.1,
                    pmdc = 0.6,
                    scutoff = 0.1) {
        namemgf <- basename(file)
        sample <- MSnbase::readMgfData(file)
        prec <- MSnbase::precursorMz(sample)
        mz <- MSnbase::mz(sample)
        ins <- MSnbase::intensity(sample)
        idx <- unlist(ins) / max(unlist(ins)) > scutoff
        if (sum(idx) > 0) {
                mz <- unlist(mz)[idx]
                ins <- unlist(ins)[idx]
                ins <- ins / max(ins) * 100
                pmdt <- stats::dist(mz, method = "manhattan")

                if (is.list(db)) {
                        pmdt <-
                                unique(round(as.numeric(pmdt), digits = 2))
                        range <- cbind(db$mz - prems, db$mz + prems)
                        mza <- data.table::as.data.table(cbind.data.frame(mzmin=db$mz - prems,mzmax= db$mz + prems))
                        mzb <- data.table::as.data.table(cbind.data.frame(mzmin=prec - ppm / prec*1e-06,mzmax= prec + ppm / prec*1e-06))
                        colnames(mza) <- colnames(mzb) <- c("min","max")
                        data.table::setkey(mzb, min, max)
                        overlapms <- data.table::foverlaps(mza, mzb, which = TRUE)
                        over <- overlapms[stats::complete.cases(overlapms)&!duplicated(overlapms),]
                        re <- data.frame(over)
                        if (nrow(re) > 0) {
                                msmsd <- db$msms[re$xid]
                                name <- db$name[re$xid]
                                mz2 <- db$mz[re$xid]
                                msmsraw <- db$msmsraw[re$xid]

                                result <- vapply(msmsd, function(x)
                                        sum(pmdt %in% unique(x)),1)
                                result <- c(result)

                                if (sum(result > length(pmdt) * pmdc) == 0) {
                                        return(NULL)
                                } else{
                                        re0 <- result[result > length(pmdt) * pmdc]
                                        order <- order(re0,decreasing = TRUE)
                                        t <-
                                                list(
                                                        name = name[result > length(pmdt) * pmdc][order],
                                                        mz = mz2[result > length(pmdt) * pmdc][order],
                                                        msms = msmsd[result > length(pmdt) * pmdc][order],
                                                        msmsraw = msmsraw[result > length(pmdt) * pmdc][order],
                                                        dmz = mz,
                                                        dprc = prec,
                                                        dins = ins,
                                                        file = namemgf
                                                )
                                        return(t)
                                }
                        } else{
                                return(NULL)
                        }
                } else{
                        return(NULL)
                }
        } else{
                return(NULL)
        }
}
