#' Perform MS/MS pmd annotation for mgf file
#' @param file mgf file generated from MS/MS data
#' @param db database could be 'qtof', 'qqq', 'orb' data base from hmdb or list object from `getms2pmd`
#' @param ppm mass accuracy, default 10
#' @param prems precersor mass range, default 1.1 to include M+H or M-H
#' @param pmdc pmd length percentage cutoff for annotation. 0.6(default) means 60 percentage of the pmds in your sample could be found in certain compound pmd database
#' @return list with HMDB annotation results
#' @examples
#' file <- system.file("extdata", "challenge-msms.mgf", package = "pmd")
#' anno <- pmdanno(file)
#' @export
pmdanno <- function(file,db='qtof',ppm = 10,prems = 1.1,pmdc=0.6){
        namemgf <- basename(file)
        sample <- MSnbase::readMgfData(file)
        prec <- MSnbase::precursorMz(sample)
        mz <- MSnbase::mz(sample)
        ins <- MSnbase::intensity(sample)
        idx <- unlist(ins)/max(unlist(ins))>0
        if (sum(idx)>0){
                mz <- unlist(mz)[idx]
                ins <- unlist(ins)[idx]
                pmdt <- stats::dist(mz,method = "manhattan")

                if(db == 'qtof'){
                        qtof <- get("qtof")
                        pmdt <- unique(round(as.numeric(pmdt), digits = 2))
                        range <- cbind(qtof$mz-prems,qtof$mz+prems)
                        idxmz <- enviGCMS::getoverlapmass(range,matrix(c(prec-ppm/prec,prec+ppm/prec),nrow = 1))
                        if (sum(idxmz)>0){
                                msms <- qtof$msms[idxmz]
                                name <- qtof$name[idxmz]
                                mz2 <- qtof$mz[idxmz]
                                msmsraw <- qtof$msmsraw[idxmz]
                                i <- function(x) sum(pmdt %in% unique(x))
                                result <- sapply(msms,i)
                                t <- list(name = name[result>length(pmdt)*pmdc],mz = mz2[result>length(pmdt)*pmdc],msms = msms[result>length(pmdt)*pmdc],msmsraw = msmsraw[result>length(pmdt)*pmdc],dmz = mz, dprc=prec, dins = ins/max(ins)*100, file = namemgf)
                                if(sum(result)==0){
                                        return(NULL)
                                }else{
                                        return(t)
                                }
                        }else{
                                return(NULL)
                        }
                }else if(db == 'orb'){
                        orb <- get("orb")
                        pmdt <- unique(round(as.numeric(pmdt), digits = 2))
                        range <- cbind(orb$mz-prems,orb$mz+prems)
                        idxmz <- enviGCMS::getoverlapmass(range,matrix(c(prec-ppm/prec,prec+ppm/prec),nrow = 1))
                        if(sum(idxmz)>0){
                                msms <- orb$msms[idxmz]
                                name <- orb$name[idxmz]
                                mz2 <- orb$mz[idxmz]
                                msmsraw <- orb$msmsraw[idxmz]
                                i <- function(x) sum(pmdt %in% unique(x))
                                result <- sapply(msms,i)
                                t <- list(name = name[result>length(pmdt)*pmdc],mz = mz2[result>length(pmdt)*pmdc],msms = msms[result>length(pmdt)*pmdc],msmsraw = msmsraw[result>length(pmdt)*pmdc],dmz = mz, dprc=prec, dins = ins/max(ins)*100, file = namemgf)
                                if(sum(result)==0){
                                        return(NULL)
                                }else{
                                        return(t)
                                }
                        }else{
                                return(NULL)
                        }
                }else if (db == 'qqq'){
                        qqq <- get("qqq")
                        pmdt <- unique(round(as.numeric(pmdt), digits = 0))
                        range <- cbind(qqq$mz-prems,qqq$mz+prems)
                        idxmz <- enviGCMS::getoverlapmass(range,matrix(c(prec-ppm/prec,prec+ppm/prec),nrow = 1))
                        if(sum(idxmz)>0){
                                msms <- qqq$msms[idxmz]
                                name <- qqq$name[idxmz]
                                mz2 <- qqq$mz[idxmz]
                                msmsraw <- orb$msmsraw[idxmz]
                                i <- function(x) sum(pmdt %in% unique(x))
                                result <- sapply(msms,i)
                                t <- list(name = name[result>length(pmdt)*pmdc],mz = mz2[result>length(pmdt)*pmdc],msms = msms[result>length(pmdt)*pmdc],msmsraw = msmsraw[result>length(pmdt)*pmdc],dmz = mz, dprc=prec, dins = ins/max(ins)*100, file = namemgf)
                                if(sum(result)==0){
                                        return(NULL)
                                }else{
                                        return(t)
                                }
                        }else{
                                return(NULL)
                        }
                }
        }else{
                return(NULL)
        }
}

#' Show MS/MS pmd annotation result from `pmdanno` function
#' @param anno list from pmdanno function
#' @return NULL
#' @examples
#' file <- system.file("extdata", "challenge-msms.mgf", package = "pmd")
#' anno <- pmdanno(file)
#' plotpmdanno(anno)
#' @export
plotpmdanno <- function(anno){
        for(i in 1:length(anno$msms)){
                graphics::plot(anno$msmsraw[[i]],type='h',main=anno$name[[i]])
                graphics::points(anno$dins~anno$dmz,type='h',col='red')
        }
}

#' read in MSP file as list for ms/ms annotation
#' @param file the path to your MSP file
#' @param digits mass or mass to charge ratio accuracy for pmd, default 2
#' @param icf intensity cutoff, default 10 percentage
#' @return list a list with MSP information for MS/MS annotation
#' @export
getms2pmd <- function(file,digits=2, icf=10){
        # this part is modified from compMS2Miner's code: https://github.com/WMBEdmands/compMS2Miner/blob/ee20d3d632b11729d6bbb5b5b93cd468b097251d/R/metID.matchSpectralDB.R
        msp <- readLines(file)
        # remove empty lines
        msp <- msp[msp != '']
        ncomp <- grep('^NAME:', msp, ignore.case = TRUE)
        splitFactorTmp <- rep(1:length(ncomp), diff(c(ncomp, length(msp) + 1)))

        li <- split(msp,f = splitFactorTmp)
        getmsp <- function(x){
                namet <- x[grep('^NAME:',x, ignore.case=TRUE)]
                name <- gsub('^NAME: ','',namet, ignore.case=TRUE)
                prect <- x[grep('^PRECURSORMZ: |^PRECURSOR M/Z: |^PRECURSOR MZ: |^PEPMASS: ',x, ignore.case=TRUE)]
                prec <- as.numeric(gsub('^PRECURSORMZ: |^PRECURSOR M/Z: |^PRECURSOR MZ: |^PEPMASS: ','',prect, ignore.case=TRUE))
                npt <- x[grep('^Num Peaks: ',x, ignore.case=TRUE)]
                np <- gsub('^Num Peaks: ','',npt,ignore.case = TRUE)
                if(as.numeric(np)>0){
                        # matrix of masses and intensities
                        massIntIndx <- which(grepl('^[0-9]', x) & !grepl(': ', x))
                        massesInts <- unlist(strsplit(x[massIntIndx], '\t| '))
                        massesInts <- as.numeric(massesInts[grep('^[0-9].*[0-9]$|^[0-9]$', massesInts)])
                        # if any NAs remove from indx
                        mz <-  massesInts[seq(1, length(massesInts), 2)]
                        ins <-  massesInts[seq(2, length(massesInts), 2)]
                        ins <- ins/max(ins)*100
                        msms <- cbind.data.frame(mz=mz,ins=ins)
                        msms <- msms[msms$ins>icf]
                        dis <- stats::dist(msms$mz, method = "manhattan")
                        diff <- round(as.numeric(dis), digits = digits)
                        diff <- diff[order(diff)]
                        return(list(name=name,prec=prec,msms=msms,pmd=diff))
                }else{
                        return(list(name=name,prec=prec,msms = NULL,pmd=diff))
                }

        }
        li <- lapply(li,getmsp)
        name <- sapply(li,function(x) x$name)
        mz <- sapply(li,function(x) x$prec)
        msms <- lapply(li, function(x) x$pmd)
        msmsraw <- lapply(li, function(x) x$msms)
        return(list(name = name, mz = mz, msms = msms, msmsraw = msmsraw))
}

#' read in MSP file as list for EI-MS annotation
#' @param file the path to your MSP file
#' @param digits mass or mass to charge ratio accuracy for pmd, default 0
#' @param icf intensity cutoff, default 10 percentage
#' @return list a list with MSP information for EI-MS annotation
#' @export
getmspmd <- function(file,digits=2){
        # this part is modified from compMS2Miner's code: https://github.com/WMBEdmands/compMS2Miner/blob/ee20d3d632b11729d6bbb5b5b93cd468b097251d/R/metID.matchSpectralDB.R
        msp <- readLines(file)
        # remove empty lines
        msp <- msp[msp != '']
        ncomp <- grep('^NAME:', msp, ignore.case = TRUE)
        splitFactorTmp <- rep(1:length(ncomp), diff(c(ncomp, length(msp) + 1)))

        li <- split(msp,f = splitFactorTmp)
        getmsp <- function(x){
                namet <- x[grep('^NAME:',x, ignore.case=TRUE)]
                name <- gsub('^NAME: ','',namet, ignore.case=TRUE)
                npt <- x[grep('^Num Peaks: ',x, ignore.case=TRUE)]
                np <- gsub('^Num Peaks: ','',npt,ignore.case = TRUE)
                if(as.numeric(np)>0){
                        # matrix of masses and intensities
                        massIntIndx <- which(grepl('^[0-9]', x) & !grepl(': ', x))
                        massesInts <- unlist(strsplit(x[massIntIndx], '\t| '))
                        massesInts <- as.numeric(massesInts[grep('^[0-9].*[0-9]$|^[0-9]$', massesInts)])
                        # if any NAs remove from indx
                        mz <-  massesInts[seq(1, length(massesInts), 2)]
                        ins <-  massesInts[seq(2, length(massesInts), 2)]
                        ins <- ins/max(ins)*100
                        msms <- cbind.data.frame(mz=mz,ins=ins)
                        msms <- msms[msms$ins>icf]
                        dis <- stats::dist(msms$mz, method = "manhattan")
                        diff <- round(as.numeric(dis), digits = digits)
                        diff <- diff[order(diff)]
                        return(list(name=name,msms=msms,pmd=diff))
                }else{
                        return(list(name=name,msms = NULL,pmd=diff))
                }
        }
        li <- lapply(li,getmsp)
        name <- sapply(li,function(x) x$name)
        msms <- lapply(li, function(x) x$pmd)
        msmsraw <- lapply(li, function(x) x$msms)
        return(list(name = name, msms = msms, msmsraw = msmsraw))
}
