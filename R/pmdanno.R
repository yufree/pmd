#' Perform MS/MS pmd annotation for mgf file
#' @param file mgf file generated from MS/MS data
#' @param db database could be 'qtof', 'qqq', 'orb' data base from hmdb
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
