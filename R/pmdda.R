
#' Get multiple injections index for selected retention time
#' @param rt retention time vector for peaks in seconds
#' @param drt retention time drift for targeted analysis in seconds, default 10.
#' @param n max ions numbers within retention time drift windows
#' @return index for each injection
#' @examples
#' data(spmeinvivo)
#' pmd <- getpaired(spmeinvivo)
#' std <- getstd(pmd)
#' index <- gettarget(std$rt[std$stdmassindex])
#' table(index)
#' @export
gettarget <- function(rt, drt = 10, n = 6) {
        dis <- stats::dist(rt, method = "manhattan")
        fit <- stats::hclust(dis)
        inji <- rtcluster <- stats::cutree(fit, h = drt)
        maxd <- max(table(rtcluster))
        m <- length(unique(rtcluster))
        inj <- ceiling(maxd / n)
        message(paste('You need', inj, 'injections!'))
        for (i in c(1:m)) {
                z = 1:inj
                x <- rt[rtcluster == i]
                while (length(x) > inj & length(x) > n) {
                        t <- sample(x, n)
                        w <- sample(z, 1)
                        inji[rt %in% t] <- w
                        z <- z[!(z %in% w)]
                        x <- x[!(x %in% t)]
                }
                inji[rtcluster == i &
                             rt %in% x] <-
                        sample(z, sum(rtcluster == i &
                                              rt %in% x), replace = T)
        }
        return(inji)
}

#' Link pos mode peak list with neg mode peak list by pmd.
#' @param pos a list with mzrt profile collected from positive mode.
#' @param neg a list with mzrt profile collected from negative mode.
#' @param pmd numeric or numeric vector
#' @param digits mass or mass to charge ratio accuracy for pmd, default 2
#' @return dataframe with filtered postive and negative peak list
#' @export
getposneg <- function(pos,neg, pmd = 2.02, digits = 2){
        df <- NULL
        x <- rep(NA,length(pos$mz))
        for(i in 1:length(pos$mz)){
                if(sum(round((pos$mz[i]-neg$mz),digits) %in% pmd) != 0){
                        index <- round((pos$mz[i]-neg$mz),digits) %in% pmd
                        if(sum(index)>1){
                                cor <- apply(neg$data[index,],1,function(x) suppressWarnings(cor(as.numeric(x),as.numeric(pos$data[i,]))))
                        }else{
                                cor <- suppressWarnings(cor(as.numeric(pos$data[i,]),as.numeric(neg$data[index,])))
                        }

                        t <- cbind.data.frame(pos=pos$mz[i],rt = pos$rt[i],neg=neg$mz[index],rt=neg$rt[index],diffmz=pos$mz[i]-neg$mz[index],diffrt=pos$rt[i]-neg$rt[index],cor=cor)
                        df <- rbind.data.frame(df,t)
                }
        }
        return(df)
}

