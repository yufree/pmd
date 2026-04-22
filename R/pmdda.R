
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
        rtcluster <- .rt_clusters(rt, drt)
        inji <- rtcluster
        maxd <- max(table(rtcluster))
        m <- length(unique(rtcluster))
        inj <- ceiling(maxd / n)
        message(paste('You need', inj, 'injections!'))
        for (i in seq_len(m)) {
                z <- 1:inj
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
                                              rt %in% x), replace = TRUE)
        }
        return(inji)
}

#' Link pos mode peak list with neg mode peak list by pmd.
#' @param pos a list with mzrt profile collected from positive mode.
#' @param neg a list with mzrt profile collected from negative mode.
#' @param pmd numeric or numeric vector
#' @param digits mass or mass to charge ratio accuracy for pmd, default 2
#' @return dataframe with filtered positive and negative peak list
#' @export
getposneg <- function(pos, neg, pmd = 2.02, digits = 2) {
        np <- length(pos$mz)
        nn <- length(neg$mz)
        if (np == 0L || nn == 0L) return(NULL)

        # All pairwise m/z differences, rounded; find pairs matching any pmd.
        d  <- outer(pos$mz, neg$mz, `-`)
        rd <- round(d, digits)
        hits <- which(rd %in% pmd)                       # linear indices
        if (length(hits) == 0L) return(NULL)

        ip <- ((hits - 1L) %% np) + 1L                   # pos row index
        in_ <- ((hits - 1L) %/% np) + 1L                 # neg row index

        # Row-wise correlation between pos$data[ip,] and neg$data[in_,] via
        # z-score trick: cor = (1/(k-1)) * sum(zA * zB).
        zscore <- function(M) {
                M  <- as.matrix(M)
                mu <- rowMeans(M)
                sdv <- sqrt(rowSums((M - mu)^2) / pmax(ncol(M) - 1L, 1L))
                sdv[sdv == 0] <- NA_real_                # guard constant rows
                (M - mu) / sdv
        }
        Zp <- suppressWarnings(zscore(pos$data))
        Zn <- suppressWarnings(zscore(neg$data))
        k  <- ncol(Zp)
        if (ncol(Zn) != k)
                stop("getposneg: pos$data and neg$data must have the same number of samples.",
                     call. = FALSE)
        cor_vec <- rowSums(Zp[ip, , drop = FALSE] *
                                   Zn[in_, , drop = FALSE]) / max(k - 1L, 1L)

        data.frame(
                pos    = pos$mz[ip],
                rt     = pos$rt[ip],
                neg    = neg$mz[in_],
                rt.1   = neg$rt[in_],
                diffmz = d[hits],
                diffrt = pos$rt[ip] - neg$rt[in_],
                cor    = cor_vec,
                check.names = FALSE,
                stringsAsFactors = FALSE
        )
}

