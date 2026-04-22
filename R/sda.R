# Suppress R CMD check NOTEs for data.table non-standard-evaluation column
# references inside getchainseq (and its helpers).  These are not real R
# variables at parse time, so the codetools walker flags them.
utils::globalVariables(c(
        ".", ".SD", ".end", ".len",
        ":=", "cor", "cor_mean",
        "from", "to",
        "path_len", "pmd", "pmd_r", "rtd"
))

.datatable.aware <- TRUE

# Build pairwise PMD data.table from mz (and optional rt, rtg, data).
# Returns a data.table with ms1, ms2, diff, and optional rt1/rt2/diffrt,
# rtg1/rtg2/rtgdiff, cor, md, diff2 columns.
# Uses direct index generation (not stats::dist + lower.tri) to avoid the
# O(n^2) matrix allocation when n is large.
# Order of pairs matches the lower triangle of an n x n matrix in
# column-major order (same as stats::dist), so ms1 comes from the later
# position and ms2 from the earlier position.
.pair_df <- function(mz, rt = NULL, rtg = NULL, data = NULL,
                     digits = NULL, include_md = FALSE) {
        n <- length(mz)
        if (n < 2L) {
                return(data.table::data.table(ms1 = numeric(0), ms2 = numeric(0),
                                              diff = numeric(0)))
        }
        # column-major lower-triangular indices (i > j)
        j <- rep.int(seq_len(n - 1L), (n - 1L):1L)
        i <- sequence.default((n - 1L):1L, from = seq.int(2L, n))
        diff <- abs(mz[i] - mz[j])
        cols <- list(ms1 = mz[i], ms2 = mz[j], diff = diff)
        if (!is.null(rt)) {
                cols$rt1 <- rt[i]; cols$rt2 <- rt[j]
                cols$diffrt <- abs(rt[i] - rt[j])
        }
        if (!is.null(rtg)) {
                cols$rtg1 <- rtg[i]; cols$rtg2 <- rtg[j]
                cols$rtgdiff <- abs(rtg[i] - rtg[j])
        }
        if (!is.null(data)) {
                cormat <- stats::cor(t(data))
                cols$cor <- cormat[cbind(i, j)]
        }
        if (include_md)      cols$md    <- diff %% 1
        if (!is.null(digits)) cols$diff2 <- round(diff, digits)
        data.table::as.data.table(cols)
}

#' Perform structure/reaction directed analysis for peaks list.
#' @param list a list with mzrt profile
#' @param rtcutoff cutoff of the distances in retention time hierarchical clustering analysis, default 10
#' @param corcutoff cutoff of the correlation coefficient, default NULL
#' @param digits mass or mass to charge ratio accuracy for pmd, default 2
#' @param accuracy measured mass or mass to charge ratio in digits, default 4
#' @param freqcutoff pmd frequency cutoff for structures or reactions, default NULL. This cutoff will be found by PMD network analysis when it is NULL.
#' @return list with tentative isotope, adducts, and neutral loss peaks' index, retention time clusters.
#' @examples
#' data(spmeinvivo)
#' pmd <- getpaired(spmeinvivo)
#' std <- getstd(pmd)
#' sda <- getsda(std)
#' @seealso \code{\link{getpaired}},\code{\link{getstd}},\code{\link{plotpaired}}
#' @export
getsda <-
        function(list,
                 rtcutoff = 10,
                 corcutoff = NULL,
                 digits = 2,
                 accuracy = 4,
                 freqcutoff = NULL) {
                .check_mzrt(list, "getsda")
                if (!is.null(list$stdmass) && is.null(list$stdmassindex))
                        stop("getsda: `stdmass` present but `stdmassindex` missing; ",
                             "did you skip getstd()?", call. = FALSE)
                if (!is.null(list$paired) && is.null(list$pairedindex))
                        stop("getsda: `paired` present but `pairedindex` missing; ",
                             "did you skip getpaired()?", call. = FALSE)
                if (is.null(list$stdmass) & is.null(list$paired)) {
                        mz <- list$mz
                        rt <- list$rt
                        data <- list$data
                        rtg <- .rt_clusters(rt, rtcutoff)
                } else if (is.null(list$stdmass)) {
                        mz <- list$mz[list$pairedindex]
                        rt <- list$rt[list$pairedindex]
                        data <- list$data[list$pairedindex, , drop = FALSE]
                        rtg <- list$rtcluster[list$pairedindex]
                } else {
                        mz <- list$mz[list$stdmassindex]
                        rt <- list$rt[list$stdmassindex]
                        data <- list$data[list$stdmassindex, , drop = FALSE]
                        rtg <- list$rtcluster[list$stdmassindex]
                }
                # PMD analysis (pairwise table)
                dt <- .pair_df(mz, rt = rt, rtg = rtg, data = data,
                               digits = digits, include_md = TRUE)
                dt <- dt[dt$rtgdiff > 0, ]
                df <- as.data.frame(dt)
                # use unique isomers
                index <-
                        !duplicated(paste0(round(df$ms1, accuracy),
                                           round(df$ms2, accuracy)))
                diff <- df$diff2[index]
                freq <-
                        sort(table(diff), decreasing = TRUE)
                if (is.null(freqcutoff)) {
                        dis <- c()
                        for (i in seq_len(min(length(freq), 100))) {
                                pmd <- as.numeric(names(freq))[1:i]
                                dfx <- df[df$diff2 %in% pmd, c(1, 2)]
                                net <-
                                        igraph::graph_from_data_frame(dfx, directed = FALSE)
                                dis[i] <- igraph::mean_distance(net)
                        }
                        n <- which.max(dis)
                        freqt <- freq[n - 1]
                        if (n == 100) {
                                warning(
                                        "Average distance is still increasing, you need to check manually for frequency cutoff."
                                )
                        }
                        if (sum(df$diff2 == 0) > freqt &
                            0 %in% as.numeric(names(freq))) {
                                list$sda <-
                                        df[df$diff2 %in% c(0, as.numeric(names(freq[freq > freqt]))), ]
                        } else{
                                list$sda <- df[df$diff2 %in% as.numeric(names(freq[freq > freqt])), ]
                        }
                        message(
                                paste(
                                        "PMD frequency cutoff is",
                                        freqt,
                                        'by PMD network analysis with largest network average distance',
                                        round(max(dis), 2),
                                        '.'
                                )
                        )
                        # i <- dis <- t <- 1
                        # while(n>=t){
                        #         pmd <- as.numeric(names(freq[freq>i]))
                        #         dfx <- df[df$diff2 %in% pmd,c(1,2)]
                        #         net <- igraph::graph_from_data_frame(dfx,directed = FALSE)
                        #         t <- n
                        #         n <- length(igraph::groups(igraph::components(net)))
                        #         i <- i+1
                        # }
                        # message(paste("PMD frequency cutoff is", i-1, 'by PMD network analysis with',t,'clusters.'))
                        # if(sum(df$diff2 == 0)>(i-1) & 0 %in% as.numeric(names(freq))){
                        #         list$sda <- df[df$diff2 %in% c(0,as.numeric(names(freq[freq>(i-1)]))),]
                        # }else{
                        #         list$sda <- df[df$diff2 %in% as.numeric(names(freq[freq>(i-1)])),]
                        # }
                } else{
                        if (sum(df$diff2 == 0) > freqcutoff &
                            0 %in% as.numeric(names(freq))) {
                                list$sda <- df[(df$diff2 %in% c(0, as.numeric(names(
                                        freq[freq >=
                                                     freqcutoff]
                                )))), ]
                        } else{
                                list$sda <- df[(df$diff2 %in% c(as.numeric(names(
                                        freq[freq >=
                                                     freqcutoff]
                                )))), ]
                        }
                }

                if (!is.null(corcutoff) & !is.null(data)) {
                        list$sda <- list$sda[abs(list$sda$cor) >= corcutoff,]
                }
                # show message about std mass
                sub <- names(table(list$sda$diff2))
                n <- length(sub)
                message(paste(n, "groups were found as high frequency PMD group."))
                message(paste(sub, "was found as high frequency PMD.",
                              "\n"))
                return(list)
        }
#' Perform structure/reaction directed analysis for mass only.
#' @param mz numeric vector for independent mass or mass to charge ratio. Mass to charge ratio from GlobalStd algorithm is suggested. Isomers would be excluded automated
#' @param pmd a specific paired mass distance or a vector of pmds, default NULL
#' @param freqcutoff pmd frequency cutoff for structures or reactions, default 10
#' @param digits mass or mass to charge ratio accuracy for pmd, default 3
#' @param top top n pmd frequency cutoff when the freqcutoff is too small for large data set
#' @param formula vector for formula when you don't have mass or mass to charge ratio data
#' @param mdrange mass defect range to ignore. Default c(0.25,0.9) to retain the possible reaction related paired mass
#' @param verbose logic, if TURE, return will be llist with paired mass distances table. Default FALSE.
#' @return logical matrix with row as the same order of mz or formula and column as high  frequency pmd group when verbose is FALSE
#' @examples
#' data(spmeinvivo)
#' pmd <- getpaired(spmeinvivo)
#' std <- getstd(pmd)
#' sda <- getrda(spmeinvivo$mz[std$stdmassindex])
#' sda <- getrda(spmeinvivo$mz, pmd = c(2.016,15.995,18.011,14.016))
#' @seealso \code{\link{getsda}}
#' @export
getrda <-
        function(mz,
                 pmd = NULL,
                 freqcutoff = 10,
                 digits = 3,
                 top = 20,
                 formula = NULL,
                 mdrange = c(0.25,0.9),
                 verbose = FALSE) {
                if (!is.null(formula)) {
                        mz <- unlist(Map(enviGCMS::getmass, formula))
                }
                mz <- unique(mz)
                df <- as.data.frame(.pair_df(mz, digits = digits, include_md = TRUE))
                if(is.null(pmd[1])){
                        if(!is.null(mdrange)){
                                df <- df[df$md<mdrange[1]|df$md>mdrange[2],]
                        }
                        freq <-
                                sort(table(df$diff2), decreasing = TRUE)
                        message(paste(length(freq), 'pmd found.'))
                        if (!is.null(top)) {
                                freq <- utils::head(freq, top)
                        }
                        sda <-
                                df[(df$diff2 %in% c(as.numeric(names(freq[freq >= freqcutoff])))), ]
                } else{
                        sda <- df[(df$diff2 %in% round(pmd,digits = digits)), ]
                }
                pmd <- unique(sda$diff2)[order(unique(sda$diff2))]
                message(paste(length(pmd), 'pmd used.'))
                df <- NULL

                split <- split.data.frame(sda, sda$diff2)
                rtpmd <- function(bin, i) {
                        mass <- unique(c(bin$ms1[bin$diff2 == i], bin$ms2[bin$diff2 == i]))
                        index <- mz %in% mass
                        return(index)
                }
                result <-
                        mapply(rtpmd, split, as.numeric(names(split)))
                rownames(result) <- mz

                if(verbose){
                        return(list(mdt = as.data.frame(sda),result=result))
                }else{
                        return(as.data.frame(result))
                }
        }
#' Perform correlation directed analysis for peaks list.
#' @param list a list with mzrt profile
#' @param rtcutoff cutoff of the distances in retention time hierarchical clustering analysis, default 10
#' @param corcutoff cutoff of the correlation coefficient, default NULL
#' @param accuracy measured mass or mass to charge ratio in digits, default 4
#' @return list with correlation directed analysis results
#' @examples
#' data(spmeinvivo)
#' cluster <- getpseudospectrum(spmeinvivo)
#' cbp <- enviGCMS::getfilter(cluster,rowindex = cluster$stdmassindex2)
#' cda <- getcda(cbp)
#' @seealso \code{\link{getsda}},\code{\link{getrda}}
#' @export
getcda <- function(list,
                   corcutoff = 0.9,
                   rtcutoff = 10,
                   accuracy = 4) {
        .check_mzrt(list, "getcda", need_data = TRUE)
        mz <- list$mz
        rt <- list$rt
        data <- list$data
        rtg <- .rt_clusters(rt, rtcutoff)
        dt <- .pair_df(mz, rt = rt, rtg = rtg, data = data, include_md = TRUE)
        list$cda <- as.data.frame(dt[abs(dt$cor) >= corcutoff, ])
        return(list)
}

#' Get pmd for specific reaction
#' @param list a list with mzrt profile
#' @param pmd a specific paired mass distance or a vector of pmds
#' @param rtcutoff cutoff of the distances in retention time hierarchical clustering analysis, default 10
#' @param corcutoff cutoff of the correlation coefficient, default NULL
#' @param digits mass or mass to charge ratio accuracy for pmd, default 2
#' @param accuracy measured mass or mass to charge ratio in digits, default 4
#' @return list with paired peaks for specific pmd or pmds.
#' @examples
#' data(spmeinvivo)
#' pmd <- getpmd(spmeinvivo,pmd=15.99)
#' @seealso \code{\link{getpaired}},\code{\link{getstd}},\code{\link{getsda}},\code{\link{getrda}}
#' @export
getpmd <-
        function(list,
                 pmd,
                 rtcutoff = 10,
                 corcutoff = NULL,
                 digits = 2,
                 accuracy = 4) {
                mz <- list$mz
                data <- list$data
                has_rt <- !is.null(list$rt)
                rt <- if (has_rt) list$rt else NULL
                rtg <- if (has_rt) .rt_clusters(rt, rtcutoff) else NULL

                dt <- .pair_df(mz, rt = rt, rtg = rtg, data = data,
                               digits = digits)

                if (!is.null(corcutoff)) dt <- dt[abs(dt$cor) >= corcutoff, ]
                dt <- if (has_rt) dt[dt$rtgdiff > 0 & dt$diff2 %in% pmd, ]
                      else            dt[dt$diff2 %in% pmd, ]

                df <- as.data.frame(dt)
                list$pmd <- df
                # high (hi) = larger mz in the pair; low (lo) = smaller mz
                hi <- pmax(df$ms1, df$ms2)
                lo <- pmin(df$ms1, df$ms2)
                if (has_rt) {
                        swap <- df$ms1 <= df$ms2
                        rtg_hi <- ifelse(swap, df$rtg2, df$rtg1)
                        rtg_lo <- ifelse(swap, df$rtg1, df$rtg2)
                        indexh <- unique(paste(round(hi, accuracy), rtg_hi))
                        indexl <- unique(paste(round(lo, accuracy), rtg_lo))
                        index0 <- paste(round(list$mz, accuracy), rtg)
                } else {
                        indexh <- unique(round(hi, accuracy))
                        indexl <- unique(round(lo, accuracy))
                        index0 <- round(list$mz, accuracy)
                }
                list$pmdindex  <- index0 %in% unique(c(indexh, indexl))
                list$pmdindexh <- index0 %in% indexh
                list$pmdindexl <- index0 %in% indexl
                return(list)
        }

#' Get pmd details for specific reaction after the removal of isotopouge.
#' @param mz a vector of mass to charge ratio.
#' @param group mass to charge ratio group from either retention time or mass spectrometry imaging segmentation.

#' @param pmd a specific paired mass distance or a vector of pmds
#' @param digits mass or mass to charge ratio accuracy for pmd, default 2.
#' @param mdrange mass defect range to ignore. Default c(0.25,0.9) to retain the possible reaction related paired mass.
#' @return dataframe with paired peaks for specific pmd or pmds. When group is provided, a column named net will be generated to show if certain pmd will be local(within the same group) or global(across the groups)
#' @examples
#' data(spmeinvivo)
#' pmddf <- getpmddf(spmeinvivo$mz,pmd=15.99)
#' @seealso \code{\link{getpaired}},\code{\link{getstd}},\code{\link{getsda}},\code{\link{getrda}}
#' @export
getpmddf <- function(mz,group=NULL,pmd=NULL,digits=2,mdrange=c(0.25,0.9)){
        df <- as.data.frame(.pair_df(mz, digits = digits, include_md = TRUE))
        # getpmddf traditionally returns ms1 <= ms2 (pmin/pmax)
        swap <- df$ms1 > df$ms2
        if (any(swap)) {
                tmp <- df$ms1[swap]; df$ms1[swap] <- df$ms2[swap]; df$ms2[swap] <- tmp
        }
        if(!is.null(group)){
                df$group1 <- group[match(df$ms1,mz)]
                df$group2 <- group[match(df$ms2,mz)]
        }else{
                df$group1 <- df$group2 <- rep(1,length(df$md))
        }

        idx <- df$md<mdrange[1]|df$md>mdrange[2]
        df <- df[idx,]
        isoindex <- (round(df$diff, digits) != 0) & ((
                df$diff %% 1 < 0.01 &
                        df$diff >= 1 &
                        df$diff < 2
        ) | (
                df$diff %% 2 < 0.01 &
                        df$diff >= 2 &
                        df$diff < 3
        ) | (
                df$diff %% 1 > 0.99 &
                        df$diff >= 1 &
                        df$diff < 2
        ) | (
                df$diff %% 1 > 0.99 &
                        df$diff >= 0 &
                        df$diff < 1
        )
        )
        dfiso <- df[isoindex,]
        dfdeiso <- df[!isoindex,]

        if(!is.null(pmd)){
                dfdeisopmd <- dfdeiso[dfdeiso$diff2 %in% pmd,]
                dfdeisopmd$net <- ifelse(abs(dfdeisopmd$group1-dfdeisopmd$group2)==0,'local','global')
                return(dfdeisopmd)
        }else{
                dfdeiso$net <- ifelse(abs(dfdeiso$group1-dfdeiso$group2)==0,'local','global')
                return(dfdeiso)
        }
}

#' Get reaction chain for specific mass to charge ratio
#' @param list a list with mzrt profile
#' @param diff paired mass distance(s) of interests
#' @param mass a specific mass for known compound or a vector of masses. You could also input formula for certain compounds
#' @param digits mass or mass to charge ratio accuracy for pmd, default 2
#' @param accuracy measured mass or mass to charge ratio in digits, default 4
#' @param rtcutoff cutoff of the distances in retention time hierarchical clustering analysis, default 10
#' @param corcutoff cutoff of the correlation coefficient, default 0.6
#' @param ppm all the peaks within this mass accuracy as seed mass or formula
#' @return a list with mzrt profile and reaction chain dataframe
#' @examples
#' data(spmeinvivo)
#' # check metabolites of C18H39NO
#' pmd <- getchain(spmeinvivo,diff = c(2.02,14.02,15.99),mass = 286.3101)
#' # remove the retention time for mass only data
#' spmeinvivo$rt <- NULL
#' pmd <- getchain(spmeinvivo,diff = c(2.02,14.02,15.99),mass = 286.3101)
#' @export
getchain <-
        function(list,
                 diff,
                 mass,
                 digits = 2,
                 accuracy = 4,
                 rtcutoff = 10,
                 corcutoff = 0.6,
                 ppm = 25) {
                if (is.character(mass)) {
                        mass <- unlist(Map(enviGCMS::getmass, mass))
                }
                massup <- mass + mass * ppm / 1e6
                massdown <- mass - mass * ppm / 1e6
                updown <- vapply(Map(function(x)
                        x < massup & x > massdown, list$mz), function(x)
                                sum(x & TRUE) > 0, TRUE)
                mass <- list$mz[updown]
                mass <- unique(round(mass, accuracy))

                mz <- list$mz
                data <- list$data
                has_rt <- !is.null(list$rt)
                rt <- if (has_rt) list$rt else NULL
                rtg <- if (has_rt) .rt_clusters(rt, rtcutoff) else NULL

                dt <- .pair_df(mz, rt = rt, rtg = rtg, data = data,
                               digits = digits)
                keep <- dt$diff2 %in% diff
                if (has_rt) keep <- keep & dt$rtgdiff > 0
                dt <- dt[keep, ]
                if (!is.null(corcutoff)) dt <- dt[abs(dt$cor) >= corcutoff, ]
                df <- as.data.frame(dt)

                seed <- NULL
                ms1 <- round(df$ms1, digits = accuracy)
                ms2 <- round(df$ms2, digits = accuracy)
                if (length(mass) == 1) {
                        mass <- round(mass, accuracy)
                        sdat <-
                                unique(c(mass, ms2[ms1 %in% mass], ms1[ms2 %in% mass]))
                        while (!identical(sdat, seed)) {
                                seed <- sdat
                                sdat <-
                                        unique(c(sdat, ms2[ms1 %in% sdat], ms1[ms2 %in% sdat]))
                        }
                        list$sdac <- df[ms1 %in% sdat | ms2 %in% sdat , ]
                        return(list)
                } else if (length(mass) == 0) {
                        warning(
                                'No mass input and all mass in the list will be used for reaction chain construction!'
                        )
                        sdac <- NULL
                        mass <- round(list$mz, accuracy)
                        for (i in seq_along(mass)) {
                                sdat <-
                                        unique(c(mass[i], ms2[ms1 %in% mass[i]], ms1[ms2 %in% mass[i]]))
                                if (length(sdat) != 1) {
                                        while (!identical(sdat, seed)) {
                                                seed <- sdat
                                                sdat <-
                                                        unique(c(sdat, ms2[ms1 %in% sdat], ms1[ms2 %in% sdat]))
                                        }
                                        sdact <-
                                                df[ms1 %in% sdat |
                                                           ms2 %in% sdat , ]
                                        sdact$mass <- mass[i]
                                        sdac <-
                                                rbind.data.frame(sdac, sdact)
                                }
                        }
                        list$sdac <- sdac[!duplicated(sdac),]
                        return(list)
                } else{
                        sdac <- NULL
                        mass <- round(mass, accuracy)
                        for (i in seq_along(mass)) {
                                sdat <-
                                        unique(c(mass[i], ms2[ms1 %in% mass[i]], ms1[ms2 %in% mass[i]]))
                                if (length(sdat) != 1) {
                                        while (!identical(sdat, seed)) {
                                                seed <- sdat
                                                sdat <-
                                                        unique(c(sdat, ms2[ms1 %in% sdat], ms1[ms2 %in% sdat]))
                                        }
                                        sdact <-
                                                df[ms1 %in% sdat | ms2 %in% sdat, ]
                                        sdact$mass <- mass[i]
                                        sdac <-
                                                rbind.data.frame(sdac, sdact)
                                }
                        }
                        list$sdac <- sdac[!duplicated(sdac),]
                        return(list)
                }

        }

#' Get quantitative paired peaks list for specific reaction/pmd
#' @param list a list with mzrt profile and data
#' @param pmd a specific paired mass distances
#' @param rtcutoff cutoff of the distances in retention time hierarchical clustering analysis, default 10
#' @param digits mass or mass to charge ratio accuracy for pmd, default 2
#' @param accuracy measured mass or mass to charge ratio in digits, default 4
#' @param cvcutoff ratio or intensity cv cutoff for quantitative paired peaks, default 30
#' @param method quantification method can be 'static' or 'dynamic'. See details.
#' @param outlier logical, if true, outlier of ratio will be removed, default False.
#' @param ... other parameters for getpmd
#' @details PMD based reaction quantification methods have two options: 'static' will only consider the stable mass pairs across samples and such reactions will be limited by the enzyme or other factors than substrates. 'dynamic' will consider the unstable paired masses by normalization the relatively unstable peak with stable peak between paired masses and such reactions will be limited by one or both peaks in the paired masses.
#' @return list with quantitative paired peaks.
#' @examples
#' data(spmeinvivo)
#' pmd <- getreact(spmeinvivo,pmd=15.99)
#' @seealso \code{\link{getpaired}},\code{\link{getstd}},\code{\link{getsda}},\code{\link{getrda}},\code{\link{getpmd}},
#' @export
getreact <-
        function(list,
                 pmd,
                 rtcutoff = 10,
                 digits = 2,
                 accuracy = 4,
                 cvcutoff = 30,
                 outlier = FALSE,
                 method = 'static',
                 ...) {
                p <-
                        pmd::getpmd(
                                list,
                                pmd = pmd,
                                rtcutoff = rtcutoff,
                                digits = digits,
                                accuracy = accuracy,
                                ...
                        )
                # Vectorised row-wise %RSD = sd/mean*100 (NA safe)
                row_rsd <- function(m) {
                        mu <- rowMeans(m, na.rm = TRUE)
                        sdv <- sqrt(rowSums((m - mu)^2, na.rm = TRUE) /
                                    pmax(rowSums(!is.na(m)) - 1L, 1L))
                        sdv / mu * 100
                }
                if (sum(p$pmdindex) > 0) {
                        list <- enviGCMS::getfilter(p, p$pmdindex)
                        data <- list$data
                        pmd <- list$pmd
                        has_rt <- !is.null(list$rt)
                        # Build lookup keys: one per feature, one per pmd pair endpoint
                        key   <- if (has_rt) paste(list$mz, list$rt) else as.character(list$mz)
                        keys1 <- if (has_rt) paste(pmd$ms1, pmd$rt1) else as.character(pmd$ms1)
                        keys2 <- if (has_rt) paste(pmd$ms2, pmd$rt2) else as.character(pmd$ms2)
                        idx1 <- match(keys1, key)
                        idx2 <- match(keys2, key)

                        rr1 <- data[idx1, , drop = FALSE]
                        rr2 <- data[idx2, , drop = FALSE]
                        ratios <- rr1 / rr2

                        if (outlier) {
                                # Per-row outlier removal via boxplot.stats
                                r_vec <- vapply(seq_len(nrow(ratios)), function(k) {
                                        rv <- ratios[k, ]
                                        out <- grDevices::boxplot.stats(rv)$out
                                        rv <- rv[!rv %in% out]
                                        stats::sd(rv, na.rm = TRUE) / mean(rv, na.rm = TRUE) * 100
                                }, numeric(1))
                                list$pmd$r <- r_vec
                        } else {
                                list$pmd$r <- row_rsd(ratios)
                        }
                        list$pmd$rh <- row_rsd(rr1)
                        list$pmd$rl <- row_rsd(rr2)
                        list$pmdindex <- list$pmdindexh <- list$pmdindexl <- NULL
                        if (method == 'static') {
                                list$pmd <- list$pmd[list$pmd$r < cvcutoff & (list$pmd$rh>cvcutoff | list$pmd$rl>cvcutoff),]
                                list$pmd <-
                                        list$pmd[stats::complete.cases(list$pmd), ]
                                if (nrow(list$pmd) > 0&!is.null(list$rt)) {

                                        idx <- paste(list$mz, list$rt)
                                        idx2 <- unique(c(
                                                paste(list$pmd$ms1, list$pmd$rt1),
                                                paste(list$pmd$ms2, list$pmd$rt2)
                                        ))
                                        list <-
                                                enviGCMS::getfilter(list, idx %in% idx2)
                                        idx <-
                                                paste(list$mz, list$rt)
                                        pmdh <-
                                                list$data[match(paste(list$pmd$ms1, list$pmd$rt1),
                                                                idx), , drop = FALSE]
                                        pmdl <-
                                                list$data[match(paste(list$pmd$ms2, list$pmd$rt2),
                                                                idx), , drop = FALSE]
                                        list$pmddata <- pmdh + pmdl
                                        return(list)
                                } else if (nrow(list$pmd) > 0&is.null(list$rt)){
                                        mzl <- list$pmd$rh<cvcutoff & list$pmd$rl<cvcutoff
                                        idx <- unique(c(list$pmd$ms1, list$pmd$ms2))
                                        list <-
                                                enviGCMS::getfilter(list, list$mz %in% idx)
                                        pmdh <-
                                                list$data[match(list$pmd$ms1,list$mz), , drop = FALSE]
                                        pmdl <-
                                                list$data[match(list$pmd$ms2,list$mz), , drop = FALSE]
                                        list$pmddata <- pmdh + pmdl
                                        return(list)
                                } else {
                                        message('No static quantitative peaks could be used.')
                                }
                        } else if (method == 'dynamic') {
                                list$pmd <-
                                        list$pmd[!(list$pmd$r < cvcutoff)& (list$pmd$rh>cvcutoff | list$pmd$rl>cvcutoff),]
                                list$pmd <-
                                        list$pmd[stats::complete.cases(list$pmd), ]
                                if (nrow(list$pmd) > 0&!is.null(list$rt)) {
                                        idx <- paste(list$mz, list$rt)
                                        idx2 <- unique(c(
                                                paste(list$pmd$ms1, list$pmd$rt1),
                                                paste(list$pmd$ms2, list$pmd$rt2)
                                        ))
                                        list <-
                                                enviGCMS::getfilter(list, idx %in% idx2)
                                        idx <-
                                                paste(list$mz, list$rt)
                                        idy <- list$pmd$rh > list$pmd$rl
                                        pmddata <-
                                                as.data.frame(matrix(
                                                        nrow = nrow(list$pmd),
                                                        ncol = ncol(list$data)
                                                ))
                                        pmddata[idy, ] <-
                                                list$data[match(paste(
                                                        list$pmd$ms2[idy],
                                                        list$pmd$rt2[idy]
                                                ),
                                                idx), , drop = FALSE]/list$data[match(paste(
                                                        list$pmd$ms1[idy],
                                                        list$pmd$rt1[idy]
                                                ),
                                                idx), , drop = FALSE]

                                        pmddata[!idy, ] <-
                                                list$data[match(paste(
                                                        list$pmd$ms1[!idy],
                                                        list$pmd$rt1[!idy]
                                                ),
                                                idx), , drop = FALSE]/list$data[match(paste(
                                                        list$pmd$ms2[!idy],
                                                        list$pmd$rt2[!idy]
                                                ),
                                                idx), , drop = FALSE]
                                        list$pmddata <- pmddata
                                        colnames(list$pmddata) <-
                                                colnames(list$data)
                                        return(list)
                                } else if (nrow(list$pmd) > 0&is.null(list$rt)){
                                        idx <- unique(c(list$pmd$ms1, list$pmd$ms2))
                                        list <-
                                                enviGCMS::getfilter(list, list$mz %in% idx)
                                        idy <- list$pmd$rh > list$pmd$rl
                                        pmddata <-
                                                as.data.frame(matrix(
                                                        nrow = nrow(list$pmd),
                                                        ncol = ncol(list$data)
                                                ))
                                        pmddata[idy, ] <-
                                                list$data[match(
                                                        list$pmd$ms2[idy],
                                                        list$mz), , drop = FALSE]/list$data[match(
                                                                list$pmd$ms1[idy],
                                                                list$mz), , drop = FALSE]
                                        pmddata[!idy, ] <-
                                                list$data[match(list$pmd$ms1[!idy],
                                                                list$mz), , drop = FALSE]/list$data[match(list$pmd$ms2[!idy],
                                                                                            list$mz), , drop = FALSE]
                                        list$pmddata <- pmddata
                                        colnames(list$pmddata) <-
                                                colnames(list$data)
                                        return(list)

                                }
                                else{
                                        message(
                                                'No dynamic quantitative peak could be used.'
                                        )
                                }
                        }
                }else{
                        message(
                                'No pmd peaks could be found for quantitative analysis.'
                        )
                }
        }

#' Parse a PMD pattern string into a list of step specs
#'
#' Converts a compact string grammar into the list-of-steps format accepted by
#' \code{\link{getchainseq}}. Steps are separated by semicolons. Each step is a
#' numeric PMD (or \code{*} for wildcard) optionally followed by a quantifier:
#' \itemize{
#'   \item \code{+}  one or more (min=1, max=Inf)
#'   \item \code{*}  zero or more (min=0, max=Inf)
#'   \item \code{?}  zero or one  (min=0, max=1)
#'   \item \code{\{n\}}    exactly n
#'   \item \code{\{n,m\}}  n to m times
#'   \item \code{\{n,\}}   at least n times
#' }
#' Whitespace is ignored.
#'
#' @param s a single character string, e.g. \code{"162.0528; -18.0106{0,3}; 14.0157+"}
#' @return a list of step specs suitable for \code{getchainseq(pattern = ...)}
#' @examples
#' parse_pmd_pattern("162.0528; -18.0106{0,3}")
#' parse_pmd_pattern("*; 14.0157+")
#' @export
parse_pmd_pattern <- function(s) {
        if (!is.character(s) || length(s) != 1L)
                stop("s must be a single character string.")
        steps <- trimws(strsplit(s, ";", fixed = TRUE)[[1]])
        steps <- steps[nzchar(steps)]
        lapply(steps, function(tok) {
                tok <- gsub("\\s+", "", tok)
                m <- regmatches(tok, regexec(
                        "^(\\*|-?[0-9]+\\.?[0-9]*)(\\+|\\*|\\?|\\{[0-9]+,?[0-9]*\\})?$",
                        tok))[[1]]
                if (length(m) == 0L)
                        stop(sprintf("Cannot parse step: '%s'", tok))
                head <- m[2]; quant <- m[3]
                pmd_val <- if (head == "*") NA_real_ else as.numeric(head)

                if (!nzchar(quant)) {
                        mn <- 1L; mx <- 1L
                } else if (quant == "+") {
                        mn <- 1L; mx <- Inf
                } else if (quant == "*") {
                        mn <- 0L; mx <- Inf
                } else if (quant == "?") {
                        mn <- 0L; mx <- 1L
                } else {
                        inner <- substr(quant, 2L, nchar(quant) - 1L)
                        parts <- strsplit(inner, ",", fixed = TRUE)[[1]]
                        if (length(parts) == 1L) {
                                mn <- mx <- as.integer(parts[1])
                        } else if (length(parts) == 2L) {
                                mn <- as.integer(parts[1])
                                mx <- if (nzchar(parts[2])) as.integer(parts[2]) else Inf
                        } else {
                                stop(sprintf("Bad quantifier: '%s'", quant))
                        }
                }
                list(pmd = pmd_val, min = mn, max = mx)
        })
}


#' Get reaction pathway chains matching an ordered PMD pattern with quantifiers
#'
#' Searches a feature network for directed paths whose successive mass differences
#' follow a user-specified PMD pattern. Unlike \code{\link{getchain}}, which
#' extracts the connected component reachable by any PMD in a set, this function
#' matches an \emph{ordered} sequence of PMDs, with support for wildcards and
#' regex-style quantifiers. This unifies three use cases: fixed reaction
#' sequences (e.g. glycosylation followed by dehydration), homologous series
#' (repeated \code{+CH2}, PEG units, etc.), and paths with unknown intermediates
#' (wildcard steps).
#'
#' Pattern grammar (any of):
#' \itemize{
#'   \item Numeric vector: each element is one fixed-PMD step. Example:
#'     \code{c(162.0528, -18.0106)} = glycosylation followed by one dehydration.
#'   \item List of step specs: each step is a list with \code{pmd} (numeric, or
#'     \code{NA} for wildcard), \code{min} (default 1), \code{max}
#'     (default = \code{min}; \code{Inf} allowed).
#'   \item Character string parsed by \code{\link{parse_pmd_pattern}}.
#' }
#'
#' Signs of PMDs are significant: \code{+162} = mass gain, \code{-162} =
#' mass loss. This separates e.g. glycosylation from deglycosylation.
#'
#' @param list a pmd-style list with \code{mz}, \code{rt} (optional), \code{data}
#' @param pattern numeric vector, list of step specs, or DSL string
#'   (see \code{\link{parse_pmd_pattern}})
#' @param mass optional seed mass(es) or formula(s); only paths starting within
#'   \code{ppm} of a seed are returned
#' @param digits PMD matching precision, default 4
#' @param rtcutoff RT hierarchical-clustering cutoff for isomer grouping, default 10
#' @param corcutoff correlation cutoff between linked features, default 0.6;
#'   pass \code{NULL} to disable
#' @param ppm seed mass ppm tolerance, default 25
#' @param rtdir RT direction per edge: \code{"any"}, \code{"increasing"},
#'   \code{"decreasing"}, or a numeric vector (one entry per unrolled edge,
#'   values in \code{\{-1, 0, 1\}})
#' @param max_paths safety cap on returned paths, default 1e5
#' @param allow_cycles if \code{FALSE} (default), paths may not revisit a node
#' @return the input \code{list} with two new elements: \code{sdacseq}
#'   (a data.table of matched paths) and \code{pattern} (the normalized pattern).
#'   \code{sdacseq} columns: \code{n1..nK} (node indices; NA-padded for shorter
#'   paths in variable-length matches), \code{path_len}, \code{mz_1..mz_K},
#'   \code{rt_1..rt_K} (if RT present), \code{pmd_1..pmd_\{K-1\}} (observed
#'   PMDs), and \code{cor_mean}.
#' @examples
#' \dontrun{
#' data(spmeinvivo)
#'
#' # Fixed sequence: glycosylation then dehydration
#' r1 <- getchainseq(spmeinvivo, c(162.0528, -18.0106))
#'
#' # Homologous series: 2+ CH2 extensions
#' r2 <- getchainseq(spmeinvivo,
#'                   list(list(pmd = 14.0157, min = 2, max = Inf)))
#'
#' # DSL string with wildcards and quantifiers
#' r3 <- getchainseq(spmeinvivo,
#'                   "162.0528; -18.0106{0,3}; 14.0157+")
#'
#' # Anchor to a known compound
#' r4 <- getchainseq(spmeinvivo, c(162.0528, -18.0106), mass = 286.3101)
#' }
#' @seealso \code{\link{getchain}}, \code{\link{gethomolog}},
#'   \code{\link{parse_pmd_pattern}}
#' @export
getchainseq <- function(list,
                        pattern,
                        mass = NULL,
                        digits = 4,
                        rtcutoff = 10,
                        corcutoff = 0.6,
                        ppm = 25,
                        rtdir = "any",
                        max_paths = 1e5,
                        allow_cycles = FALSE) {

        DT <- data.table::data.table

        # ---- Normalize pattern ----
        if (is.character(pattern) && length(pattern) == 1L) {
                pattern <- parse_pmd_pattern(pattern)
        }
        norm_pattern <- function(p) {
                if (is.numeric(p)) {
                        return(lapply(p, function(v)
                                list(pmd = v, min = 1L, max = 1L)))
                }
                if (is.list(p)) {
                        return(lapply(p, function(s) {
                                pv <- if (is.null(s$pmd)) NA_real_ else s$pmd
                                mn <- as.integer(if (is.null(s$min)) 1L else s$min)
                                mx <- if (is.null(s$max)) mn
                                else if (is.infinite(s$max)) Inf
                                else as.integer(s$max)
                                list(pmd = pv, min = mn, max = mx)
                        }))
                }
                stop("pattern must be a numeric vector, list of step specs, or DSL string.")
        }
        pat <- norm_pattern(pattern)

        # ---- Features & RT groups ----
        mz <- list$mz
        data_mat <- list$data
        n <- length(mz)
        has_rt <- !is.null(list$rt)
        rt <- if (has_rt) list$rt else rep(NA_real_, n)

        rtg <- if (has_rt) {
                .rt_clusters(rt, rtcutoff)
        } else {
                seq_len(n)
        }

        cormat <- stats::cor(t(data_mat))

        # ---- Build global directed edge table ----
        if (n > 10000L)
                warning("Large feature count (n=", n,
                        "); consider pre-filtering with globalstd first.")

        diffmat <- outer(mz, mz, function(a, b) b - a)
        rtgmat  <- outer(rtg, rtg, `!=`)
        diag(rtgmat) <- FALSE
        idx <- which(rtgmat, arr.ind = TRUE)

        edges_global <- DT(
                from = idx[, 1],
                to   = idx[, 2],
                pmd  = diffmat[idx],
                rtd  = if (has_rt) rt[idx[, 2]] - rt[idx[, 1]] else NA_real_,
                cor  = cormat[idx]
        )
        if (!is.null(corcutoff)) {
                keep <- abs(edges_global$cor) >= corcutoff
                edges_global <- edges_global[keep, ]
        }
        edges_global[, pmd_r := round(pmd, digits)]
        data.table::setkey(edges_global, pmd_r)
        rm(diffmat, rtgmat, idx); invisible(gc())

        get_step_edges <- function(spec, edge_idx) {
                e <- if (is.na(spec$pmd)) data.table::copy(edges_global)
                else edges_global[pmd_r == round(spec$pmd, digits)]
                if (nrow(e) == 0L) return(e)
                dir_val <- if (is.numeric(rtdir)) rtdir[edge_idx]
                else if (identical(rtdir, "increasing")) 1
                else if (identical(rtdir, "decreasing")) -1
                else 0
                if (has_rt && !is.na(dir_val) && dir_val != 0)
                        e <- e[sign(rtd) == dir_val]
                e
        }

        # ---- Seed filter ----
        seed_nodes <- seq_len(n)
        if (!is.null(mass)) {
                if (is.character(mass))
                        mass <- unlist(Map(enviGCMS::getmass, mass))
                keep <- rep(FALSE, n)
                for (m in mass) {
                        tol <- m * ppm / 1e6
                        keep <- keep | (mz >= m - tol & mz <= m + tol)
                }
                seed_nodes <- which(keep)
                if (length(seed_nodes) == 0L) {
                        message("No features match seed mass within ppm.")
                        list$sdacseq <- NULL
                        return(list)
                }
        }

        # ---- Path state ----
        # n1..nK : node indices (NA-padded when variable-length)
        # .end   : current end-node index (join key)
        # .len   : current number of nodes
        paths <- DT(n1 = seed_nodes, .end = seed_nodes, .len = 1L)

        extend_once <- function(paths, edges_step, new_col_name) {
                if (nrow(paths) == 0L || nrow(edges_step) == 0L)
                        return(paths[0])

                e <- edges_step[, .(.end = from, new_to = to)]
                joined <- merge(paths, e, by = ".end", allow.cartesian = TRUE)
                if (nrow(joined) == 0L) return(joined)

                if (!allow_cycles) {
                        ncols <- grep("^n[0-9]+$", colnames(joined), value = TRUE)
                        bad <- rep(FALSE, nrow(joined))
                        for (cc in ncols) {
                                vals <- joined[[cc]]
                                bad <- bad | (!is.na(vals) & vals == joined$new_to)
                        }
                        joined <- joined[!bad]
                }
                if (nrow(joined) == 0L) return(joined)

                data.table::setnames(joined, "new_to", new_col_name)
                joined[, .end := get(new_col_name)]
                joined[, .len := .len + 1L]

                if (nrow(joined) > max_paths) {
                        warning(sprintf("max_paths (%d) reached; truncating.", max_paths))
                        joined <- joined[seq_len(max_paths)]
                }
                joined
        }

        col_idx  <- 1L   # index of the last n-column added
        edge_idx <- 0L   # index of the current edge in the unrolled pattern

        for (k in seq_along(pat)) {
                spec <- pat[[k]]

                # Required repeats
                if (spec$min > 0L) {
                        for (r in seq_len(spec$min)) {
                                edge_idx <- edge_idx + 1L
                                es <- get_step_edges(spec, edge_idx)
                                if (nrow(es) == 0L) {
                                        message(sprintf(
                                                "Step %d rep %d: no edges for pmd=%s.",
                                                k, r, as.character(spec$pmd)))
                                        list$sdacseq <- NULL
                                        return(list)
                                }
                                col_idx <- col_idx + 1L
                                paths <- extend_once(paths, es,
                                                     paste0("n", col_idx))
                                if (nrow(paths) == 0L) {
                                        message(sprintf(
                                                "No paths survived step %d rep %d.", k, r))
                                        list$sdacseq <- NULL
                                        return(list)
                                }
                        }
                }

                # Optional repeats: accumulate snapshots of every valid length
                extra_cap <- if (is.infinite(spec$max)) Inf
                else spec$max - spec$min

                if (extra_cap > 0) {
                        accumulated <- list(data.table::copy(paths))
                        current <- data.table::copy(paths)
                        r <- 0L
                        while (r < extra_cap) {
                                r <- r + 1L
                                edge_idx <- edge_idx + 1L
                                es <- get_step_edges(spec, edge_idx)
                                if (nrow(es) == 0L) break

                                col_idx <- col_idx + 1L
                                ext <- extend_once(current, es,
                                                   paste0("n", col_idx))
                                if (nrow(ext) == 0L) {
                                        col_idx <- col_idx - 1L
                                        break
                                }
                                accumulated[[length(accumulated) + 1L]] <- ext
                                current <- ext
                                if (nrow(current) > max_paths) {
                                        warning("max_paths reached in optional extension.")
                                        break
                                }
                        }
                        paths <- data.table::rbindlist(accumulated, fill = TRUE,
                                                       use.names = TRUE)
                }
        }

        if (nrow(paths) == 0L) {
                list$sdacseq <- NULL
                return(list)
        }

        # ---- Decorate output ----
        node_cols <- grep("^n[0-9]+$", colnames(paths), value = TRUE)
        node_cols <- node_cols[order(as.integer(sub("^n", "", node_cols)))]

        for (nc in node_cols) {
                k <- sub("^n", "", nc)
                paths[, paste0("mz_", k) := mz[get(nc)]]
                if (has_rt)
                        paths[, paste0("rt_", k) := rt[get(nc)]]
        }

        if (length(node_cols) >= 2L) {
                for (i in seq_len(length(node_cols) - 1L)) {
                        a <- node_cols[i]
                        b <- node_cols[i + 1L]
                        paths[, paste0("pmd_", i) := mz[get(b)] - mz[get(a)]]
                }
        }

        paths[, cor_mean := {
                idxs <- unlist(.SD)
                idxs <- idxs[!is.na(idxs)]
                if (length(idxs) < 2L) NA_real_
                else mean(abs(cormat[cbind(idxs[-length(idxs)], idxs[-1L])]))
        }, by = seq_len(nrow(paths)), .SDcols = node_cols]

        paths[, path_len := .len]
        paths[, c(".end", ".len") := NULL]

        paths <- unique(paths)
        data.table::setcolorder(paths, c(node_cols, "path_len",
                                         grep("^mz_", colnames(paths), value = TRUE),
                                         grep("^rt_", colnames(paths), value = TRUE),
                                         grep("^pmd_", colnames(paths), value = TRUE),
                                         "cor_mean"))
        list$sdacseq <- paths
        list$pattern <- pat
        return(list)
}


#' Find homologous series by a repeating PMD unit
#'
#' Convenience wrapper around \code{\link{getchainseq}} for the common case of
#' detecting homologous series: chains of features linked by repeated
#' applications of a single PMD unit. Useful for finding alkyl chain series
#' (\code{+CH2 = 14.0157}), PEG series (\code{+C2H4O = 44.0262}), polymeric
#' artifacts, and similar patterns.
#'
#' @param list a pmd-style list with \code{mz}, \code{rt} (optional), \code{data}
#' @param unit the repeating PMD, default \code{14.0157} (CH2)
#' @param min_len minimum number of nodes in the series, default 3
#' @param max_len maximum number of nodes, default \code{Inf}
#' @param ... passed through to \code{\link{getchainseq}} (e.g. \code{corcutoff},
#'   \code{rtdir}, \code{mass}, \code{ppm})
#' @return as \code{\link{getchainseq}}
#' @examples
#' \dontrun{
#' data(spmeinvivo)
#' gethomolog(spmeinvivo, unit = 14.0157, min_len = 4)
#' gethomolog(spmeinvivo, unit = 44.0262, min_len = 3)
#' }
#' @seealso \code{\link{getchainseq}}
#' @export
gethomolog <- function(list,
                       unit = 14.0157,
                       min_len = 3L,
                       max_len = Inf,
                       ...) {
        stopifnot(min_len >= 2L, max_len >= min_len)
        mn <- as.integer(min_len) - 1L
        mx <- if (is.infinite(max_len)) Inf else as.integer(max_len) - 1L
        pat <- list(list(pmd = unit, min = mn, max = mx))
        getchainseq(list, pattern = pat, ...)
}
