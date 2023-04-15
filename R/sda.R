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
                if (is.null(list$stdmass) & is.null(list$paired)) {
                        mz <- list$mz
                        rt <- list$rt
                        data <- list$data
                        dis <- stats::dist(rt, method = "manhattan")
                        fit <- stats::hclust(dis)
                        rtg <- stats::cutree(fit, h = rtcutoff)
                } else if (is.null(list$stdmass)) {
                        mz <- list$mz[list$pairedindex]
                        rt <- list$rt[list$pairedindex]
                        data <- list$data[list$pairedindex, ]
                        rtg <- list$rtcluster[list$pairedindex]
                } else {
                        mz <- list$mz[list$stdmassindex]
                        rt <- list$rt[list$stdmassindex]
                        data <- list$data[list$stdmassindex, ]
                        rtg <- list$rtcluster[list$stdmassindex]
                }
                # PMD analysis
                dis <- stats::dist(mz, method = "manhattan")
                disrt <- stats::dist(rt, method = "manhattan")
                disrtg <- stats::dist(rtg, method = "manhattan")

                if (!is.null(data)) {
                        cor <- stats::cor(t(data))
                        df <-
                                data.frame(
                                        ms1 = mz[which(lower.tri(dis), arr.ind = TRUE)[,
                                                                                       1]],
                                        ms2 = mz[which(lower.tri(dis), arr.ind = TRUE)[,
                                                                                       2]],
                                        diff = as.numeric(dis),
                                        rt1 = rt[which(lower.tri(disrt),
                                                       arr.ind = TRUE)[, 1]],
                                        rt2 = rt[which(lower.tri(disrt),
                                                       arr.ind = TRUE)[, 2]],
                                        diffrt = as.numeric(disrt),
                                        rtg1 = rtg[which(lower.tri(disrtg),
                                                         arr.ind = TRUE)[, 1]],
                                        rtg2 = rtg[which(lower.tri(disrtg),
                                                         arr.ind = TRUE)[, 2]],
                                        rtgdiff = as.numeric(disrtg),
                                        md = as.numeric(dis)%%1,
                                        cor = cor[lower.tri(cor)]
                                )
                } else{
                        df <- data.frame(
                                ms1 = mz[which(lower.tri(dis), arr.ind = TRUE)[,
                                                                               1]],
                                ms2 = mz[which(lower.tri(dis), arr.ind = TRUE)[,
                                                                               2]],
                                diff = as.numeric(dis),
                                rt1 = rt[which(lower.tri(disrt),
                                               arr.ind = TRUE)[, 1]],
                                rt2 = rt[which(lower.tri(disrt),
                                               arr.ind = TRUE)[, 2]],
                                diffrt = as.numeric(disrt),
                                rtg1 = rtg[which(lower.tri(disrtg),
                                                 arr.ind = TRUE)[, 1]],
                                rtg2 = rtg[which(lower.tri(disrtg),
                                                 arr.ind = TRUE)[, 2]],
                                rtgdiff = as.numeric(disrtg),
                                md = as.numeric(dis)%%1
                        )
                }
                df <- df[df$rtgdiff > 0, ]
                df$diff2 <- round(df$diff, digits)
                # use unique isomers
                index <-
                        !duplicated(paste0(round(df$ms1, accuracy),
                                           round(df$ms2, accuracy)))
                diff <- df$diff2[index]
                freq <-
                        sort(table(diff), decreasing = TRUE)
                if (is.null(freqcutoff)) {
                        dis <- c()
                        for (i in c(1:ifelse(length(freq) > 100, 100, length(freq)))) {
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
                if (is.null(formula)) {
                        mz <- unique(mz)
                        dis <- stats::dist(mz, method = "manhattan")
                } else{
                        mz <- unlist(Map(enviGCMS::getmass, formula))
                        mz <- unique(mz)
                        dis <- stats::dist(mz, method = "manhattan")

                }

                df <- cbind.data.frame(
                        ms1 = mz[which(lower.tri(dis), arr.ind = TRUE)[, 1]],
                        ms2 = mz[which(lower.tri(dis), arr.ind = TRUE)[, 2]],
                        diff = as.numeric(dis),
                        diff2 = round(as.numeric(dis), digits = digits),
                        md = as.numeric(dis)%%1
                )
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
#' cluster <- getcorcluster(spmeinvivo)
#' cbp <- enviGCMS::getfilter(cluster,rowindex = cluster$stdmassindex2)
#' cda <- getcda(cbp)
#' @seealso \code{\link{getsda}},\code{\link{getrda}}
#' @export
getcda <- function(list,
                   corcutoff = 0.9,
                   rtcutoff = 10,
                   accuracy = 4) {
        mz <- list$mz
        rt <- list$rt
        data <- list$data
        dis <- stats::dist(rt, method = "manhattan")
        fit <- stats::hclust(dis)
        rtg <- stats::cutree(fit, h = rtcutoff)
        cor <- stats::cor(t(data))
        disrt <- stats::dist(rt, method = "manhattan")
        disrtg <- stats::dist(rtg, method = "manhattan")
        df <-
                data.frame(
                        ms1 = mz[which(lower.tri(dis), arr.ind = TRUE)[,
                                                                       1]],
                        ms2 = mz[which(lower.tri(dis), arr.ind = TRUE)[,
                                                                       2]],
                        diff = as.numeric(dis),
                        rt1 = rt[which(lower.tri(disrt),
                                       arr.ind = TRUE)[, 1]],
                        rt2 = rt[which(lower.tri(disrt),
                                       arr.ind = TRUE)[, 2]],
                        diffrt = as.numeric(disrt),
                        rtg1 = rtg[which(lower.tri(disrtg),
                                         arr.ind = TRUE)[, 1]],
                        rtg2 = rtg[which(lower.tri(disrtg),
                                         arr.ind = TRUE)[, 2]],
                        rtgdiff = as.numeric(disrtg),
                        md = as.numeric(dis)%%1,
                        cor = cor[lower.tri(cor)]
                )
        list$cda <- df[abs(df$cor) >= corcutoff,]
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
                # PMD analysis
                # remove isomers
                if(!is.null(list$rt)){
                        rt <- list$rt
                        dis <- stats::dist(rt, method = "manhattan")
                        fit <- stats::hclust(dis)
                        rtg <- stats::cutree(fit, h = rtcutoff)
                        disrt <- stats::dist(rt, method = "manhattan")
                        disrtg <- stats::dist(rtg, method = "manhattan")
                        dis <- stats::dist(mz, method = "manhattan")
                        cor <- stats::cor(t(data))
                        df <- data.frame(
                                ms1 = mz[which(lower.tri(dis), arr.ind = TRUE)[,
                                                                               1]],
                                ms2 = mz[which(lower.tri(dis), arr.ind = TRUE)[,
                                                                               2]],
                                diff = as.numeric(dis),
                                rt1 = rt[which(lower.tri(disrt),
                                               arr.ind = TRUE)[, 1]],
                                rt2 = rt[which(lower.tri(disrt),
                                               arr.ind = TRUE)[, 2]],
                                diffrt = as.numeric(disrt),
                                rtg1 = rtg[which(lower.tri(disrtg),
                                                 arr.ind = TRUE)[, 1]],
                                rtg2 = rtg[which(lower.tri(disrtg),
                                                 arr.ind = TRUE)[, 2]],
                                rtgdiff = as.numeric(disrtg),
                                cor = cor[lower.tri(cor)]
                        )
                        if (!is.null(corcutoff)) {
                                df <- df[abs(df$cor) >= corcutoff,]
                        }

                        df$diff2 <- round(df$diff, digits)

                        df <- df[df$rtgdiff > 0 & df$diff2 %in% pmd, ]
                        ms1 <- ifelse(df$ms1 > df$ms2, df$ms1, df$ms2)
                        ms2 <- ifelse(df$ms1 > df$ms2, df$ms2, df$ms1)
                        rtg1 <- ifelse(df$ms1 > df$ms2, df$rtg1, df$rtg2)
                        rtg2 <- ifelse(df$ms1 > df$ms2, df$rtg2, df$rtg1)
                        list$pmd <- df

                        index <-
                                c(paste(round(ms1, accuracy), rtg1), paste(round(ms2, accuracy), rtg2))
                        index <- unique(index)
                        indexh <- paste(round(ms1, accuracy), rtg1)
                        indexh <- unique(indexh)
                        indexl <- paste(round(ms2, accuracy), rtg2)
                        indexl <- unique(indexl)

                        index0 <- paste(round(list$mz, accuracy), rtg)
                        list$pmdindex <- index0 %in% index
                        list$pmdindexh <- index0 %in% indexh
                        list$pmdindexl <- index0 %in% indexl
                }else{
                        dis <- stats::dist(mz, method = "manhattan")
                        cor <- stats::cor(t(data))
                        df <- data.frame(
                                ms1 = mz[which(lower.tri(dis), arr.ind = TRUE)[,
                                                                               1]],
                                ms2 = mz[which(lower.tri(dis), arr.ind = TRUE)[,
                                                                               2]],
                                diff = as.numeric(dis),
                                cor = cor[lower.tri(cor)]
                        )
                        if (!is.null(corcutoff)) {
                                df <- df[abs(df$cor) >= corcutoff,]
                        }

                        df$diff2 <- round(df$diff, digits)

                        df <- df[df$diff2 %in% pmd, ]
                        ms1 <- ifelse(df$ms1 > df$ms2, df$ms1, df$ms2)
                        ms2 <- ifelse(df$ms1 > df$ms2, df$ms2, df$ms1)
                        list$pmd <- df
                        index <-
                                c(round(ms1, accuracy), round(ms2, accuracy))
                        index <- unique(index)
                        indexh <- unique(round(ms1, accuracy))
                        indexl <- unique(round(ms2, accuracy))

                        index0 <- round(list$mz, accuracy)
                        list$pmdindex <- index0 %in% index
                        list$pmdindexh <- index0 %in% indexh
                        list$pmdindexl <- index0 %in% indexl
                }
                return(list)
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
                massup <- mass + mass * ppm / 10e6
                massdown <- mass - mass * ppm / 10e6
                updown <- vapply(Map(function(x)
                        x < massup & x > massdown, list$mz), function(x)
                                sum(x & TRUE) > 0, TRUE)
                mass <- list$mz[updown]
                mass <- unique(round(mass, accuracy))

                mz <- list$mz
                rt <- list$rt
                data <- list$data
                dis <- stats::dist(rt, method = "manhattan")
                fit <- stats::hclust(dis)
                rtg <- stats::cutree(fit, h = rtcutoff)
                # PMD analysis
                # remove isomers
                dis <- stats::dist(mz, method = "manhattan")
                disrt <- stats::dist(rt, method = "manhattan")
                disrtg <- stats::dist(rtg, method = "manhattan")
                diffx <- as.numeric(dis)
                diff2 <- round(diffx, digits)
                diffrt <- as.numeric(disrt)
                rtgdiff <- as.numeric(disrtg)
                cor <- stats::cor(t(data))

                idx <- diff2 %in% diff
                idx2 <- rtgdiff > 0
                idx3 <- idx&idx2

                ms1 = mz[which(lower.tri(dis), arr.ind = TRUE)[, 1]][idx3]
                ms2 = mz[which(lower.tri(dis), arr.ind = TRUE)[, 2]][idx3]
                rt1 = rt[which(lower.tri(disrt), arr.ind = TRUE)[, 1]][idx3]
                rt2 = rt[which(lower.tri(disrt),arr.ind = TRUE)[, 2]][idx3]
                rtg1 = rtg[which(lower.tri(disrtg),arr.ind = TRUE)[, 1]][idx3]
                rtg2 = rtg[which(lower.tri(disrtg),arr.ind = TRUE)[, 2]][idx3]
                diffx = diffx[idx3]
                diff2 = diff2[idx3]
                diffrt = diffrt[idx3]
                rtgdiff = rtgdiff[idx3]
                cor = cor[lower.tri(cor)][idx3]

                df <- data.frame(
                        ms1 = ms1,
                        ms2 = ms2,
                        diff = diffx,
                        rt1 = rt1,
                        rt2 = rt2,
                        diffrt = diffrt,
                        rtg1 = rtg1,
                        rtg2 = rtg2,
                        rtgdiff = rtgdiff,
                        cor = cor,
                        diff2 = diff2
                )

                if (!is.null(corcutoff)) {
                        df <- df[abs(df$cor) >= corcutoff,]
                }
                ms1 <- ifelse(df$ms1 > df$ms2, df$ms1, df$ms2)
                ms2 <- ifelse(df$ms1 > df$ms2, df$ms2, df$ms1)
                rtg1 <- ifelse(df$ms1 > df$ms2, df$rtg1, df$rtg2)
                rtg2 <- ifelse(df$ms1 > df$ms2, df$rtg2, df$rtg1)

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
                # with retention time
                getr <- function(v) {
                        ratio <- NULL
                        ratio1 <-
                                data[list$mz %in% v[1] &
                                             list$rt %in% v[4],]
                        ratio2 <-
                                data[list$mz %in% v[2] &
                                             list$rt %in% v[5],]
                        ratio <-
                                as.numeric(ratio1) / as.numeric(ratio2)
                        if (outlier) {
                                outlier_values <- grDevices::boxplot.stats(ratio)$out
                                ratio <-
                                        ratio[!ratio %in% outlier_values]
                        }
                        rsd <-
                                stats::sd(ratio, na.rm = TRUE) / mean(ratio, na.rm = TRUE) * 100
                        rsdh <-
                                stats::sd(as.numeric(ratio1), na.rm = TRUE) / mean(as.numeric(ratio1), na.rm = TRUE) * 100
                        rsdl <-
                                stats::sd(as.numeric(ratio2), na.rm = TRUE) / mean(as.numeric(ratio2), na.rm = TRUE) * 100
                        return(list(
                                rsd = unlist(rsd),
                                rsdh = unlist(rsdh),
                                rsdl = unlist(rsdl)
                        ))
                }
                # without retention time
                getr2 <- function(v) {
                        ratio <- NULL
                        ratio1 <-
                                data[list$mz %in% v[1],]
                        ratio2 <-
                                data[list$mz %in% v[2],]
                        ratio <-
                                as.numeric(ratio1) / as.numeric(ratio2)
                        if (outlier) {
                                outlier_values <- grDevices::boxplot.stats(ratio)$out
                                ratio <-
                                        ratio[!ratio %in% outlier_values]
                        }
                        rsd <-
                                stats::sd(ratio, na.rm = TRUE) / mean(ratio, na.rm = TRUE) * 100
                        rsdh <-
                                stats::sd(as.numeric(ratio1), na.rm = TRUE) / mean(as.numeric(ratio1), na.rm = TRUE) * 100
                        rsdl <-
                                stats::sd(as.numeric(ratio2), na.rm = TRUE) / mean(as.numeric(ratio2), na.rm = TRUE) * 100
                        return(list(
                                rsd = unlist(rsd),
                                rsdh = unlist(rsdh),
                                rsdl = unlist(rsdl)
                        ))
                }
                if (sum(p$pmdindex) > 0) {
                        list <- enviGCMS::getfilter(p, p$pmdindex)
                        data <- list$data
                        pmd <- list$pmd

                        if(!is.null(list$rt)){
                                ratio <- apply(pmd, 1, getr)
                        }else{
                                ratio <- apply(pmd, 1, getr2)
                        }
                        list$pmd$r <-
                                sapply(ratio, function(x)
                                        x$rsd)
                        list$pmd$rh <-
                                sapply(ratio, function(x)
                                        x$rsdh)
                        list$pmd$rl <-
                                sapply(ratio, function(x)
                                        x$rsdl)
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
                                                                idx), ]
                                        pmdl <-
                                                list$data[match(paste(list$pmd$ms2, list$pmd$rt2),
                                                                idx), ]
                                        list$pmddata <- pmdh + pmdl
                                        return(list)
                                } else if (nrow(list$pmd) > 0&is.null(list$rt)){
                                        mzl <- list$pmd$rh<cvcutoff & list$pmd$rl<cvcutoff
                                        idx <- unique(c(list$pmd$ms1, list$pmd$ms2))
                                        list <-
                                                enviGCMS::getfilter(list, list$mz %in% idx)
                                        pmdh <-
                                                list$data[match(list$pmd$ms1,list$mz), ]
                                        pmdl <-
                                                list$data[match(list$pmd$ms2,list$mz), ]
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
                                                idx), ]/list$data[match(paste(
                                                        list$pmd$ms1[idy],
                                                        list$pmd$rt1[idy]
                                                ),
                                                idx), ]

                                        pmddata[!idy, ] <-
                                                list$data[match(paste(
                                                        list$pmd$ms1[!idy],
                                                        list$pmd$rt1[!idy]
                                                ),
                                                idx), ]/list$data[match(paste(
                                                        list$pmd$ms2[!idy],
                                                        list$pmd$rt2[!idy]
                                                ),
                                                idx), ]
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
                                                        list$mz), ]/list$data[match(
                                                                list$pmd$ms1[idy],
                                                                list$mz), ]
                                        pmddata[!idy, ] <-
                                                list$data[match(list$pmd$ms1[!idy],
                                                                list$mz), ]/list$data[match(list$pmd$ms2[!idy],
                                                                                            list$mz), ]
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
