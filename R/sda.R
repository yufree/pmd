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
                        data <- list$data[list$pairedindex,]
                        rtg <- list$rtcluster[list$pairedindex]
                } else {
                        mz <- list$mz[list$stdmassindex]
                        rt <- list$rt[list$stdmassindex]
                        data <- list$data[list$stdmassindex,]
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
                                        ms1 = mz[which(lower.tri(dis), arr.ind = T)[,
                                                                                    1]],
                                        ms2 = mz[which(lower.tri(dis), arr.ind = T)[,
                                                                                    2]],
                                        diff = as.numeric(dis),
                                        rt1 = rt[which(lower.tri(disrt),
                                                       arr.ind = T)[, 1]],
                                        rt2 = rt[which(lower.tri(disrt),
                                                       arr.ind = T)[, 2]],
                                        diffrt = as.numeric(disrt),
                                        rtg1 = rtg[which(lower.tri(disrtg),
                                                         arr.ind = T)[, 1]],
                                        rtg2 = rtg[which(lower.tri(disrtg),
                                                         arr.ind = T)[, 2]],
                                        rtgdiff = as.numeric(disrtg),
                                        cor = cor[lower.tri(cor)]
                                )
                }else{
                        df <- data.frame(
                                ms1 = mz[which(lower.tri(dis), arr.ind = T)[,
                                                                            1]],
                                ms2 = mz[which(lower.tri(dis), arr.ind = T)[,
                                                                            2]],
                                diff = as.numeric(dis),
                                rt1 = rt[which(lower.tri(disrt),
                                               arr.ind = T)[, 1]],
                                rt2 = rt[which(lower.tri(disrt),
                                               arr.ind = T)[, 2]],
                                diffrt = as.numeric(disrt),
                                rtg1 = rtg[which(lower.tri(disrtg),
                                                 arr.ind = T)[, 1]],
                                rtg2 = rtg[which(lower.tri(disrtg),
                                                 arr.ind = T)[, 2]],
                                rtgdiff = as.numeric(disrtg)
                        )
                }
                df <- df[df$rtgdiff > 0,]
                df$diff2 <- round(df$diff, digits)
                # use unique isomers
                index <-
                        !duplicated(paste0(
                                round(df$ms1, accuracy),
                                round(df$ms2, accuracy)
                        ))
                diff <- df$diff2[index]
                freq <-
                        sort(table(diff),decreasing = T)
                if (is.null(freqcutoff)) {
                        dis <- c()
                        for (i in c(1:ifelse(length(freq)>100,100,length(freq)))){
                                pmd <- as.numeric(names(freq))[1:i]
                                dfx <- df[df$diff2 %in% pmd,c(1,2)]
                                net <- igraph::graph_from_data_frame(dfx,directed = F)
                                dis[i] <- igraph::mean_distance(net)
                        }
                        n <- which.max(dis)
                        freqt <- freq[n-1]
                        if(n==100){
                                warning("Average distance is still increasing, you need to check manually for frequency cutoff.")
                        }
                        if(sum(df$diff2 == 0)>freqt & 0 %in% as.numeric(names(freq))){
                                list$sda <- df[df$diff2 %in% c(0,as.numeric(names(freq[freq>freqt]))),]
                        }else{
                                list$sda <- df[df$diff2 %in% as.numeric(names(freq[freq>freqt])),]
                        }
                        message(paste("PMD frequency cutoff is", freqt, 'by PMD network analysis with largest network average distance',round(max(dis),2),'.'))
                        # i <- dis <- t <- 1
                        # while(n>=t){
                        #         pmd <- as.numeric(names(freq[freq>i]))
                        #         dfx <- df[df$diff2 %in% pmd,c(1,2)]
                        #         net <- igraph::graph_from_data_frame(dfx,directed = F)
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
                }else{
                        if (sum(df$diff2 == 0) > freqcutoff &
                            0 %in% as.numeric(names(freq))) {
                                list$sda <- df[(df$diff2 %in% c(0, as.numeric(names(
                                        freq[freq >=
                                                     freqcutoff]
                                )))),]
                        } else{
                                list$sda <- df[(df$diff2 %in% c(as.numeric(names(
                                        freq[freq >=
                                                     freqcutoff]
                                )))),]
                        }
                }

                if (!is.null(corcutoff)&!is.null(data)) {
                        list$sda <- list$sda[abs(list$sda$cor) >= corcutoff, ]
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
#' @param mz numeric vector for independant mass or mass to charge ratio. Mass to charge ratio from GlobalStd algorithm is suggested. Isomers would be excluded automately
#' @param freqcutoff pmd freqency cutoff for structures or reactions, default 10
#' @param digits mass or mass to charge ratio accuracy for pmd, default 3
#' @param top top n pmd freqency cutoff when the freqcutoff is too small for large data set
#' @param formula vector for formula when you don't have mass or mass to charge ratio data
#' @return logical matrix with row as the same order of mz or formula and column as high freqency pmd group
#' @examples
#' data(spmeinvivo)
#' pmd <- getpaired(spmeinvivo)
#' std <- getstd(pmd)
#' sda <- getrda(spmeinvivo$mz[std$stdmassindex])
#' @seealso \code{\link{getsda}}
#' @export
getrda <-
        function(mz,
                 freqcutoff = 10,
                 digits = 3,
                 top = 20,
                 formula = NULL) {

                if (is.null(formula)) {
                        mz <- unique(mz)
                        dis <- stats::dist(mz, method = "manhattan")
                } else{
                        mz <- unlist(Map(enviGCMS::getmass, formula))
                        mz <- unique(mz)
                        dis <- stats::dist(mz, method = "manhattan")

                }

                df <- cbind.data.frame(
                        ms1 = mz[which(lower.tri(dis), arr.ind = T)[, 1]],
                        ms2 = mz[which(lower.tri(dis), arr.ind = T)[, 2]],
                        diff = as.numeric(dis),
                        diff2 = round(as.numeric(dis), digits = digits)
                )
                freq <-
                        sort(table(df$diff2),decreasing = T)
                message(paste(length(freq), 'pmd found.'))
                if (!is.null(top)) {
                        freq <- utils::head(freq, top)
                }
                sda <-
                        df[(df$diff2 %in% c(as.numeric(names(freq[freq >= freqcutoff])))),]
                # mz <- unique(c(sda$ms1,sda$ms2))
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
                return(as.data.frame(result))
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
                   accuracy = 4){

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
                        ms1 = mz[which(lower.tri(dis), arr.ind = T)[,
                                                                    1]],
                        ms2 = mz[which(lower.tri(dis), arr.ind = T)[,
                                                                    2]],
                        diff = as.numeric(dis),
                        rt1 = rt[which(lower.tri(disrt),
                                       arr.ind = T)[, 1]],
                        rt2 = rt[which(lower.tri(disrt),
                                       arr.ind = T)[, 2]],
                        diffrt = as.numeric(disrt),
                        rtg1 = rtg[which(lower.tri(disrtg),
                                         arr.ind = T)[, 1]],
                        rtg2 = rtg[which(lower.tri(disrtg),
                                         arr.ind = T)[, 2]],
                        rtgdiff = as.numeric(disrtg),
                        cor = cor[lower.tri(cor)]
                )
        list$cda <- df[abs(df$cor) >= corcutoff, ]
        return(list)
}

#' Get pmd for specific reaction
#' @param list a list with mzrt profile
#' @param pmd a specific paired mass distance
#' @param rtcutoff cutoff of the distances in retention time hierarchical clustering analysis, default 10
#' @param corcutoff cutoff of the correlation coefficient, default NULL
#' @param digits mass or mass to charge ratio accuracy for pmd, default 2
#' @param accuracy measured mass or mass to charge ratio in digits, default 4
#' @return list with paired peaks for specific pmd.
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

                cor <- stats::cor(t(data))

                df <- data.frame(
                        ms1 = mz[which(lower.tri(dis), arr.ind = T)[,
                                                                    1]],
                        ms2 = mz[which(lower.tri(dis), arr.ind = T)[,
                                                                    2]],
                        diff = as.numeric(dis),
                        rt1 = rt[which(lower.tri(disrt),
                                       arr.ind = T)[, 1]],
                        rt2 = rt[which(lower.tri(disrt),
                                       arr.ind = T)[, 2]],
                        diffrt = as.numeric(disrt),
                        rtg1 = rtg[which(lower.tri(disrtg),
                                         arr.ind = T)[, 1]],
                        rtg2 = rtg[which(lower.tri(disrtg),
                                         arr.ind = T)[, 2]],
                        rtgdiff = as.numeric(disrtg),
                        cor = cor[lower.tri(cor)]
                )

                if (!is.null(corcutoff)) {
                        df <- df[abs(df$cor) >= corcutoff, ]
                }

                df$diff2 <- round(df$diff, digits)

                df <- df[df$rtgdiff > 0 & df$diff2 == pmd,]
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
                return(list)
        }
#' Get reaction chain for specific mass to charge ratio
#' @param list a list with mzrt profile
#' @param diff paired mass distance(s) of interests
#' @param mass a specific mass for known compound or a vector of masses. You could also input formula for certain compounds
#' @param digits mass or mass to charge ratio accuracy for pmd, default 2
#' @param accuracy measured mass or mass to charge ratio in digits, default 4
#' @param rtcutoff cutoff of the distances in retention time hierarchical clustering analysis, default 10
#' @param corcutoff cutoff of the correlation coefficient, default NULL
#' @param ppm all the peaks within this mass accuracy as seed mass or formula
#' @return a list with mzrt profile and reaction chain dataframe
#' @examples
#' data(spmeinvivo)
#' # check metabolites of C18H39NO
#' pmd <- getchain(spmeinvivo,diff = c(2.02,14.02,15.99),mass = 286.3101)
#' @export
getchain <- function(list, diff, mass, digits = 2, accuracy = 4, rtcutoff= 10, corcutoff=0.6,ppm=25) {
        if (is.character(mass)) {
                mass <- unlist(Map(enviGCMS::getmass, mass))
        }
        massup <- mass+mass*ppm/10e6
        massdown <- mass-mass*ppm/10e6
        updown <- sapply(Map(function(x)
                x < massup & x > massdown, list$mz), function(x)
                        sum(x&T)>0)
        mass <- list$mz[updown]
        mass <- unique(round(mass,accuracy))

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
        cor <- stats::cor(t(data))

        df <- data.frame(
                ms1 = mz[which(lower.tri(dis), arr.ind = T)[,
                                                            1]],
                ms2 = mz[which(lower.tri(dis), arr.ind = T)[,
                                                            2]],
                diff = as.numeric(dis),
                rt1 = rt[which(lower.tri(disrt),
                               arr.ind = T)[, 1]],
                rt2 = rt[which(lower.tri(disrt),
                               arr.ind = T)[, 2]],
                diffrt = as.numeric(disrt),
                rtg1 = rtg[which(lower.tri(disrtg),
                                 arr.ind = T)[, 1]],
                rtg2 = rtg[which(lower.tri(disrtg),
                                 arr.ind = T)[, 2]],
                rtgdiff = as.numeric(disrtg),
                cor = cor[lower.tri(cor)]
        )

        if (!is.null(corcutoff)) {
                df <- df[abs(df$cor) >= corcutoff, ]
        }

        df$diff2 <- round(df$diff, digits)

        df <- df[df$rtgdiff > 0,]
        ms1 <- ifelse(df$ms1 > df$ms2, df$ms1, df$ms2)
        ms2 <- ifelse(df$ms1 > df$ms2, df$ms2, df$ms1)
        rtg1 <- ifelse(df$ms1 > df$ms2, df$rtg1, df$rtg2)
        rtg2 <- ifelse(df$ms1 > df$ms2, df$rtg2, df$rtg1)

        sda <- df[df$diff2 %in% diff,]
        seed <- NULL
        ms1 <- round(sda$ms1, digits = accuracy)
        ms2 <- round(sda$ms2, digits = accuracy)
        if (length(mass) == 1) {
                mass <- round(mass, accuracy)
                sdat <-
                        unique(c(mass, ms2[ms1 %in% mass], ms1[ms2 %in% mass]))
                while (!identical(sdat, seed)) {
                        seed <- sdat
                        sdat <-
                                unique(c(sdat, ms2[ms1 %in% sdat], ms1[ms2 %in% sdat]))
                }
                list$sdac <- sda[ms1 %in% sdat | ms2 %in% sdat ,]
                return(list)
        } else if (length(mass) == 0) {
                warning(
                        'No mass input and all mass in the list will be used for reaction chain construction!'
                )
                sdac <- NULL
                mass <- round(list$mz, accuracy)
                for (i in 1:length(mass)) {
                        sdat <-
                                unique(c(mass[i], ms2[ms1 %in% mass[i]], ms1[ms2 %in% mass[i]]))
                        if (length(sdat) != 1) {
                                while (!identical(sdat, seed)) {
                                        seed <- sdat
                                        sdat <-
                                                unique(c(sdat, ms2[ms1 %in% sdat], ms1[ms2 %in% sdat]))
                                }
                                sdact <-
                                        sda[ms1 %in% sdat | ms2 %in% sdat ,]
                                sdact$mass <- mass[i]
                                sdac <- rbind.data.frame(sdac, sdact)
                        }
                }
                list$sdac <- sdac[!duplicated(sdac), ]
                return(list)
        } else{
                sdac <- NULL
                mass <- round(mass, accuracy)
                for (i in 1:length(mass)) {
                        sdat <-
                                unique(c(mass[i], ms2[ms1 %in% mass[i]], ms1[ms2 %in% mass[i]]))
                        if (length(sdat) != 1) {
                                while (!identical(sdat, seed)) {
                                        seed <- sdat
                                        sdat <-
                                                unique(c(sdat, ms2[ms1 %in% sdat], ms1[ms2 %in% sdat]))
                                }
                                sdact <- sda[ms1 %in% sdat | ms2 %in% sdat,]
                                sdact$mass <- mass[i]
                                sdac <- rbind.data.frame(sdac, sdact)
                        }
                }
                list$sdac <- sdac[!duplicated(sdac), ]
                return(list)
        }

}

#' Get quantitative paired peaks list for specific reaction/pmd
#' @param list a list with mzrt profile and data
#' @param pmd a specific paired mass distances
#' @param rtcutoff cutoff of the distances in retention time hierarchical clustering analysis, default 10
#' @param digits mass or mass to charge ratio accuracy for pmd, default 2
#' @param accuracy measured mass or mass to charge ratio in digits, default 4
#' @param ratiocv ratio cv cutoff for quantitative paired peaks, default 30
#' @param outlier logical, if true, outliar of ratio will be removed, default False.
#' @param ... other parameters for getpmd
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
                 ratiocv = 30,
                 outlier = F,
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
                getr <- function(v) {
                        ratio <- NULL
                        ratio1 <-
                                data[list$mz %in% v[1] & list$rt %in% v[4], ]
                        ratio2 <-
                                data[list$mz %in% v[2] & list$rt %in% v[5], ]
                        ratio <- as.numeric(ratio1) / as.numeric(ratio2)
                        if(outlier){
                                outlier_values <- grDevices::boxplot.stats(ratio)$out
                                ratio <- ratio[!ratio%in%outlier_values]
                        }
                        rsd <-
                                stats::sd(ratio, na.rm = T) / mean(ratio, na.rm = T) * 100
                        return(rsd)
                }
                if(sum(p$pmdindex)>0){
                        list <- enviGCMS::getfilter(p, p$pmdindex)
                        data <- list$data
                        pmd <- list$pmd
                        list$pmd$r <- apply(pmd, 1, getr)
                        list$pmd <- list$pmd[list$pmd$r < ratiocv , ]
                        list$pmd <- list$pmd[stats::complete.cases(list$pmd),]
                        if(nrow(list$pmd)>0){
                                idx <- paste(list$mz, list$rt)
                                idx2 <-
                                        unique(c(
                                                paste(list$pmd$ms1, list$pmd$rt1),
                                                paste(list$pmd$ms2, list$pmd$rt2)
                                        ))
                                list <- enviGCMS::getfilter(list, idx %in% idx2)
                                return(list)
                        }else{
                                message('No quantitative peaks could be used.')
                        }

                }else{
                        message('No pmd peaks could be found.')
                }
        }

