#' Filter ions/peaks based on retention time hierarchical clustering, paired mass distances(PMD) and PMD frequency analysis.
#' @param list a list with mzrt profile
#' @param rtcutoff cutoff of the distances in retention time hierarchical clustering analysis, default 10
#' @param ng cutoff of global PMD's retention time group numbers, default NULL
#' @param digits mass or mass to charge ratio accuracy for pmd, default 2
#' @param accuracy measured mass or mass to charge ratio in digits, default 4
#' @return list with tentative isotope, multi-chargers, adducts, and neutral loss peaks' index, retention time clusters.
#' @examples
#' data(spmeinvivo)
#' pmd <- getpaired(spmeinvivo)
#' @seealso \code{\link{getstd}},\code{\link{getsda}},\code{\link{plotpaired}}
#' @export
getpaired <-
        function(list,
                 rtcutoff = 10,
                 ng = NULL,
                 digits = 2,
                 accuracy = 4) {
                # paired mass diff analysis
                dis <- stats::dist(list$rt, method = "manhattan")
                fit <- stats::hclust(dis)
                rtcluster <- stats::cutree(fit, h = rtcutoff)
                n <- length(unique(rtcluster))
                message(paste(n, "retention time cluster found."))
                # automate ng selection when ng is NULL
                ng <- ifelse(is.null(ng), round(n * 0.2), ng)
                if (!is.null(list$data)) {
                        groups <- cbind.data.frame(mz = list$mz,
                                                   rt = list$rt,
                                                   list$data)
                } else {
                        groups <- cbind.data.frame(mz = list$mz,
                                                   rt = list$rt)
                        message('No intensity data!')
                }

                split <- split.data.frame(groups, rtcluster)

                rtpmd <- function(bin, i) {
                        medianrtxi <- stats::median(bin$rt)
                        if (ncol(bin) > 2) {
                                if (nrow(bin) > 1) {
                                        # get mz diff
                                        dis <-
                                                stats::dist(bin$mz, method = "manhattan")
                                        # get intensity cor
                                        cor <-
                                                stats::cor(t(bin[,-c(1, 2)]))

                                        df <-
                                                data.frame(
                                                        ms1 = bin$mz[which(lower.tri(dis),
                                                                           arr.ind = T)[, 1]],
                                                        ms2 = bin$mz[which(lower.tri(dis),
                                                                           arr.ind = T)[, 2]],
                                                        diff = as.numeric(dis),
                                                        rt = medianrtxi,
                                                        rtg = i,
                                                        cor = cor[lower.tri(cor)]
                                                )
                                        # remove multi chargers
                                        multiindex <-
                                                (round(df$diff %% 1, 1) == 0.5)
                                        multiindex2 <-
                                                (round(df$diff, 1) == 0.5)
                                        mass <-
                                                unique(df[multiindex, 1], df[multiindex, 2])
                                        mass2 <-
                                                unique(df[multiindex2, 1], df[multiindex2, 2])
                                        multimass <-
                                                mass[round(mass %% 1, 1) == 0.5 |
                                                             round(mass, accuracy) %in% round(mass2, accuracy)]
                                        dfmulti <- df[multiindex,]
                                        if (nrow(dfmulti) > 0) {
                                                df <-
                                                        df[!(df[, 1] %in% multimass) &
                                                                   !(df[, 2] %in% multimass), ]
                                        }
                                        # remove isotope
                                        isoindex <-
                                                (round(df$diff, digits) != 0) &
                                                ((
                                                        df$diff %% 1 < 0.01 &
                                                                df$diff >= 1 &
                                                                df$diff < 2
                                                ) | (
                                                        df$diff %% 2 < 0.01 &
                                                                df$diff >= 2 &
                                                                df$diff < 3
                                                )
                                                )
                                        massstd <-
                                                apply(df[isoindex, ], 1, function(x)
                                                        min(x[1], x[2]))
                                        massstdmax <-
                                                apply(df[isoindex, ], 1, function(x)
                                                        max(x[1], x[2]))
                                        isomass <-
                                                unique(c(massstd[(massstd %in% massstdmax)], massstdmax))
                                        dfiso <- df[isoindex,]
                                        if (nrow(dfiso) > 0) {
                                                df <-
                                                        df[!(df[, 1] %in% isomass) &
                                                                   !(df[, 2] %in% isomass),]
                                        }
                                        dfdiff <- df
                                        return(
                                                list(
                                                        dfmulti = dfmulti,
                                                        dfiso = dfiso,
                                                        dfdiff = dfdiff,
                                                        solo = NULL
                                                )
                                        )
                                } else {
                                        solo <- cbind(bin[, c(1:2)],
                                                      rtg = i,
                                                      cor = 1)
                                        return(
                                                list(
                                                        dfmulti = NULL,
                                                        dfiso = NULL,
                                                        dfdiff = NULL,
                                                        solo = solo
                                                )
                                        )
                                }
                        } else{
                                if (nrow(bin) > 1) {
                                        # get mz diff
                                        dis <-
                                                stats::dist(bin$mz, method = "manhattan")
                                        df <-
                                                data.frame(
                                                        ms1 = bin$mz[which(lower.tri(dis),
                                                                           arr.ind = T)[, 1]],
                                                        ms2 = bin$mz[which(lower.tri(dis),
                                                                           arr.ind = T)[, 2]],
                                                        diff = as.numeric(dis),
                                                        rt = medianrtxi,
                                                        rtg = i
                                                )
                                        # remove multi chargers
                                        multiindex <-
                                                (round(df$diff %% 1, 1) == 0.5)
                                        multiindex2 <-
                                                (round(df$diff, 1) == 0.5)
                                        mass <-
                                                unique(df[multiindex, 1], df[multiindex, 2])
                                        mass2 <-
                                                unique(df[multiindex2, 1], df[multiindex2, 2])
                                        multimass <-
                                                mass[round(mass %% 1, 1) == 0.5 |
                                                             round(mass, accuracy) %in% round(mass2, accuracy)]
                                        dfmulti <- df[multiindex,]
                                        if (nrow(dfmulti) > 0) {
                                                df <-
                                                        df[!(df[, 1] %in% multimass) &
                                                                   !(df[, 2] %in% multimass), ]
                                        }
                                        # remove isotope
                                        isoindex <-
                                                (round(df$diff, digits) != 0) &
                                                ((
                                                        df$diff %% 1 < 0.01 &
                                                                df$diff >= 1 &
                                                                df$diff < 2
                                                ) | (
                                                        df$diff %% 2 < 0.01 &
                                                                df$diff >= 2 &
                                                                df$diff < 3
                                                )
                                                )
                                        massstd <-
                                                apply(df[isoindex, ], 1, function(x)
                                                        min(x[1], x[2]))
                                        massstdmax <-
                                                apply(df[isoindex, ], 1, function(x)
                                                        max(x[1], x[2]))
                                        isomass <-
                                                unique(c(massstd[(massstd %in% massstdmax)], massstdmax))
                                        dfiso <- df[isoindex,]
                                        if (nrow(dfiso) > 0) {
                                                df <-
                                                        df[!(df[, 1] %in% isomass) &
                                                                   !(df[, 2] %in% isomass),]
                                        }
                                        dfdiff <- df
                                        return(
                                                list(
                                                        dfmulti = dfmulti,
                                                        dfiso = dfiso,
                                                        dfdiff = dfdiff,
                                                        solo = NULL
                                                )
                                        )
                                } else {
                                        solo <- cbind(bin, rtg = i)
                                        return(
                                                list(
                                                        dfmulti = NULL,
                                                        dfiso = NULL,
                                                        dfdiff = NULL,
                                                        solo = solo
                                                )
                                        )
                                }
                        }
                }

                rtpmdtemp <-
                        mapply(rtpmd,
                                               split,
                                               as.numeric(names(split)),
                                               SIMPLIFY = F)
                result <- do.call(Map, c(rbind, rtpmdtemp))

                # filter the list get the rt cluster
                list$rtcluster <- rtcluster
                result$dfdiff$diff2 <- round(result$dfdiff$diff, digits)
                # speed up
                pmd <-
                        as.numeric(names(table(result$dfdiff$diff2)[table(result$dfdiff$diff2) > ng]))
                # pmd <- unique(result$dfdiff$diff2)
                idx <- NULL
                for (i in 1:length(pmd)) {
                        l <- length(unique(result$dfdiff[result$dfdiff$diff2 == pmd[i], 'rtg']))
                        idx <-
                                c(idx, ifelse(l > ng, T, F))
                }
                pmd2 <- pmd[idx]
                list$paired <-
                        result$dfdiff[result$dfdiff$diff2 %in% pmd2,]

                if (nrow(result$dfdiff) > 0) {
                        list$pairedindex <-
                                paste(round(list$mz, accuracy), list$rtcluster) %in%
                                paste(c(
                                        round(list$paired$ms1, accuracy),
                                        round(list$paired$ms2,
                                              accuracy)
                                ),
                                c(list$paired$rtg, list$paired$rtg))
                }
                # get the data index by rt groups with high frequences
                # PMD
                list$diffindex <-
                        paste(round(list$mz, accuracy), list$rtcluster) %in%
                        paste(c(
                                round(result$dfdiff$ms1, accuracy),
                                round(result$dfdiff$ms2,
                                      accuracy)
                        ),
                        c(result$dfdiff$rtg, result$dfdiff$rtg))
                list$diff <- result$dfdiff
                # get the data index by rt groups with single ions
                if (!is.null(result$solo)) {
                        list$soloindex <- paste(round(list$mz, accuracy), list$rtcluster) %in%
                                paste(round(result$solo$mz, accuracy), result$solo$rtg)
                        list$solo <- result$solo
                }
                # get the data index by rt groups with isotope ions
                if (!is.null(result$dfiso)) {
                        list$isoindex <- paste(round(list$mz, accuracy), list$rtcluster) %in%
                                paste(c(
                                        round(result$dfiso$ms1, accuracy),
                                        round(result$dfiso$ms2, accuracy)
                                ),
                                c(result$dfiso$rtg, result$dfiso$rtg))
                        list$iso <- result$dfiso
                }
                # get the data index by rt groups with multi charger ions
                if (!is.null(result$dfmulti)) {
                        list$multiindex <- paste(round(list$mz, accuracy), list$rtcluster) %in%
                                paste(c(
                                        round(result$dfmulti$ms1, accuracy),
                                        round(result$dfmulti$ms2,
                                              accuracy)
                                ),
                                c(result$dfmulti$rtg, result$dfmulti$rtg))
                        list$multi <- result$dfmulti
                }

                # show message about std mass
                message(paste(sum(list$pairedindex), "paired masses found "))
                message(
                        paste(
                                length(unique(list$paired$diff2)),
                                "unique within RT clusters high frequency PMD(s) used for further investigation."
                        )
                )
                message(paste(
                        sum(list$isoindex),
                        "isotopologue(s) related paired mass found."
                ))
                message(paste(
                        sum(list$multiindex),
                        "multi-charger(s) related paired mass found."
                ))

                # return results
                return(list)
        }

#' Find the independent ions for each retention time hierarchical clustering based on PMD relationship within each retention time cluster and isotope and return the index of the std data for each retention time cluster.
#' @param list a list from getpaired function
#' @param corcutoff cutoff of the correlation coefficient, default NULL
#' @param digits mass or mass to charge ratio accuracy for pmd, default 2
#' @param accuracy measured mass or mass to charge ratio in digits, default 4
#' @return list with std mass index
#' @examples
#' data(spmeinvivo)
#' pmd <- getpaired(spmeinvivo)
#' std <- getstd(pmd)
#' @seealso \code{\link{getpaired}},\code{\link{getsda}},\code{\link{plotstd}}
#' @export
getstd <- function(list, corcutoff = NULL, digits = 2, accuracy = 4) {
        resultstd2A <- resultstd2B1 <- resultstd2B2 <- resultstd2B3 <- NULL
        # filter high freq ions and find std mass
        resultdiff <- list$paired
        resultiso <- list$iso

        if (!is.null(corcutoff)) {
                resultdiff <- resultdiff[abs(resultdiff$cor) > corcutoff, ]
                resultiso <- resultiso[resultiso$cor > corcutoff, ]
        }
        # filter the mass from mass pairs within retention time
        # group group 1: RT groups with solo peak
        resultstd1 <- NULL
        if (!is.null(list$solo)) {
                resultstd1 <- cbind(list$solo$mz, list$solo$rt,
                                    list$solo$rtg)
        }
        # group 2: RT groups with multiple peaks group 2A: RT
        # groups with multiple peaks while no isotope/paired
        # relationship
        index2A <-
                !(
                        unique(list$rtcluster) %in% unique(resultdiff$rtg) |
                                unique(list$rtcluster) %in% unique(resultiso$rtg)
                ) &
                !(unique(list$rtcluster) %in% unique(list$solo$rtg))
        rtg2A <- unique(list$rtcluster)[index2A]
        message(
                paste(
                        sum(index2A),
                        'group(s) with multiple peaks while no isotope/paired relationship'
                )
        )
        # print(rtg2A)
        if (sum(index2A) > 0) {
                for (i in 1:length(rtg2A)) {
                        mass <- list$mz[list$rtcluster == rtg2A[i]]
                        rt <- list$rt[list$rtcluster == rtg2A[i]]
                        mass <- max(mass)
                        suppressWarnings(resultstdtemp <-
                                                 c(mass, stats::median(rt),
                                                   rtg2A[i]))
                        suppressWarnings(resultstd2A <- rbind(resultstd2A,
                                                              resultstdtemp))
                }
        }
        # group 2B: RT groups with multiple peaks with
        # isotope/paired relationship index2B <-
        # (unique(list$rtcluster) %in%
        # unique(resultdiff$rtg))|(unique(list$rtcluster) %in%
        # unique(resultiso$rtg)) group 2B1: RT groups with
        # multiple peaks with isotope without paired
        # relationship
        index2B1 <-
                !(unique(list$rtcluster) %in% unique(resultdiff$rtg)) &
                unique(list$rtcluster) %in% unique(resultiso$rtg) &
                !(unique(list$rtcluster) %in% unique(list$solo$rtg))
        rtg2B1 <- unique(list$rtcluster)[index2B1]
        # print(rtg2B1)
        message(
                paste(
                        sum(index2B1),
                        'group(s) with multiple peaks with isotope without paired relationship'
                )
        )
        if (sum(index2B1) > 0) {
                for (i in 1:length(rtg2B1)) {
                        # filter the isotope peaks
                        dfiso <-
                                resultiso[resultiso$rtg == rtg2B1[i], ]
                        if (nrow(dfiso) > 0) {
                                massstd <- apply(dfiso, 1, function(x)
                                        min(x[1],
                                            x[2]))
                                massstdmax <-
                                        apply(dfiso, 1, function(x)
                                                max(x[1],
                                                    x[2]))
                                mass <-
                                        unique(massstd[!(massstd %in% massstdmax)])
                                suppressWarnings(
                                        resultstdtemp <- cbind(
                                                mz = ifelse(is.na(mass), NULL, c(mass)),
                                                rt = dfiso$rt,
                                                rtg = dfiso$rtg
                                        )
                                )
                                resultstd2B1 <-
                                        rbind(resultstd2B1, resultstdtemp)
                        } else{
                                resultstd2B1 <- NULL
                        }
                }
        }

        # group 2B2: RT groups with multiple peaks with paired
        # relationship without isotope
        index2B2 <-
                (unique(list$rtcluster) %in% unique(resultdiff$rtg)) &
                !(unique(list$rtcluster) %in% unique(resultiso$rtg)) &
                !(unique(list$rtcluster) %in% unique(list$solo$rtg))
        rtg2B2 <- unique(list$rtcluster)[index2B2]
        # print(rtg2B2)
        message(paste(
                sum(index2B2),
                'group(s) with paired relationship without isotope'
        ))
        if (sum(index2B2) > 0) {
                for (i in 1:length(rtg2B2)) {
                        # filter the paired peaks
                        df <-
                                resultdiff[resultdiff$rtg == rtg2B2[i], ]
                        if (nrow(df) > 0) {
                                mass <- apply(df, 1, function(x)
                                        min(x[1],
                                            x[2]))
                                mass <- unique(mass)
                                suppressWarnings(
                                        resultstdtemp <-
                                                cbind(
                                                        mz = ifelse(is.na(mass), NULL, c(mass)),
                                                        rt = df$rt,
                                                        rtg = df$rtg
                                                )
                                )
                                suppressWarnings(
                                        resultstd2B2 <- rbind(resultstd2B2,
                                                              resultstdtemp)
                                )
                        }
                }
        }

        # group 2B3: RT groups with multiple peaks with paired
        # relationship and isotope
        index2B3 <-
                (unique(list$rtcluster) %in% unique(resultdiff$rtg)) &
                (unique(list$rtcluster) %in% unique(resultiso$rtg)) &
                !(unique(list$rtcluster) %in% unique(list$solo$rtg))
        rtg2B3 <- unique(list$rtcluster)[index2B3]
        # print(rtg2B3)
        message(paste(
                sum(index2B3),
                'group(s) with paired relationship and isotope'
        ))
        if (sum(index2B3) > 0) {
                for (i in 1:length(rtg2B3)) {
                        # filter the isotope peaks
                        dfiso <- resultiso[resultiso$rtg == rtg2B3[i], ]
                        dfpaired <- resultdiff[resultdiff$rtg == rtg2B3[i], ]
                        if (nrow(dfiso) > 0 & nrow(dfpaired) > 0) {
                                # remove peaks with more than one isotopes
                                massstd <- apply(dfiso, 1, function(x)
                                        min(x[1],
                                            x[2]))
                                massstdmax <- apply(dfiso, 1, function(x)
                                        max(x[1],
                                            x[2]))
                                massstd <-
                                        unique(massstd[!(massstd %in% massstdmax)])
                                dis <-
                                        stats::dist(massstd, method = "manhattan")
                                df <- data.frame(
                                        ms1 = massstd[which(lower.tri(dis),
                                                            arr.ind = T)[, 1]],
                                        ms2 = massstd[which(lower.tri(dis),
                                                            arr.ind = T)[, 2]],
                                        diff = round(as.numeric(dis),
                                                     digits)
                                )
                                # remove the adducts
                                if (sum((df$diff %in% dfpaired$diff2)) > 0) {
                                        massstd <- unique(apply(df[df$diff %in%
                                                                           dfpaired$diff2, ], 1, function(x)
                                                                                   min(x[1],
                                                                                       x[2])))
                                        massused <-
                                                unique(c(df$ms1, df$ms2))

                                        massadd <-
                                                unique(c(df$ms1[df$diff %in%
                                                                        dfpaired$diff2], df$ms2[df$diff %in%
                                                                                                        dfpaired$diff2]))
                                        massextra <-
                                                massused[!(massused %in% massadd)]
                                        mass <- c(massextra, massstd)
                                } else {
                                        mass <- massstd
                                }
                                suppressWarnings(
                                        resultstdtemp <- cbind(
                                                mz = ifelse(is.na(mass), NULL, c(mass)),
                                                rt = dfiso$rt,
                                                rtg = dfiso$rtg
                                        )
                                )
                                suppressWarnings(
                                        resultstd2B3 <- rbind(resultstd2B3,
                                                              resultstdtemp)
                                )
                        } else if (nrow(dfiso) > 0) {
                                # remove peaks with more than one peaks
                                massstd <- apply(dfiso, 1, function(x)
                                        min(x[1],
                                            x[2]))
                                massstdmax <- apply(dfiso, 1, function(x)
                                        max(x[1],
                                            x[2]))
                                mass <-
                                        unique(massstd[!(massstd %in% massstdmax)])
                                suppressWarnings(
                                        resultstdtemp <- cbind(
                                                mz = ifelse(is.na(mass), NULL, c(mass)),
                                                rt = dfiso$rt,
                                                rtg = dfiso$rtg
                                        )
                                )
                                suppressWarnings(
                                        resultstd2B3 <- rbind(resultstd2B3,
                                                              resultstdtemp)
                                )
                        } else if (nrow(dfpaired) > 0) {
                                mass <- apply(dfpaired, 1, function(x)
                                        min(x[1],
                                            x[2]))
                                mass <- unique(mass)
                                suppressWarnings(
                                        resultstdtemp <- cbind(
                                                mz = ifelse(is.na(mass), NULL, c(mass)),
                                                rt = dfpaired$rt,
                                                rtg = dfpaired$rtg
                                        )
                                )
                                suppressWarnings(
                                        resultstd2B3 <- rbind(resultstd2B3,
                                                              resultstdtemp)
                                )
                        } else{
                                resultstd2B3 <- NULL
                        }

                }
        }

        # Combine the peaks from rt groups with single ion
        resultstd <- rbind(resultstd1,
                           resultstd2A,
                           resultstd2B3,
                           resultstd2B2,
                           resultstd2B1)
        resultstd <- unique(resultstd)

        colnames(resultstd) <- c("mz", "rt", "rtg")
        resultstd <- as.data.frame(resultstd)

        # return the data
        list$stdmassindex <-
                paste(round(list$mz, accuracy), list$rtcluster) %in%
                paste(round(resultstd$mz, accuracy), resultstd$rtg)
        list$stdmass <- resultstd
        # use correlation to refine peaks within the same retention groups
        if (!is.null(corcutoff)) {
                mzo <- NULL
                for (i in 1:length(unique(list$rtcluster))) {
                        resulttemp <- list$data[list$rtcluster == i & list$stdmassindex, ]
                        mz <-
                                list$mz[list$rtcluster == i & list$stdmassindex]
                        cor2 <- stats::cor(t(resulttemp))
                        df <- data.frame(
                                ms1 = mz[which(lower.tri(cor2),
                                               arr.ind = T)[, 1]],
                                ms2 = mz[which(lower.tri(cor2),
                                               arr.ind = T)[, 2]],
                                cor = cor2[lower.tri(cor2)]
                        )
                        df2 <-
                                apply(df, 1, function(x)
                                        ifelse(abs(x[3]) > corcutoff, x[1], NA))
                        df2 <- unique(stats::na.omit(df2))
                        mz2 <-
                                paste0(round(mz[(round(mz, accuracy) %in% round(df2, accuracy))], accuracy), '@', i)
                        mzo <- c(mzo, mz2)
                }
                list$stdmassindex <-
                        list$stdmassindex &
                        (!(paste0(
                                round(list$mz, accuracy), '@', list$rtcluster
                        ) %in% mzo))
                # show message about std mass
                n <- sum(list$stdmassindex)
                message(paste(n, "std mass found."))
        } else{
                # show message about std mass
                n <- nrow(resultstd)
                message(paste(n, "std mass found."))
        }
        message(paste(
                sum(list$soloindex),
                "retention group(s) have single peaks."
        ))
        return(list)
}

#' Perform structure/reaction directed analysis for peaks list.
#' @param list a list with mzrt profile
#' @param rtcutoff cutoff of the distances in retention time hierarchical clustering analysis, default 10
#' @param freqcutoff cutoff of frequency of PMDs between RT cluster for peaks, default 10
#' @param top top n pmd freqency cutoff when the freqcutoff is too small for large data set, default 50
#' @param corcutoff cutoff of the correlation coefficient, default NULL
#' @param digits mass or mass to charge ratio accuracy for pmd, default 2
#' @param accuracy measured mass or mass to charge ratio in digits, default 4
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
                 freqcutoff = 10,
                 top = 50,
                 corcutoff = NULL,
                 digits = 2,
                 accuracy = 4) {
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
                        rtg <- list$rtcluster[list$stdmassindex]
                } else {
                        mz <- list$mz[list$stdmassindex]
                        rt <- list$rt[list$stdmassindex]
                        data <- list$data[list$stdmassindex, ]
                        rtg <- list$rtcluster[list$stdmassindex]
                }
                # PMD analysis
                # remove isomers
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
                        df <- df[df$rtgdiff > 0, ]
                        df$diff2 <- round(df$diff, digits)
                        # use unique isomers
                        index <-
                                !duplicated(paste0(round(df$ms1, accuracy), round(df$ms2, accuracy)))
                        diff <- df$diff2[index]
                        freq <-
                                table(diff)[order(table(diff), decreasing = T)]
                        if (!is.null(top)) {
                                freq <- utils::head(freq, top)
                                message(
                                        paste(
                                                "Top",
                                                top,
                                                "high frequency PMD groups were remained.",
                                                "\n"
                                        )
                                )
                        }

                        if (sum(df$diff2 == 0) > freqcutoff) {
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
                        if (!is.null(corcutoff)) {
                                list$sda <- list$sda[abs(list$sda$cor) > corcutoff,]
                        }

                } else{
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

                        df <- df[df$rtgdiff > 0, ]
                        df$diff2 <- round(df$diff, digits)

                        # use unique isomers
                        index <-
                                !duplicated(paste0(round(df$ms1, accuracy), round(df$ms2, accuracy)))
                        diff <- df$diff2[index]
                        freq <-
                                table(diff)[order(table(diff), decreasing = T)]
                        if (!is.null(top)) {
                                freq <- utils::head(freq, top)
                                message(
                                        paste(
                                                "Top",
                                                top,
                                                "high frequency PMD groups were remained.",
                                                "\n"
                                        )
                                )
                        }

                        if (sum(df$diff2 == 0) > freqcutoff) {
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


                # show message about std mass
                sub <- names(table(list$sda$diff2))
                n <- length(sub)
                message(paste(n, "groups were found as high frequency PMD group.",
                              "\n"))
                message(paste(sub, "were found as high frequency PMD.",
                              "\n"))
                return(list)
        }

#' GlobalStd algorithm with structure/reaction directed analysis
#' @param list a peaks list with mass to charge, retention time and intensity data
#' @param rtcutoff cutoff of the distances in cluster, default 10
#' @param ng cutoff of global PMD's retention time group numbers
#' @param top top n pmd freqency cutoff when the freqcutoff is too small for large data set, default 50
#' @param corcutoff cutoff of the correlation coefficient, default NULL
#' @param freqcutoff cutoff of frequency of PMDs between RT cluster for independent peaks, default 10
#' @param digits mass or mass to charge ratio accuracy for pmd, default 2
#' @param accuracy measured mass or mass to charge ratio in digits, default 4
#' @return list with GlobalStd algorithm processed data.
#' @examples
#' data(spmeinvivo)
#' re <- globalstd(spmeinvivo)
#' @seealso \code{\link{getpaired}},\code{\link{getstd}},\code{\link{getsda}},\code{\link{plotstd}},\code{\link{plotstdsda}},\code{\link{plotstdrt}}
#' @export
globalstd <- function(list,
                      rtcutoff = 10,
                      ng = 10,
                      corcutoff = NULL,
                      freqcutoff = 10,
                      top = 50,
                      digits = 2, accuracy = 4) {
        list <-
                getpaired(list,
                          rtcutoff = rtcutoff,
                          ng = ng,digits = digits, accuracy = accuracy)
        if (sum(list$pairedindex) > 0) {
                list2 <- getstd(list, corcutoff = corcutoff, digits = digits,accuracy = accuracy)
                list3 <-
                        getsda(
                                list2,
                                rtcutoff = rtcutoff,
                                freqcutoff = freqcutoff,
                                corcutoff = corcutoff,
                                digits = digits
                        )
        } else{
                message('no paired relationship, directly go to structure directed analysis.')
                list3 <-
                        getsda(
                                list,
                                rtcutoff = rtcutoff,
                                freqcutoff = freqcutoff,
                                corcutoff = corcutoff,
                                top = top,
                                digits = digits
                        )
        }

        return(list3)
}

#' Get Pseudo-Spectrum as peaks cluster based on correlation analysis.
#' @param list a list with peaks intensity
#' @param corcutoff cutoff of the correlation coefficient, default 0.9
#' @param rtcutoff cutoff of the distances in cluster, default 10
#' @param accuracy measured mass or mass to charge ratio in digits, default 4
#' @return list with Pseudo-Spectrum index
#' @examples
#' data(spmeinvivo)
#' cluster <- getcorcluster(spmeinvivo)
#' @export
getcorcluster <- function(list,
                          corcutoff = 0.9,
                          rtcutoff = 10,
                          accuracy = 4) {
        mz <- list$mz
        rt <- list$rt

        dis <- stats::dist(list$rt, method = "manhattan")
        fit <- stats::hclust(dis)
        rtcluster <- stats::cutree(fit, h = rtcutoff)
        n <- length(unique(rtcluster))
        message(paste(n, "retention time cluster found."))
        cluster <- mzs <- mzo <- NULL
        data <- list$data

        for (i in 1:length(unique(rtcluster))) {
                # find the mass within RT
                rtxi <- rt[rtcluster == i]
                bin <- data[rtcluster == i, ]
                if (is.matrix(bin)) {
                        msdata <- apply(bin, 1, mean)
                } else {
                        msdata <- mean(bin)
                }
                mzt <- mz[rtcluster == i]
                cor2 <- stats::cor(t(bin))
                df <- data.frame(ms1 = mzt[which(lower.tri(cor2),
                                                 arr.ind = T)[, 1]],
                                 ms2 = mzt[which(lower.tri(cor2),
                                                 arr.ind = T)[, 2]],
                                 cor = cor2[lower.tri(cor2)])
                df2 <-
                        apply(df, 1, function(x)
                                ifelse(abs(x[3]) > corcutoff, x[1], NA))
                df2 <- unique(stats::na.omit(df2))
                mzi <- mzt[!(mzt %in% df2)]
                clustert <- NULL
                for (j in 1:length(mzi)) {
                        mzic <- df$ms1[round(df$ms2, accuracy) == round(mzi[j], accuracy)]
                        mzic <- unique(c(mzic, mzi[j]))
                        ins <-
                                msdata[round(mzt, accuracy) %in% round(mzic, accuracy)]
                        tdf <- cbind.data.frame(mzic, j, i, ins)
                        clustert <- rbind(clustert, tdf)
                }
                cluster <- rbind(cluster, clustert)
                mz1 <-
                        stats::aggregate(clustert, by = list(clustert$j), max)
                mz1 <- unique(mz1[, 2])
                mz2 <-
                        paste0(round(mz[(round(mz, accuracy) %in% round(df2, accuracy))], accuracy), '@', i)
                mzo <- c(mzo, mz2)
                mzs <- c(mzs, mz1)
        }
        list$stdmassindex <-
                !(paste0(round(list$mz, accuracy), '@', rtcluster) %in% mzo)
        if (!is.null(mzs)) {
                list$stdmassindex2 <- round(list$mz, accuracy) %in% round(mzs, accuracy)
        }
        list$rtcluster <- cluster
        return(list)
}

#' Get Pseudo-Spectrum as peaks cluster based on pmd analysis.
#' @param list a list from getstd function
#' @param corcutoff cutoff of the correlation coefficient, default NULL
#' @param accuracy measured mass or mass to charge ratio in digits, default 4
#' @return list with Pseudo-Spectrum index
#' @examples
#' data(spmeinvivo)
#' re <- getpaired(spmeinvivo)
#' re <- getstd(re)
#' cluster <- getcluster(re)
#' @seealso \code{\link{getpaired}},\code{\link{getstd}},\code{\link{plotstd}}
#' @export
getcluster <- function(list, corcutoff = NULL, accuracy = 4) {
        mz <- list$mz[list$stdmassindex]
        rt <- list$rt[list$stdmassindex]
        rtg <- list$rtcluster[list$stdmassindex]
        if (is.null(list$data)) {
                message('You need intensity data to use corcutoff and export pseudospectra')
                msdata <- NULL
        } else{
                data <- list$data
                if (is.matrix(data)) {
                        msdata <- apply(data, 1, mean)
                } else {
                        msdata <- mean(data)
                }
        }

        stdg <- rep('stdgroup', length(list$rtcluster))
        # filter high freq ions and find std mass
        resultdiff <- list$paired
        resultiso <- list$iso
        resultmulti <- list$multi

        if (!is.null(corcutoff)) {
                resultdiff <- resultdiff[abs(resultdiff$cor) > corcutoff, ]
                resultiso <- resultiso[resultiso$cor > corcutoff, ]
                resultmulti <-
                        resultmulti[resultmulti$cor > corcutoff,]
        }

        # multi-charger
        index1 <-
                paste0(round(resultmulti$ms1, accuracy), '@', resultmulti$rtg)
        index2 <-
                paste0(round(resultmulti$ms2, accuracy), '@', resultmulti$rtg)
        # isotope
        index3 <- paste0(round(resultiso$ms1, accuracy), '@', resultiso$rtg)
        index4 <- paste0(round(resultiso$ms2, accuracy), '@', resultiso$rtg)
        # diff
        index5 <- paste0(round(resultdiff$ms1, accuracy), '@', resultdiff$rtg)
        index6 <- paste0(round(resultdiff$ms2, accuracy), '@', resultdiff$rtg)

        index00 <- paste0(round(list$mz, accuracy), '@', list$rtcluster)


        cluster <- mzs <- NULL
        for (i in 1:sum(list$stdmassindex)) {
                mzt <- mz[i]
                rtgt <- rtg[i]
                indexstd <- paste0(round(mzt, accuracy), '@', rtgt)

                multiover <-
                        unique(c(resultmulti$ms2[index1 %in% indexstd], resultmulti$ms1[index2 %in% indexstd]))
                isoover <-
                        unique(c(resultiso$ms2[index3 %in% indexstd], resultiso$ms1[index4 %in% indexstd]))
                diffover <-
                        unique(c(resultdiff$ms2[index5 %in% indexstd], resultdiff$ms1[index6 %in% indexstd]))
                stdmassg <- c(mzt, multiover, isoover, diffover)
                stdg[round(list$mz, accuracy) %in% round(stdmassg, accuracy) &
                             list$rtcluster == rtgt] <-
                        paste0(stdg[round(list$mz, accuracy) %in% round(stdmassg, accuracy) &
                                            list$rtcluster == rtgt], '@', i)
                if (!is.null(msdata)) {
                        index <- paste0(round(stdmassg, accuracy), '@', rtgt)
                        ins <- msdata[unique(index00) %in% index]
                        tdf <- cbind.data.frame(stdmassg, i, rtgt, ins)
                        mzst <- index[which.max(tdf$ins)]
                        mzs <- c(mzs, mzst)
                } else{
                        tdf <- cbind.data.frame(stdmassg, i, rtgt)
                }

                cluster <- rbind.data.frame(cluster, tdf)
        }
        list$stdg <- stdg
        list$cluster <- cluster
        if (!is.null(mzs)) {
                list$stdmassindex2 <-
                        paste0(round(list$mz, accuracy), '@', list$rtcluster) %in% mzs
        }
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
                getmass <- function(data) {
                        if (grepl('-', data)) {
                                name <- unlist(strsplit(data, '-'))
                                iso1 <-
                                        rcdk::get.isotopes.pattern(rcdk::get.formula(name[1]))
                                iso2 <-
                                        rcdk::get.isotopes.pattern(rcdk::get.formula(name[2]))
                                cus <-
                                        as.numeric(iso1[max(iso1[, 2]), 1]) - as.numeric(iso2[max(iso2[, 2]), 1])
                        } else{
                                iso <- rcdk::get.isotopes.pattern(rcdk::get.formula(data))
                                cus <-
                                        as.numeric(iso[max(iso[, 2]), 1])
                        }
                        return(cus)
                }
                if (is.null(formula)) {
                        dis <- stats::dist(mz, method = "manhattan")
                } else{
                        mz <- unlist(Map(getmass, formula))
                        dis <- stats::dist(mz, method = "manhattan")

                }
                mz <- unique(mz)

                df <- cbind.data.frame(
                        ms1 = mz[which(lower.tri(dis), arr.ind = T)[, 1]],
                        ms2 = mz[which(lower.tri(dis), arr.ind = T)[, 2]],
                        diff = as.numeric(dis),
                        diff2 = round(as.numeric(dis), digits = digits)
                )
                freq <-
                        table(df$diff2)[order(table(df$diff2), decreasing = T)]
                message(paste(length(freq), 'pmd found.'))
                if (!is.null(top)) {
                        freq <- utils::head(freq, top)
                }
                sda <- df[(df$diff2 %in% c(as.numeric(names(freq[freq >= freqcutoff])))), ]
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
                result <- mapply(rtpmd, split, as.numeric(names(split)))
                rownames(result) <- mz
                return(as.data.frame(result))
        }
#' Get pmd for specific reaction
#' @param list a list with mzrt profile
#' @param pmd a specific paired mass distances
#' @param rtcutoff cutoff of the distances in retention time hierarchical clustering analysis, default 10
#' @param digits mass or mass to charge ratio accuracy for pmd, default 2
#' @param accuracy measured mass or mass to charge ratio in digits, default 4
#' @return list with paired peaks for specific pmd.
#' @examples
#' data(spmeinvivo)
#' pmd <- getpmd(spmeinvivo,pmd=15.99)
#' @seealso \code{\link{getpaired}},\code{\link{getstd}},\code{\link{getsda}},\code{\link{getrda}}
#' @export
getpmd <- function(list, pmd, rtcutoff = 10, digits = 2, accuracy = 4) {
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

        df$diff2 <- round(df$diff, digits)

        df <- df[df$rtgdiff > 0 & df$diff2 == pmd, ]
        ms1 <- ifelse(df$ms1>df$ms2,df$ms1,df$ms2)
        ms2 <- ifelse(df$ms1>df$ms2,df$ms2,df$ms1)
        rtg1 <- ifelse(df$ms1>df$ms2,df$rtg1,df$rtg2)
        rtg2 <- ifelse(df$ms1>df$ms2,df$rtg2,df$rtg1)
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
#' @param mass a specific mass for known compound or a vector of masses
#' @param accuracy measured mass or mass to charge ratio in digits, default 4
#' @param ... other parameters for getpmd
#' @return a list with mzrt profile and reaction chain dataframe
#' @examples
#' data(spmeinvivo)
#' # check metabolites of C18H39NO
#' pmd <- getchain(spmeinvivo,diff = c(2.02,14.02,15.99,58.04,13.98),mass = 286.3101)
#' @export
getchain <- function(list,diff, mass, accuracy = 4,...){
        sda <- getpmd(list, unique(diff)[1])$pmd
        for(i in 2:length(unique(diff))){
                masst <- getpmd(list, unique(diff)[i])
                sda <- rbind.data.frame(sda,masst$pmd)
        }
        seed <- NULL
        ms1 <- round(sda$ms1,digits = accuracy)
        ms2 <- round(sda$ms2,digits = accuracy)
        if(length(mass)==1){
                mass <- round(mass,accuracy)
                sdat <- unique(c(mass,ms2[ms1 %in% mass],ms1[ms2 %in% mass]))
                while(!identical(sdat,seed)){
                        seed <- sdat
                        sdat <- unique(sdat,c(ms2[ms1 %in% sdat],ms1[ms2 %in% sdat]))
                }
                list$sdac <- sda[ms1 %in% sdat|ms2 %in% sdat , ]
                return(list)
        }else if(length(mass)==0){
                warning('No mass input and all mass in the list will be used for reaction chain construction!')
                sdac <- NULL
                mass <- round(list$mz,accuracy)
                for(i in 1:length(mass)){
                        sdat <- unique(c(mass[i],ms2[ms1 %in% mass[i]],ms1[ms2 %in% mass[i]]))
                        if(length(sdat)!=1){
                                while(!identical(sdat,seed)){
                                        seed <- sdat
                                        sdat <- unique(c(sdat,ms2[ms1 %in% sdat],ms1[ms2 %in% sdat]))
                                }
                                sdact <- sda[ms1 %in% sdat|ms2 %in% sdat , ]
                                sdact$mass <- mass[i]
                                sdac <- rbind.data.frame(sdac,sdact)
                        }
                }
                list$sdac <- sdac[!duplicated(sdac),]
                return(list)
        }else{
                sdac <- NULL
                mass <- round(mass,accuracy)
                for(i in 1:length(mass)){
                        sdat <- unique(c(mass[i],ms2[ms1 %in% mass[i]],ms1[ms2 %in% mass[i]]))
                        if(length(sdat)!=1){
                        while(!identical(sdat,seed)){
                                seed <- sdat
                                sdat <- unique(c(sdat,ms2[ms1 %in% sdat],ms1[ms2 %in% sdat]))
                        }
                        sdact <- sda[ms1 %in% sdat|ms2 %in% sdat, ]
                        sdact$mass <- mass[i]
                        sdac <- rbind.data.frame(sdac,sdact)
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
#' @param ratiocv ratio cv cutoff for quantitative paired peaks, default 30
#' @return list with quantitative paired peaks.
#' @examples
#' data(spmeinvivo)
#' pmd <- getreact(spmeinvivo,pmd=15.99)
#' @seealso \code{\link{getpaired}},\code{\link{getstd}},\code{\link{getsda}},\code{\link{getrda}},\code{\link{getpmd}},
#' @export
getreact <- function(list, pmd, rtcutoff = 10, digits = 2, accuracy = 4, ratiocv = 30){
        p <- pmd::getpmd(list,pmd=pmd,rtcutoff = rtcutoff,digits = digits, accuracy = accuracy)
        list <- enviGCMS::getfilter(p,p$pmdindex)
        data <- list$data
        pmd <- list$pmd
        if(NCOL(list$group)>1){
                group <- apply(list$group, 1 , paste0, collapse = "")
                nlv <- unique(group)
        }else{
                group <- c(t(list$group))
                nlv <- unique(group)
        }

        getr <- function(v,nlv){
                ratio <- NULL
                for(i in 1:length(nlv)){
                        ratio1 <- sum(data[list$mz%in%v[1]&list$rt%in%v[4],grepl(nlv[i],group)])
                        ratio2 <- sum(data[list$mz%in%v[2]&list$rt%in%v[5],grepl(nlv[i],group)])
                        ratioi <- ratio1/ratio2
                        ratio <- c(ratio,ratioi)
                }
                rsd <- stats::sd(ratio,na.rm = T)/mean(ratio,na.rm = T)*100
                return(rsd)
        }
        list$pmd$r <- apply(pmd,1,getr,nlv)
        list$pmd <- list$pmd[list$pmd$r<ratiocv,]
        idx <- paste(list$mz,list$rt)
        idx2 <- unique(c(paste(list$pmd$ms1,list$pmd$rt1),paste(list$pmd$ms2,list$pmd$rt2)))
        list <- enviGCMS::getfilter(list,idx%in%idx2)
        return(list)
}

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
                             rt %in% x] <- sample(z, sum(rtcluster == i &
                                                                 rt %in% x), replace = T)
        }
        return(inji)
}
