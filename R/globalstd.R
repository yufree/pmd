#' Filter ions/peaks based on retention time hierarchical clustering, paired mass distances(PMD) and PMD frequency analysis.
#' @param list a list with mzrt profile
#' @param rtcutoff cutoff of the distances in retention time hierarchical clustering analysis, default 10
#' @param ng cutoff of global PMD's retention time group numbers, default 10. If ng = NULL, 20 percent of RT cluster will be used as ng.
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
                 ng = 10,
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
                                                stats::cor(t(bin[, -c(1, 2)]))

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
                                        # higher is C13
                                        # massstd <-
                                        #         apply(df[isoindex,], 1, function(x)
                                        #                 min(x[1], x[2]))
                                        massstdmax <-
                                                apply(df[isoindex,], 1, function(x)
                                                        max(x[1], x[2]))
                                        isomass <-
                                                unique(massstdmax)
                                        dfiso <- df[isoindex, ]
                                        if (nrow(dfiso) > 0) {
                                                df <-
                                                        df[!(df[, 1] %in% isomass) &
                                                                   !(df[, 2] %in% isomass), ]
                                        }
                                        # remove multi chargers
                                        multiindex <-
                                         (round(df$diff %% 1, 1) == 0.5)
                                        dfmulti <- df[multiindex, ]
                                        mass <-
                                         unique(df[multiindex, 1], df[multiindex, 2])
                                        # From HMDB lowest mass with 0.5 is 394.4539, remove those ions related pmd
                                        multimass <-
                                         mass[round(mass %% 1, 1) == 0.5 & mass<350]

                                        if (nrow(dfmulti) > 0) {
                                                 df <-
                                                         df[!(df[, 1] %in% multimass) & !(df[, 2] %in% multimass),]
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
                                        # higher is C13
                                        # massstd <-
                                        #         apply(df[isoindex,], 1, function(x)
                                        #                 min(x[1], x[2]))
                                        massstdmax <-
                                                apply(df[isoindex,], 1, function(x)
                                                        max(x[1], x[2]))
                                        isomass <-
                                                unique(massstdmax)
                                        dfiso <- df[isoindex, ]
                                        if (nrow(dfiso) > 0) {
                                                df <-
                                                        df[!(df[, 1] %in% isomass) &
                                                                   !(df[, 2] %in% isomass), ]
                                        }
                                        # remove multi chargers
                                        multiindex <-
                                                (round(df$diff %% 1, 1) == 0.5)
                                        dfmulti <- df[multiindex, ]
                                        mass <-
                                                unique(df[multiindex, 1], df[multiindex, 2])
                                        # From HMDB lowest mass with 0.5 is 394.4539, remove those ions related pmd
                                        multimass <-
                                                mass[round(mass %% 1, 1) == 0.5 & mass<350]

                                        if (nrow(dfmulti) > 0) {
                                                df <-
                                                        df[!(df[, 1] %in% multimass) & !(df[, 2] %in% multimass),]
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
                result$dfdiff$diff2 <-
                        round(result$dfdiff$diff, digits)
                # speed up
                pmd <-
                        as.numeric(names(table(result$dfdiff$diff2)[table(result$dfdiff$diff2) > ng]))
                # pmd <- unique(result$dfdiff$diff2)
                idx <- NULL
                for (i in 1:length(pmd)) {
                        l <-
                                length(unique(result$dfdiff[result$dfdiff$diff2 == pmd[i], 'rtg']))
                        idx <-
                                c(idx, ifelse(l > ng, T, F))
                }
                pmd2 <- pmd[idx]
                list$paired <-
                        result$dfdiff[result$dfdiff$diff2 %in% pmd2, ]

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
                        list$soloindex <-
                                paste(round(list$mz, accuracy), list$rtcluster) %in%
                                paste(round(result$solo$mz, accuracy),
                                      result$solo$rtg)
                        list$solo <- result$solo
                }

                # get the data index by rt groups with isotope ions
                if (!is.null(result$dfiso)) {
                        list$isoindex <-
                                paste(round(list$mz, accuracy), list$rtcluster) %in%
                                paste(c(
                                        round(result$dfiso$ms1, accuracy),
                                        round(result$dfiso$ms2, accuracy)
                                ),
                                c(result$dfiso$rtg, result$dfiso$rtg))
                        list$iso <- result$dfiso
                }
                # get the data index by rt groups with multi charger ions
                if (!is.null(result$dfmulti)) {
                        list$multiindex <-
                                paste(round(list$mz, accuracy), list$rtcluster) %in%
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
getstd <-
        function(list,
                 corcutoff = NULL,
                 digits = 2,
                 accuracy = 4) {
                resultstd2A <- resultstd2B1 <- resultstd2B2 <- resultstd2B3 <- NULL
                # filter high freq ions and find std mass
                resultdiff <- list$paired
                resultiso <- list$iso

                if (!is.null(corcutoff)) {
                        resultdiff <- resultdiff[abs(resultdiff$cor) >= corcutoff,]
                        resultiso <- resultiso[resultiso$cor >= corcutoff,]
                }
                # filter the mass from mass pairs within retention time
                # group group 1: RT groups with solo peak
                resultstd1 <- NULL
                if (!is.null(list$solo)) {
                        resultstd1 <- cbind(list$solo$mz, list$solo$rt,
                                            list$solo$rtg)
                }
                message(paste(
                        c(
                                sum(list$soloindex),
                                "retention group(s) have single peaks.",
                                list$rtcluster[list$soloindex]
                        ),
                        collapse = ' '
                ))
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
                message(paste(
                        c(
                                sum(index2A),
                                'group(s) with multiple peaks while no isotope/paired relationship',
                                rtg2A
                        ),
                        collapse = " "
                ))
                # print(rtg2A)
                if (sum(index2A) > 0) {
                        for (i in 1:length(rtg2A)) {
                                mass <- list$mz[list$rtcluster == rtg2A[i]]
                                rt <- list$rt[list$rtcluster == rtg2A[i]]
                                mass <- max(mass)
                                suppressWarnings(resultstdtemp <-
                                                         c(mass, stats::median(rt),
                                                           rtg2A[i]))
                                suppressWarnings(resultstd2A <-
                                                         rbind(resultstd2A,
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
                message(paste(
                        c(
                                sum(index2B1),
                                'group(s) with multiple peaks with isotope without paired relationship',
                                rtg2B1
                        ),
                        collapse = " "
                ))
                if (sum(index2B1) > 0) {
                        for (i in 1:length(rtg2B1)) {
                                # filter the isotope peaks
                                dfiso <-
                                        resultiso[resultiso$rtg == rtg2B1[i],]
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

                message(paste(
                        c(
                                sum(index2B2),
                                'group(s) with paired relationship without isotope',
                                rtg2B2
                        ),
                        collapse = ' '
                ))
                # print(rtg2B2)
                if (sum(index2B2) > 0) {
                        for (i in 1:length(rtg2B2)) {
                                # filter the paired peaks
                                df <-
                                        resultdiff[resultdiff$rtg == rtg2B2[i],]
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

                message(paste(
                        c(
                                sum(index2B3),
                                'group(s) with paired relationship and isotope',
                                rtg2B3
                        ),
                        collapse = ' '
                ))
                # print(rtg2B3)
                if (sum(index2B3) > 0) {
                        for (i in 1:length(rtg2B3)) {
                                # filter the isotope peaks
                                dfiso <-
                                        resultiso[resultiso$rtg == rtg2B3[i],]
                                dfpaired <-
                                        resultdiff[resultdiff$rtg == rtg2B3[i],]
                                if (nrow(dfiso) > 0 & nrow(dfpaired) > 0) {
                                        # remove peaks with more than one isotopes
                                        massstd <-
                                                apply(dfiso, 1, function(x)
                                                        min(x[1],
                                                            x[2]))
                                        massstdmax <-
                                                apply(dfiso, 1, function(x)
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
                                                                                   dfpaired$diff2,], 1, function(x)
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
                                                mass <-
                                                        c(massextra, massstd)
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
                                        massstd <-
                                                apply(dfiso, 1, function(x)
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
                                resulttemp <- list$data[list$rtcluster == i & list$stdmassindex,]
                                mz <-
                                        list$mz[list$rtcluster == i &
                                                        list$stdmassindex]
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
                                                ifelse(abs(x[3]) >= corcutoff, x[1], NA))
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

                return(list)
        }

#' GlobalStd algorithm with structure/reaction directed analysis
#' @param list a peaks list with mass to charge, retention time and intensity data
#' @param rtcutoff cutoff of the distances in cluster, default 10
#' @param ng cutoff of global PMD's retention time group numbers, default 10. If ng = NULL, 20 percent of RT cluster will be used as ng.
#' @param corcutoff cutoff of the correlation coefficient, default NULL
#' @param digits mass or mass to charge ratio accuracy for pmd, default 2
#' @param accuracy measured mass or mass to charge ratio in digits, default 4
#' @param freqcutoff pmd freqency cutoff for structures or reactions, default NULL. This cutoff will be found by PMD network analysis when it is NULL.
#' @param sda logical, option to perform structure/reaction directed analysis, default T.
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
                      digits = 2,
                      accuracy = 4,
                      freqcutoff = NULL,
                      sda = T) {
        list <-
                getpaired(
                        list,
                        rtcutoff = rtcutoff,
                        ng = ng,
                        digits = digits,
                        accuracy = accuracy
                )
        if (sum(list$pairedindex) > 0) {
                list2 <-
                        getstd(
                                list,
                                corcutoff = corcutoff,
                                digits = digits,
                                accuracy = accuracy
                        )
                if(sda){
                list3 <-
                        getsda(
                                list2,
                                rtcutoff = rtcutoff,
                                corcutoff = corcutoff,
                                digits = digits,
                                freqcutoff = freqcutoff
                        )
                }
        } else{
                if(sda){
                        message('no paired relationship, directly go to structure directed analysis.')
                        list3 <-
                                getsda(
                                        list,
                                        rtcutoff = rtcutoff,
                                        corcutoff = corcutoff,
                                        digits = digits,
                                        freqcutoff = freqcutoff
                                )
                }else{
                        message('no paired relationship found.')
                }

        }
        if(sda){
                return(list3)
        }else{
                return(list2)
        }
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
                if(sum(rtcluster==i)>1){
                        bin <- data[rtcluster == i,]
                        if (is.matrix(bin)) {
                                msdata <- apply(bin, 1, mean)
                        } else {
                                msdata <- mean(bin)
                        }
                        mzt <- round(mz[rtcluster == i],digits = accuracy)
                        cor2 <- stats::cor(t(bin))
                        df <- data.frame(ms1 = mzt[which(lower.tri(cor2),
                                                         arr.ind = T)[, 1]],
                                         ms2 = mzt[which(lower.tri(cor2),
                                                         arr.ind = T)[, 2]],
                                         cor = cor2[lower.tri(cor2)])
                        dfc <- df[df$cor>=corcutoff,]
                        # select larger ions
                        df2 <- apply(df,1,function(x) ifelse(abs(x[3]) >= corcutoff, x[1], NA))
                        df2 <- unique(stats::na.omit(df2))
                        if(nrow(dfc)>0){
                                mztn <- unique(c(dfc$ms1,dfc$ms2))
                                nt <- igraph::graph_from_data_frame(dfc)
                                x <- igraph::components(nt)$membership
                                ins <-
                                        msdata[match(as.numeric(igraph::V(nt)$name),mzt)]
                                mztnl <- stats::aggregate(ins,by=list(x),max)
                                mzc <- mzt[ins %in% mztnl$x]
                                mzi <-  mzt[!(mzt %in% mztn)]
                                clustert <- NULL
                                for(j in 1:length(unique(x))){
                                        mzic <- round(as.numeric(igraph::V(nt)$name), accuracy)[x==j]
                                        inst <- ins[x==j]
                                        tdf <- cbind.data.frame(mz=mzic, i=j, rtgt=i, ins=inst)
                                        clustert <- rbind(clustert, tdf)
                                }
                                if(length(mzi)>0){
                                        clusterp <- cbind.data.frame(mz=mzi, i=(length(unique(x))+1):(length(unique(x))+length(mzi)), rtgt=i, ins=msdata[match(mzi,mzt)])
                                        cluster <- rbind(cluster, clustert,clusterp)
                                }else{
                                        cluster <- rbind(cluster, clustert)
                                }
                                mzo <- c(mzo, paste0(df2, '@', i))
                                mzs <- c(mzs, mzc,mzi)
                        }else{
                                clustert <- cbind.data.frame(mz=mzt, i=1:length(mzt), rtgt=i, ins=msdata)
                                cluster <- rbind(cluster,clustert)
                                mzo <- c(mzo, paste0(mzt, '@', i))
                                mzs <- c(mzs, mzt)

                        }
                }else{
                        mzt <- round(mz[rtcluster == i],digits = accuracy)
                        clustert <- cbind.data.frame(mz=mzt, i=1:length(mz[rtcluster==i]), rtgt=i, ins=mean(data[rtcluster==i,]))
                        cluster <- rbind(cluster,clustert)
                        mzo <- c(mzo, paste0(mzt, '@', i))
                        mzs <- c(mzs, mzt)
                }

        }
        list$stdmassindex <-
                !(paste0(round(list$mz, accuracy), '@', rtcluster) %in% mzo)
        if (!is.null(mzs)) {
                list$stdmassindex2 <-
                        round(list$mz, accuracy) %in% round(mzs, accuracy)
        }
        cluster$largei <- paste0(cluster$i,'@',cluster$rtgt)
        cluster <- cluster[match(round(mz,accuracy),cluster$mz),]
        list$rtcluster <- rtcluster
        list$cluster <- cluster
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
getcluster <- function(list,
                       corcutoff = NULL,
                       accuracy = 4) {
        mz <- list$mz[list$stdmassindex]
        rt <- list$rt[list$stdmassindex]
        rtg <- list$rtcluster[list$stdmassindex]
        if (is.null(list$data)) {
                message('You need intensity data to use corcutoff and export pseudospectra')
                msdata <- NULL
        } else{
                data <- list$data
                if (is.matrix(data)) {
                        msdata <- apply(data, 1, sum)
                } else {
                        msdata <-sum(data)
                }
        }

        stdg <- rep('stdgroup', length(list$rtcluster))
        # filter high freq ions and find std mass
        resultdiff <- list$paired
        resultiso <- list$iso
        resultmulti <- list$multi

        if (!is.null(corcutoff)) {
                resultdiff <- resultdiff[abs(resultdiff$cor) >= corcutoff,]
                resultiso <- resultiso[resultiso$cor >= corcutoff,]
                resultmulti <-
                        resultmulti[resultmulti$cor >= corcutoff, ]
        }

        # multi-charger
        index1 <-
                paste0(round(resultmulti$ms1, accuracy), '@', resultmulti$rtg)
        index2 <-
                paste0(round(resultmulti$ms2, accuracy), '@', resultmulti$rtg)

        # isotope
        index3 <-
                paste0(round(resultiso$ms1, accuracy), '@', resultiso$rtg)
        index4 <-
                paste0(round(resultiso$ms2, accuracy), '@', resultiso$rtg)
        # diff
        index5 <-
                paste0(round(resultdiff$ms1, accuracy), '@', resultdiff$rtg)
        index6 <-
                paste0(round(resultdiff$ms2, accuracy), '@', resultdiff$rtg)

        index00 <-
                paste0(round(list$mz, accuracy), '@', list$rtcluster)
        msdata <- msdata[!duplicated(index00)]
        index000 <- unique(index00)

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
                mzx <- list$mz[list$mz %in% stdmassg]
                if (!is.null(msdata)) {
                        index <- paste0(round(mzx, accuracy), '@', rtgt)
                        ins <- msdata[index000 %in% unique(index)]

                        tdf <-
                                cbind.data.frame(mz = mzx[!duplicated(index)], i, rtgt, ins)
                } else{
                        tdf <- cbind.data.frame(stdmassg, i, rtgt)
                }

                cluster <- rbind.data.frame(cluster, tdf)

        }
        mergegroup <- function(temp){
                t <- any(table(temp$mz)>1)
                if(t){
                        cluster2 <- temp[!duplicated(temp$mz),]
                        ti<-igraph::components(igraph::graph_from_data_frame(temp))$membership[1:length(cluster2$mz)]
                        cluster2$largei <- paste0(ti,'@',cluster2$rtgt)
                        return(cluster2)
                }else{
                        temp$largei <-  paste0(temp$i,'@',temp$rtgt)
                        return(temp)
                }
        }
        cl <- split(cluster,cluster$rtgt)
        cluster <- lapply(cl, mergegroup)
        list$cluster <- do.call("rbind", cluster)
        if (!is.null(msdata)) {
                for(i in unique(list$cluster$largei)){
                        t <- list$cluster[list$cluster$largei==i,]
                        mzst <- paste0(unique(round(t$mz[which.max(t$ins)],accuracy)),'@',t$rtgt[1])
                        mzs <- c(mzs, mzst)
                }
        }
        if (!is.null(mzs)) {
                list$stdmassindex2 <-
                        paste0(round(list$mz, accuracy), '@', list$rtcluster) %in% mzs
        }
        return(list)
}

#' Compare matrices using PCA similarity factor
#'
#' @param x Matrix with sample in column and features in row
#' @param y Matrix is compared to x.
#' @param dim number of retained dimensions in the comparison. Defaults to all.
#' @return Ratio of projected variance to total variance
#' @references Singhal, A. and Seborg, D. E. (2005), Clustering multivariate time-series data. J. Chemometrics, 19: 427-438. doi: 10.1002/cem.945
#' @author Edgar Zanella Alvarenga
#' @export
#' @examples
#' c1 <- matrix(rnorm(16),nrow=4)
#' c2 <- matrix(rnorm(16),nrow=4)
#' pcasf(c1, c2)
#'
pcasf <- function(x, y, dim = NULL) {
        cov.x <- stats::cov(x)
        cov.y <- stats::cov(y)

        if (is.null(dim)){
                dim = dim(cov.x)[1]
        }

        eg.x <- eigen(cov.x)
        eg.y <- eigen(cov.y)
        eg.x.values <- eg.x$values[1:dim]
        eg.y.values <- eg.y$values[1:dim]
        eg.x.vectors <- eg.x$vectors[,1:dim]
        eg.y.vectors <- eg.y$vectors[,1:dim]

        total_var <- eg.x.values %*% eg.y.values

        return (c(pcasf = sum((eg.x.values %o% eg.y.values) * ((t(eg.x.vectors) %*% (eg.y.vectors))**2))/total_var))
}
