#' Filter ions/peaks based on retention time hierarchical clustering, paired mass distances(PMD) and PMD frequency analysis.
#' @param list a peaks list with mass to charge, retention time and intensity data
#' @param rtcutoff cutoff of the distances in retention time hierarchical clustering analysis, default 10
#' @param ng cutoff of global PMD's retention time group numbers, If ng = NULL, 20 percent of RT cluster will be used as ng, default NULL.
#' @param digits mass or mass to charge ratio accuracy for pmd, default 2
#' @param accuracy measured mass or mass to charge ratio in digits, default 4
#' @param corcutoff cutoff of the correlation coefficient, 0.6 is suggested, default NULL
#' @return list with tentative isotope, multi-chargers, adducts, and neutral loss peaks' index, retention time clusters.
#' @examples
#' data(spmeinvivo)
#' pmd <- getpaired(spmeinvivo)
#' @seealso \code{\link{getstd}},\code{\link{getsda}},\code{\link{plotpaired}}
#' @export
getpaired <- function(list, rtcutoff = 10, ng = NULL, digits = 2,
                      accuracy = 4, corcutoff = NULL) {
        # Calculate retention time clusters
        rt_clusters <- function(rt, cutoff) {
                dis <- stats::dist(rt, method = "manhattan")
                fit <- stats::hclust(dis)
                stats::cutree(fit, h = cutoff)
        }

        # Core PMD analysis function
        analyze_pmd <- function(cluster_data, rt_group, didx, corcutoff) {
                # Handle empty clusters
                if (is.null(cluster_data) || nrow(cluster_data) < 1) return(list())

                median_rt <- stats::median(cluster_data$rt)
                results <- list(
                        iso = NULL,
                        multiiso = NULL,
                        multi = NULL,
                        diff = NULL,
                        solo = NULL,
                        iso_mz = numeric(0),
                        multiiso_mz = numeric(0),
                        rtg = rt_group
                )

                if (nrow(cluster_data) == 1) {
                        results$solo <- cbind(cluster_data[, c("mz", "rt")],
                                              rtg = rt_group)
                        return(results)
                } else {
                        # Calculate pairwise differences
                        mz_pairs <- t(utils::combn(cluster_data$mz, 2))
                        diffs <- abs(mz_pairs[, 2] - mz_pairs[, 1])

                        # Create base dataframe
                        df <- data.frame(
                                ms1 = pmin(mz_pairs[, 1], mz_pairs[, 2]),
                                ms2 = pmax(mz_pairs[, 1], mz_pairs[, 2]),
                                diff = diffs,
                                rt = median_rt,
                                rtg = rt_group
                        )
                        # Calculate correlations if intensity data exists
                        if (!is.null(data)) {
                                cormat <- stats::cor(t(cluster_data[, -(1:2)]))
                                df$cor <- cormat[lower.tri(cormat)]
                        }
                        # Identify isotopes
                        iso_mask <- find_isotopes(df, digits = digits,didx=didx, corcutoff=corcutoff)
                        edges <- df[iso_mask, c("ms1", "ms2")]
                        iso_mz <- if (nrow(edges) > 0) {
                                g <- igraph::graph_from_data_frame(edges, directed = FALSE)
                                clusters <- igraph::components(g)
                                vertex_names <- names(clusters$membership)
                                keep_masses <- tapply(vertex_names, clusters$membership, function(x) min(as.numeric(x)))
                                remove_vertices <- setdiff(vertex_names, as.character(keep_masses))
                                as.numeric(remove_vertices)
                        } else {
                                numeric(0)
                        }
                        df_iso <- df[!(df$ms1%in%iso_mz|df$ms2%in%iso_mz), ]
                        results$iso <- df[iso_mask, ]
                        # Identify multi-charged isotope ions
                        multiiso_mask <- find_multiiso(df_iso, didx=didx, corcutoff=corcutoff)
                        results$multiiso <- df_iso[multiiso_mask, ]
                        multiiso_mz <- unique(c(results$multiiso$ms1,results$multiiso$ms2))
                        df_iso_multiiso <- df_iso[!(df_iso$ms1%in%multiiso_mz|df_iso$ms2%in%multiiso_mz), ]
                        # Identify multi-charged ions
                        multi_mask <- find_multicharged(df_iso_multiiso, didx=didx, corcutoff=corcutoff)
                        results$multi <- df_iso_multiiso[multi_mask, ]
                        df_iso_multiiso_multi <- df_iso_multiiso[!multi_mask,]
                        results$multimz <- unique(c(results$multi$ms1,results$multi$ms2))

                        results$diff <- df_iso_multiiso_multi
                        results$iso_mz <- iso_mz
                        results$multiiso_mz <- multiiso_mz

                        results
                }
        }

        # Helper functions
        find_isotopes <- function(df, digits, didx, corcutoff) {
                non_zero <- round(df$diff, digits) != 0

                cond_m1_low <- (df$diff %% 1 < 0.01) & (df$diff >= 1) & (df$diff < 2)
                cond_m2_core <- (df$diff %% 2 < 0.01) & (df$diff >= 2) & (df$diff < 3)
                cond_m1_high <- (df$diff %% 1 > 0.99) & (df$diff >= 1) & (df$diff < 2)
                cond_m0_round <- (df$diff %% 1 > 0.99) & (df$diff >= 0) & (df$diff < 1)

                iso_mask <- non_zero & (
                        cond_m1_low | cond_m2_core | cond_m1_high | cond_m0_round
                )
                if(didx){
                        iso_mask&df$cor>corcutoff
                }else{
                        iso_mask
                }
        }

        find_multiiso <- function(df, digits = 1, didx, corcutoff) {
                cond_z2_isotope <- (round(df$diff %% 1, 1) == 0.5) & (df$diff < 1)
                cond_z3_isotope <- (round(df$diff %% 1, 1) == 0.3) & (df$diff < 1)
                if(didx){
                        (cond_z2_isotope|cond_z3_isotope)&df$cor>corcutoff
                }else{
                        cond_z2_isotope|cond_z3_isotope
                }
        }

        find_multicharged <- function(df, digits = 1, didx, corcutoff) {
                multicharge <- round(df$diff %% 1, digits) == 0.5
                if(didx){
                        multicharge&df$cor>corcutoff
                }else{
                        multicharge
                }
        }

        # Combine results safely
        safe_extract <- function(element) {
                do.call(rbind, lapply(cluster_results, function(x) x[[element]]))
        }

        # Main processing
        mz <- round(list$mz,accuracy)
        rt <- list$rt
        data <- list$data
        rt_cluster <- rt_clusters(list$rt, rtcutoff)
        n_clusters <- length(unique(rt_cluster))
        message(n_clusters, " retention time clusters found.")

        ng <- if (is.null(ng)) round(n_clusters * 0.2) else ng
        message("Using ng = ", ng)

        # Prepare data as dataframe
        processed_data <- data.frame(mz = mz, rt = rt)

        if (!is.null(data)) {
                processed_data <- cbind(processed_data, as.data.frame(data))
        }
        # enable corcutoff filtering
        didx <- ifelse(!is.null(data)&!is.null(corcutoff),T,F)

        # Split data and process clusters
        split_list <- split(processed_data, rt_cluster)
        cluster_results <- Map(
                function(x, i, didx, corcutoff) {
                        if (nrow(x) > 0) analyze_pmd(x, i, didx, corcutoff) else list()
                },
                split_list,
                as.numeric(names(split_list)),
                didx,
                corcutoff = ifelse(didx,corcutoff,0)
        )

        list$iso <- safe_extract("iso")
        list$multiiso = safe_extract("multiiso")
        list$multi = safe_extract("multi")
        list$solo = safe_extract("solo")
        diff = safe_extract("diff")

        # Post-processing
        list$rtcluster <- rt_cluster

        # Cor filtering
        if(didx){
                diff <- diff[diff$cor>corcutoff,]
        }

        # Frequency filtering
        diff2 <- round(diff$diff, digits)
        difflist <- split(diff2,diff$rtg)
        diff3 <- sapply(difflist,function(x) unique(x))
        diff3 <- do.call(c,diff3)
        pmd_freq <- table(as.numeric(diff3))
        common_pmd <- as.numeric(names(pmd_freq[pmd_freq > ng]))

        list$paired <- diff[
                round(diff$diff, digits) %in% round(common_pmd,digits),
        ]

        # Collect isotope and multicharged mz values
        iso_mz_all <- unlist(lapply(cluster_results, function(x) paste0(x$iso_mz,'@',x$rtg)))
        multiiso_mz_all <- unlist(lapply(cluster_results, function(x) paste0(x$multiiso_mz,'@',x$rtg)))
        multi_mz_all <- unlist(lapply(cluster_results, function(x) paste0(x$multimz,'@',x$rtg)))
        paired_mz_all <- unique(c(paste0(list$paired$ms1,'@',list$paired$rtg),paste0(list$paired$ms2,'@',list$paired$rtg)))
        multi_mz_all2 <- multi_mz_all[!(multi_mz_all%in%paired_mz_all)]
        solo_mz <- paste0(list$solo$mz,'@',list$solo$rtg)

        # Create logical vectors for isotope and multicharged
        mzrtg <- paste0(round(mz, accuracy),'@',list$rtcluster)
        list$soloindex <- mzrtg %in% solo_mz
        list$isoindex <- mzrtg %in% iso_mz_all
        list$multiisoindex <- mzrtg %in% multiiso_mz_all
        list$multiindex <- mzrtg %in% multi_mz_all2
        list$pairedindex <- mzrtg %in% paired_mz_all
        # Reporting
        message(length(unique(round(list$paired$diff,digits))), " unique PMDs retained.")
        message(paste(c("The unique within RT clusters high frequency PMD(s) is(are) ", unique(round(list$paired$diff,digits))), collapse=" "), '.')
        message(sum(list$isoindex), " isotope peaks found.")
        message(sum(list$multiisoindex), " multiple charged isotope peaks found.")
        message(sum(list$multiindex), " multiple charged peaks found.")
        message(sum(list$pairedindex), " paired peaks found.")

        return(list)
}

#' Identify standard ions through retention time clustering and PMD relationships
#'
#' @param list A list object from getpaired() containing paired features
#' @param digits Rounding digits for mass differences
#' @param accuracy Mass accuracy for standard ion identification
#' @return List with added standard ion indices and metadata
#' @examples
#' data(spmeinvivo)
#' pmd <- getpaired(spmeinvivo)
#' std <- getstd(pmd)
#' @seealso \code{\link{getpaired}},\code{\link{getsda}},\code{\link{plotstd}}
#' @export
getstd <- function(list, digits = 2, accuracy = 4) {

        initialize_results <- function() {
                list(
                        solo = NULL,
                        groupA = NULL,
                        groupB1 = NULL,
                        groupB2 = NULL,
                        groupB3 = NULL
                )
        }
        group_descriptions <- c(
                solo = "group(s) have single peaks",
                A = "group(s) with multiple peaks while no isotope/paired relationship",
                B1 = "group(s) with isotope without paired relationship",
                B2 = "group(s) with paired without isotope relationship",
                B3 = "group(s) with both paired and isotope relationship"
        )

        process_group <- function(rt_group, type = c("solo", "A", "B1", "B2", "B3")) {
                type <- match.arg(type)
                switch(type,
                       "solo" = process_solo(rt_group),
                       "A" = process_groupA(rt_group),
                       "B1" = process_groupB1(rt_group),
                       "B2" = process_groupB2(rt_group),
                       "B3" = process_groupB3(rt_group)
                )
        }

        # group group 1: RT groups with solo peak
        process_solo <- function(rtg) {
                list$solo[list$solo$rtg == rtg, c("mz", "rt", "rtg")]
        }
        # group 2: RT groups with multiple peaks group 2A: RT groups with multiple peaks while no isotope/paired relationship
        process_groupA <- function(rtg) {
                mz <- max(list$mz[list$rtcluster == rtg])
                rt <- stats::median(list$rt[list$rtcluster == rtg])
                cbind(mz, rt, rtg)
        }
        # group 2B: RT groups with multiple peaks with isotope/paired relationship
        # group 2B1: RT groups with multiple peaks with isotope without paired relationship
        process_groupB1 <- function(rtg) {
                df_iso <- resultiso[resultiso$rtg == rtg, ]
                if (nrow(df_iso) == 0) return(NULL)

                iso_mz <- if (nrow(df_iso) > 0) {
                        g <- igraph::graph_from_data_frame(df_iso, directed = FALSE)
                        clusters <- igraph::components(g)
                        vertex_names <- names(clusters$membership)
                        keep_masses <- tapply(vertex_names, clusters$membership, function(x) min(as.numeric(x)))
                        as.numeric(keep_masses)
                } else {
                        numeric(0)
                }

                cbind(mz = unique(iso_mz), rt = df_iso$rt[1], rtg = df_iso$rtg[1])
        }
        # group 2B2: RT groups with multiple peaks with paired relationship without isotope
        process_groupB2 <- function(rtg) {
                df_paired <- resultdiff[resultdiff$rtg == rtg, ]
                if (nrow(df_paired) == 0) return(NULL)

                mass_std <- if (nrow(df_paired) > 0) {
                        g <- igraph::graph_from_data_frame(df_paired, directed = FALSE)
                        clusters <- igraph::components(g)
                        vertex_names <- names(clusters$membership)
                        keep_masses <- tapply(vertex_names, clusters$membership, function(x) min(as.numeric(x)))
                        as.numeric(keep_masses)
                } else {
                        numeric(0)
                }

                cbind(mz = unique(mass_std), rt = df_paired$rt[1], rtg = df_paired$rtg[1])
        }
        # group 2B3: RT groups with multiple peaks with paired relationship and isotope
        process_groupB3 <- function(rtg) {
                df_iso <- resultiso[resultiso$rtg == rtg, ]
                df_paired <- resultdiff[resultdiff$rtg == rtg, ]
                df <- rbind.data.frame(df_iso,df_paired)

                if (nrow(df) == 0) return(NULL)

                mass_std <- if (nrow(df) > 0) {
                        g <- igraph::graph_from_data_frame(df, directed = FALSE)
                        clusters <- igraph::components(g)
                        vertex_names <- names(clusters$membership)
                        keep_masses <- tapply(vertex_names, clusters$membership, function(x) min(as.numeric(x)))
                        as.numeric(keep_masses)
                } else {
                        numeric(0)
                }

                cbind(mz = unique(mass_std),
                      rt = ifelse(nrow(df_iso) > 0, df_iso$rt[1], df_paired$rt[1]),
                      rtg = ifelse(nrow(df_iso) > 0, df_iso$rtg[1], df_paired$rtg[1]))
        }
        # main function
        resultdiff <- list$paired
        resultiso <- list$iso

        # process by group
        rt_groups <- unique(list$rtcluster)
        group_types <- list(
                solo = list$solo$rtg,
                A = rt_groups[!rt_groups %in% c(unique(resultdiff$rtg), unique(resultiso$rtg), list$solo$rtg)],
                B1 = rt_groups[rt_groups %in% setdiff(unique(resultiso$rtg), c(unique(resultdiff$rtg), list$solo$rtg))],
                B2 = rt_groups[rt_groups %in% setdiff(unique(resultdiff$rtg), c(unique(resultiso$rtg), list$solo$rtg))],
                B3 = rt_groups[rt_groups %in% intersect(unique(resultdiff$rtg), unique(resultiso$rtg))]
        )

        results <- lapply(names(group_types), function(type) {
                groups <- group_types[[type]]
                if (length(groups) == 0) return(NULL)

                full_rtg <- switch(type,
                                   "solo" = list$rtcluster[list$soloindex],
                                   "A" = rt_groups[!rt_groups %in% c(unique(resultdiff$rtg), unique(resultiso$rtg), list$solo$rtg)],
                                   "B1" = rt_groups[rt_groups %in% setdiff(unique(resultiso$rtg), c(unique(resultdiff$rtg), list$solo$rtg))],
                                   "B2" = rt_groups[rt_groups %in% setdiff(unique(resultdiff$rtg), c(unique(resultiso$rtg), list$solo$rtg))],
                                   "B3" = rt_groups[rt_groups %in% intersect(unique(resultdiff$rtg), unique(resultiso$rtg))]
                )

                message(paste(
                        c(length(groups),
                          group_descriptions[type],
                          if(length(full_rtg) <= 10) full_rtg else c(utils::head(full_rtg,5),"...",utils::tail(full_rtg,5))),
                        collapse = " "
                ))

                do.call(rbind, lapply(groups, function(g) process_group(g, type)))
        })

        resultstd <- do.call(rbind, c(list(list$solo[, c("mz", "rt", "rtg")]), results))
        resultstd <- unique(stats::na.omit(resultstd))
        colnames(resultstd) <- c("mz", "rt", "rtg")
        list$stdmassindex <- paste(round(list$mz, accuracy), list$rtcluster) %in%
                paste(round(resultstd$mz, accuracy), resultstd$rtg)

        list$stdmass <- resultstd
        message(sprintf("%d standard masses identified.", sum(list$stdmassindex)))

        return(list)
}

#' GlobalStd algorithm with structure/reaction directed analysis
#' @param list a peaks list with mass to charge, retention time and intensity data
#' @param rtcutoff cutoff of the distances in cluster, default 10
#' @param ng cutoff of global PMD's retention time group numbers, If ng = NULL, 20 percent of RT cluster will be used as ng, default NULL.
#' @param corcutoff cutoff of the correlation coefficient, default NULL
#' @param digits mass or mass to charge ratio accuracy for pmd, default 2
#' @param accuracy measured mass or mass to charge ratio in digits, default 4
#' @param freqcutoff pmd frequency cutoff for structures or reactions, default NULL. This cutoff will be found by PMD network analysis when it is NULL.
#' @param sda logical, option to perform structure/reaction directed analysis, default FALSE.
#' @return list with GlobalStd algorithm processed data.
#' @examples
#' data(spmeinvivo)
#' re <- globalstd(spmeinvivo)
#' @seealso \code{\link{getpaired}},\code{\link{getstd}},\code{\link{getsda}},\code{\link{plotstd}},\code{\link{plotstdsda}},\code{\link{plotstdrt}}
#' @export
globalstd <- function(list,
                      rtcutoff = 10,
                      ng = NULL,
                      corcutoff = NULL,
                      digits = 2,
                      accuracy = 4,
                      freqcutoff = NULL,
                      sda = FALSE) {
        list <-
                getpaired(
                        list,
                        rtcutoff = rtcutoff,
                        corcutoff = corcutoff,
                        ng = ng,
                        digits = digits,
                        accuracy = accuracy
                )
        if (sum(list$pairedindex) > 0) {
                list2 <-
                        getstd(
                                list,
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
#' @param list A list containing peak intensities, m/z values, and retention times. Must include elements: mz, rt, data.
#' @param corcutoff Cutoff value for correlation coefficient (default: 0.9).
#' @param rtcutoff Cutoff value for retention time clustering (default: 10).
#' @param accuracy Number of decimal places for m/z rounding (default: 4).
#' @param ng cutoff of global PMD's retention time group numbers, If ng = NULL, 20 percent of RT cluster will be used as ng, default NULL.
#' @param digits mass or mass to charge ratio accuracy for pmd, default 2
#' @return A list with pseudo-spectrum clustering results.
#' @examples
#' data(spmeinvivo)
#' pseudo <- getpseudospectrum(spmeinvivo)
#' @export
getpseudospectrum <- function(list,
                          corcutoff = NULL,
                          rtcutoff = 10,
                          accuracy = 4,
                          ng = NULL,
                          digits = 2) {

        # Input validation
        if (!all(c("mz", "rt", "data") %in% names(list))) {
                stop("Input list must contain 'mz', 'rt', and 'data'.")
        }
        if (length(list$mz) != nrow(list$data) || length(list$rt) != nrow(list$data)) {
                stop("Length of 'mz' and 'rt' must match the number of rows in 'data'.")
        }

        # group 2: RT groups with multiple peaks group 2A: RT groups with multiple peaks while no isotope/paired relationship
        process_groupA <- function(rtg,didx) {
                if(didx){
                        # Pairwise ions
                        mz_pairs <- t(utils::combn(mz[list$rtcluster == rtg], 2))
                        # Create base dataframe
                        df <- data.frame(
                                ms1 = pmin(mz_pairs[, 1], mz_pairs[, 2]),
                                ms2 = pmax(mz_pairs[, 1], mz_pairs[, 2]),
                                rtg = rtg
                        )
                        data <- list$data[list$rtcluster == rtg,]
                        cormat <- stats::cor(t(data))
                        df$cor <- cormat[lower.tri(cormat)]
                        dfcor <- df[df$cor>corcutoff,]
                        if(nrow(dfcor)>0){
                                edges <- dfcor[,c("ms1", "ms2")]
                                g <- igraph::graph_from_data_frame(edges, directed = FALSE)
                                clusters <- igraph::components(g)
                                mem <- clusters$membership
                                mzs <- as.numeric(names(clusters$membership))
                                mzo <- setdiff(mz[list$rtcluster == rtg],mzs)
                                if(length(mzo)>0){
                                        #pseudoA <- cbind.data.frame(mz=c(mzs,mzo),ins=msdata[match(c(mzs,mzo),mz[list$rtcluster == rtg])],rtg=rtg,sid=paste0(rtg,'@',c(mem,(length(mem)+1):(length(mzo)+length(mem)))))
                                        pseudoA <- cbind.data.frame(mz=c(mzs),ins=msdata[match(mzs,mz[list$rtcluster == rtg])],rtg=rtg,sid=paste0(rtg,'@',mem))
                                }else{
                                        pseudoA <- cbind.data.frame(mz=mz[list$rtcluster == rtg],ins=msdata[list$rtcluster == rtg],rtg=rtg,sid=paste0(rtg,'@',mem))
                                }

                        }else{
                                #pseudoA <- cbind.data.frame(mz=mz[list$rtcluster == rtg],ins=msdata[list$rtcluster == rtg],rtg=rtg,sid=paste0(rtg,'@',seq_along(mz[list$rtcluster == rtg])))
                                pseudoA <- cbind.data.frame(mz=mz[list$rtcluster == rtg],ins=msdata[list$rtcluster == rtg],rtg=rtg,sid=paste0(rtg,'@1'))
                        }

                }else{
                        #pseudoA <- cbind.data.frame(mz=mz[list$rtcluster == rtg],ins=msdata[list$rtcluster == rtg],rtg=rtg,sid=paste0(rtg,'@',seq_along(mz[list$rtcluster == rtg])))
                        pseudoA <- cbind.data.frame(mz=mz[list$rtcluster == rtg],ins=msdata[list$rtcluster == rtg],rtg=rtg,sid=paste0(rtg,'@1'))
                }
                pseudoA
        }
        # group 2B: RT groups with multiple peaks with isotope/paired relationship
        process_groupB <- function(rtg,didx) {
                df <- resultdf[resultdf$rtg == rtg, ]
                if(didx){
                        dfcor <- df[df$cor>corcutoff,]
                        if(nrow(dfcor)>0){
                                edges <- dfcor[,c("ms1", "ms2")]
                                g <- igraph::graph_from_data_frame(edges, directed = FALSE)
                                clusters <- igraph::components(g)
                                mem <- clusters$membership
                                mzs <- as.numeric(names(clusters$membership))
                                mzo <- setdiff(mz[list$rtcluster == rtg],mzs)
                                if(length(mzo)>0){
                                        #pseudoB <- cbind.data.frame(mz=c(mzs,mzo),ins=msdata[match(c(mzs,mzo),mz[list$rtcluster == rtg])],rtg=rtg,sid=paste0(rtg,'@',c(mem,(length(mem)+1):(length(mzo)+length(mem)))))
                                        pseudoB <- cbind.data.frame(mz=mzs,ins=msdata[match(mzs,mz[list$rtcluster == rtg])],rtg=rtg,sid=paste0(rtg,'@',mem))
                                }else{
                                        pseudoB <- cbind.data.frame(mz=mzs,ins=msdata[match(mzs,mz[list$rtcluster == rtg])],rtg=rtg,sid=paste0(rtg,'@',mem))
                                }
                        }else{
                                #pseudoB <- cbind.data.frame(mz=mz[list$rtcluster == rtg],ins=msdata[list$rtcluster == rtg],rtg=rtg,sid=paste0(rtg,'@',seq_along(mz[list$rtcluster == rtg])))
                                pseudoB <- cbind.data.frame(mz=mz[list$rtcluster == rtg],ins=msdata[list$rtcluster == rtg],rtg=rtg,sid=paste0(rtg,'@1'))
                        }
                }else{
                        edges <- df[,c("ms1", "ms2")]
                        g <- igraph::graph_from_data_frame(edges, directed = FALSE)
                        clusters <- igraph::components(g)
                        mem <- clusters$membership
                        mzs <- as.numeric(names(clusters$membership))
                        pseudoB <- cbind.data.frame(mz=mzs,ins=msdata[match(mzs,mz[list$rtcluster == rtg])],rtg=rtg,sid=paste0(rtg,'@',mem))
                }
                pseudoB
        }
        # main function
        # Generate peaks group
        list <-
                getpaired(
                        list,
                        rtcutoff = rtcutoff,
                        corcutoff = corcutoff,
                        ng = ng,
                        digits = digits,
                        accuracy = accuracy
                )
        resultdiff <- list$paired
        resultiso <- list$iso
        resultmultiiso <- list$multiiso
        resultmulti <- list$multi
        resultdf <- rbind.data.frame(resultdiff,resultiso,resultmultiiso,resultmulti)
        mz <- round(list$mz,accuracy)
        data <- list$data
        # enable corcutoff filtering
        didx <- ifelse(!is.null(data)&!is.null(corcutoff),T,F)
        # Generate spectrum
        if (is.null(list$data)) {
                message('You need intensity data to use corcutoff and export pseudospectra')
                msdata <- NULL
        } else{
                data <- list$data
                if (is.matrix(data)|is.data.frame(data)) {
                        msdata <- apply(data, 1, sum)
                } else {
                        msdata <-sum(data)
                }
        }
        pseudosolo <- cbind.data.frame(mz=mz[list$soloindex],ins=msdata[list$soloindex],rtg=list$rtcluster[list$soloindex],sid=paste0(list$rtcluster[list$soloindex],'@1'))

        # process by group
        rt_groups <- unique(list$rtcluster[!list$soloindex])
        A = rt_groups[!rt_groups %in% c(unique(resultdiff$rtg), unique(resultiso$rtg), list$solo$rtg)]
        B = rt_groups[rt_groups %in% c(unique(resultdiff$rtg), unique(resultiso$rtg))]
        resultA <- list()
        for(i in 1:length(A)){
                resultA[[i]] <- process_groupA(A[i],didx)
        }
        resultB <- list()
        for(i in 1:length(B)){
                resultB[[i]] <- process_groupB(B[i],didx)
        }
        resultstdA <- do.call(rbind, resultA)
        resultstdB <- do.call(rbind, resultB)

        list$pseudo <- rbind(pseudosolo,resultstdA,resultstdB)
        pseudolist <- split(list$pseudo,list$pseudo$sid)
        mzlist <- sapply(pseudolist, function(x) paste(x$mz[which.max(x$ins)],x$rtg))
        mzh <- do.call(c,mzlist)

        list$stdmassindex2 <- paste(round(list$mz, accuracy), list$rtcluster) %in% mzh
        coverage <- paste(round(list$mz, accuracy), list$rtcluster) %in%
                paste(round(list$pseudo$mz, accuracy), list$pseudo$rtg)
        message(paste(length(unique(list$pseudo$sid))),' pseudo spectrum found.')
        message(paste(round(sum(coverage)/length(coverage),2),'peaks covered by PMD relationship.'))
        return(list)
}

#' Get Pseudo-Spectrum as peaks cluster based on correlation analysis.
#' @param list a list with peaks intensity
#' @param corcutoff cutoff of the correlation coefficient, default 0.9
#' @param rtcutoff cutoff of the distances in cluster, default 10
#' @param accuracy measured mass or mass to charge ratio in digits, default 4
#' @return list with Pseudo-Spectrum index
#' @examples
#' data(spmeinvivo)
#' pseudo <- getcorpseudospectrum(spmeinvivo)
#' @export
getcorpseudospectrum <- function(list,
                          corcutoff = 0.9,
                          rtcutoff = 10,
                          accuracy = 4) {
        process_group <- function(subset_df) {
                max_b_index <- which.max(subset_df$ins)
                combined_value <- paste0(subset_df$mz[max_b_index], "@", subset_df$rtg[max_b_index])
                return(combined_value)
        }
        process_group2 <- function(subset_df) {
                max_b_index <- which.max(subset_df$mz)
                combined_value <- paste0(subset_df$mz[max_b_index], "@", subset_df$rtg[max_b_index])
                return(combined_value)
        }
        mz <- round(list$mz,accuracy)
        rt <- list$rt

        dis <- stats::dist(list$rt, method = "manhattan")
        fit <- stats::hclust(dis)
        list$rtcluster <- stats::cutree(fit, h = rtcutoff)
        n <- length(unique(list$rtcluster))
        message(paste(n, "retention time cluster found."))
        data <- list$data
        pseudo <- mzsall <- mzall <-list()
        for (i in seq_along(unique(list$rtcluster))) {
                # find the mass within RT
                if(sum(list$rtcluster==i)>1){
                        data <- list$data[list$rtcluster == i,]
                        if (is.matrix(data)) {
                                msdata <- apply(data, 1, sum)
                        } else {
                                msdata <- sum(data)
                        }
                        # Pairwise ions
                        mz_pairs <- t(utils::combn(mz[list$rtcluster == i], 2))
                        # Create base dataframe
                        df <- data.frame(
                                ms1 = pmin(mz_pairs[, 1], mz_pairs[, 2]),
                                ms2 = pmax(mz_pairs[, 1], mz_pairs[, 2]),
                                rtg = i
                        )

                        cormat <- stats::cor(t(data))
                        df$cor <- cormat[lower.tri(cormat)]
                        dfcor <- df[df$cor>corcutoff,]
                        if(nrow(dfcor)>0){
                                edges <- dfcor[,c("ms1", "ms2")]
                                g <- igraph::graph_from_data_frame(edges, directed = FALSE)
                                clusters <- igraph::components(g)
                                mem <- clusters$membership
                                mzs <- round(as.numeric(names(clusters$membership)),accuracy)
                                mzo <- setdiff(mz[list$rtcluster == i],mzs)
                                if(length(mzo)>0){
                                        #pseudo[[i]] <- dt <-cbind.data.frame(mz=c(mzs,mzo),ins=msdata[match(c(mzs,mzo),mz[list$rtcluster == i])],rtg=i,sid=paste0(i,'@',c(mem,(length(mem)+1):(length(mzo)+length(mem)))))
                                        pseudo[[i]] <- dt <-cbind.data.frame(mz=mzs,ins=msdata[match(mzs,mz[list$rtcluster == i])],rtg=i,sid=paste0(i,'@',mem))
                                        split_df <- split(dt, dt$sid)
                                        mzsall[[i]] <- sapply(split_df, process_group)
                                        mzall[[i]] <- sapply(split_df, process_group2)
                                }else{
                                        pseudo[[i]] <- dt <-cbind.data.frame(mz=mz[list$rtcluster == i],ins=msdata,rtg=i,sid=paste0(i,'@',mem))
                                        split_df <- split(dt, dt$sid)
                                        mzsall[[i]] <- sapply(split_df, process_group)
                                        mzall[[i]] <- sapply(split_df, process_group2)
                                }
                        }else{
                                #pseudo[[i]] <- dt <- cbind.data.frame(mz=mz[list$rtcluster == i],ins=msdata,rtg=i,sid=paste0(i,'@',seq_along(mz[list$rtcluster == i])))
                                pseudo[[i]] <- dt <- cbind.data.frame(mz=mz[list$rtcluster == i],ins=msdata,rtg=i,sid=paste0(i,'@1'))
                                split_df <- split(dt, dt$sid)
                                mzsall[[i]] <- sapply(split_df, process_group)
                                mzall[[i]] <- sapply(split_df, process_group2)
                        }
                }else{
                        data <- list$data[list$rtcluster == i,]
                        if (is.matrix(data)) {
                                msdata <- apply(data, 1, sum)
                        } else {
                                msdata <- sum(data)
                        }
                        mzt <- mz[list$rtcluster == i]
                        mzall[[i]] <- mzsall[[i]] <- paste0(mzt,'@',i)
                        pseudo[[i]] <- cbind.data.frame(mz=mzt, ins=msdata, rtg=i, sid=paste0(i,'@1'))
                }
        }
        list$pseudo <- do.call(rbind, pseudo)
        mzsall <- do.call(c, mzsall)
        mzall <- do.call(c, mzall)
        list$stdmassindex <- paste0(mz, '@', list$rtcluster) %in% mzall
        list$stdmassindex2 <-
                paste0(mz, '@', list$rtcluster) %in% mzsall
        coverage <- paste(round(list$mz, accuracy), list$rtcluster) %in%
                paste(round(list$pseudo$mz, accuracy), list$pseudo$rtg)
        message(paste(length(unique(list$pseudo$sid))),' pseudo spectrum found.')
        message(paste(round(sum(coverage)/length(coverage),2),'peaks covered by PMD relationship.'))
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
                dim <- dim(cov.x)[1]
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
