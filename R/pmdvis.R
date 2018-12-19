#' Plot the retention time group
#' @param list a list from getpaired function
#' @param ... other parameters for plot function
#' @return NULL
#' @seealso \code{\link{getpaired}}, \code{\link{globalstd}}
#' @examples
#' \donttest{
#' data(spmeinvivo)
#' pmd <- getpaired(spmeinvivo)
#' plotrtg(pmd)
#' }
#' @export
plotrtg <- function(list, ...) {
        std <- list$data
        col <-
                (grDevices::colorRampPalette(rev(
                        RColorBrewer::brewer.pal(11, "RdYlBu")
                )))(length(unique(list$rtcluster)))
        graphics::par(mfrow = c(2, 1), mar = c(4, 4, 2, 1) + 0.1)
        graphics::plot(
                list$rt,
                list$mz,
                xlab = 'retention time(s)',
                ylab = 'm/z',
                pch = 19,
                col = col[list$rtcluster],
                ...
        )
        graphics::barplot(
                table(list$rtcluster),
                col = col[unique(list$rtcluster)],
                ylab = 'peak numbers',
                las = 2,
                xlab = 'retention time group',
                cex.names = 0.5
        )
}
#' Plot the mass pairs and high frequency mass distances
#' @param list a list from getpaired function
#' @param index index for PMD value
#' @param ... other parameters for plot function
#' @return NULL
#' @seealso \code{\link{getpaired}}, \code{\link{globalstd}}
#' @examples
#' \donttest{
#' data(spmeinvivo)
#' pmd <- getpaired(spmeinvivo)
#' plotpaired(pmd)
#' }
#' @export
plotpaired <- function(list, index = NULL, ...) {
        paired <- list$paired
        if (is.null(index)) {
                diffgroup <- as.numeric(as.factor(paired$diff2))
                col <-
                        (grDevices::colorRampPalette(rev(
                                RColorBrewer::brewer.pal(11,
                                                         "RdYlBu")
                        )))(length(unique(paired$diff2)))
                graphics::par(mfrow = c(2, 1),
                              mar = c(4, 4, 2,
                                      1) + 0.1)
                graphics::plot(
                        range(paired$rt),
                        range(paired$ms1,
                              paired$ms2),
                        type = "n",
                        xlab = "retention time(s)",
                        ylab = "m/z",
                        ...
                )
                graphics::segments(
                        paired$rt,
                        paired$ms1,
                        paired$rt,
                        paired$ms2,
                        col = col[diffgroup],
                        lwd = 1.5
                )
                graphics::barplot(
                        table(list$paired$diff2),
                        col = col,
                        ylab = "Frequency",
                        las = 2,
                        xlab = "Paired mass distance",
                        cex.names = 0.618
                )
        } else {
                paired <- paired[index,]
                graphics::plot(
                        range(list$paired$rt),
                        range(list$paired$ms1,
                              list$paired$ms2),
                        type = "n",
                        xlab = "retention time(s)",
                        ylab = "m/z",
                        main = paste(paired$diff2[1],
                                     "group"),
                        ...
                )
                graphics::segments(
                        paired$rt,
                        paired$ms1,
                        paired$rt,
                        paired$ms2,
                        lwd = 1.5,
                        pch = 19
                )
                graphics::points(
                        paired$rt,
                        paired$ms1,
                        pch = 19,
                        col = grDevices::rgb(0, 0, 1, alpha = 0.318)
                )
                graphics::points(
                        paired$rt,
                        paired$ms2,
                        pch = 19,
                        col = grDevices::rgb(0, 0, 1, alpha = 0.318)
                )
        }

}

#' Plot the std mass from GlobalStd algorithm
#' @param list a list from getstd function
#' @return NULL
#' @seealso \code{\link{getstd}}, \code{\link{globalstd}}
#' @examples
#' \donttest{
#' data(spmeinvivo)
#' pmd <- getpaired(spmeinvivo)
#' std <- getstd(pmd)
#' plotstd(std)
#' }
#' @export
plotstd <- function(list) {
        std <- list$stdmass
        graphics::par(mfrow = c(1, 2), mar = c(4, 4, 2, 1) +
                              0.1)
        col <- grDevices::rgb(0, 0, 1, alpha = 0.318)
        graphics::plot(
                list$rt,
                list$mz,
                xlab = "retention time(s)",
                ylab = "m/z",
                pch = 19,
                col = col,
                main = "all peaks"
        )
        graphics::plot(
                std$rt,
                std$mz,
                xlab = "retention time(s)",
                ylab = "m/z",
                pch = 19,
                col = col,
                main = "GlobalStd peaks"
        )
}

#' Plot the std mass from GlobalStd algorithm in certain retention time groups
#' @param list a list from getstd function
#' @param rtcluster retention time group index
#' @param ... other parameters for plot function
#' @return NULL
#' @seealso \code{\link{getstd}}, \code{\link{globalstd}},\code{\link{plotstd}},\code{\link{plotpaired}},\code{\link{plotstdsda}}
#' @examples
#' \donttest{
#' data(spmeinvivo)
#' pmd <- getpaired(spmeinvivo)
#' std <- getstd(pmd)
#' plotstdrt(std,rtcluster = 6)
#' }
#' @export
#'
plotstdrt <- function(list, rtcluster, ...) {
        data <- list$data[list$rtcluster == rtcluster,]
        if (length(data) > ncol(list$data)) {
                msdata <- apply(data, 1, mean)
        } else {
                msdata <- mean(data)
        }
        mz <- list$mz[list$rtcluster == rtcluster]
        rt <- stats::median(list$rt[list$rtcluster == rtcluster])
        graphics::plot(
                mz,
                msdata,
                type = "h",
                xlab = paste("m/z",
                             "@", rt, "s"),
                ylab = "Intensity",
                ...
        )
        stdmz <- list$stdmass$mz[list$stdmass$rtg == rtcluster]
        index <- round(mz, 4) %in% round(stdmz, 4)

        graphics::points(mz[index],
                         msdata[index],
                         type = "h",
                         lwd = 2,
                         col = "red")
}

#' Plot the std mass from GlobalStd algorithm in structure directed analysis(SDA) groups
#' @param list a list from getsda function
#' @param index index for PMD value
#' @param ... other parameters for plot function
#' @return NULL
#' @seealso \code{\link{getstd}}, \code{\link{globalstd}},\code{\link{plotstd}},\code{\link{plotpaired}},\code{\link{plotstdrt}}
#' @examples
#' \donttest{
#' data(spmeinvivo)
#' re <- globalstd(spmeinvivo)
#' plotstdsda(re)
#' }
#' @export
plotstdsda <- function(list, index = NULL, ...) {
        sda <- list$sda
        diffgroup <- as.numeric(as.factor(sda$diff2))
        if (is.null(index)) {
                col <- (grDevices::colorRampPalette(rev(
                        RColorBrewer::brewer.pal(11,
                                                 "RdYlBu")
                )))(length(unique(diffgroup)))
                graphics::par(mfrow = c(2, 1),
                              mar = c(4, 4, 2,
                                      1) + 0.1)
                graphics::plot(
                        range(sda$rt1, sda$rt2),
                        range(sda$ms1,
                              sda$ms2),
                        type = "n",
                        xlab = "retention time(s)",
                        ylab = "m/z",
                        ...
                )
                graphics::segments(sda$rt1,
                                   sda$ms1,
                                   sda$rt2,
                                   sda$ms2,
                                   col = col[diffgroup],
                                   lwd = 1.5)
                graphics::points(
                        sda$rt1,
                        sda$ms1,
                        pch = 19,
                        col = grDevices::rgb(0,
                                             0, 1, alpha = 0.318)
                )
                graphics::points(
                        sda$rt2,
                        sda$ms2,
                        pch = 19,
                        col = grDevices::rgb(0,
                                             0, 1, alpha = 0.318)
                )
                graphics::barplot(
                        table(sda$diff2),
                        col = col,
                        ylab = "Frequency",
                        las = 2,
                        xlab = "Paired mass distances",
                        cex.names = 0.618
                )

        } else {
                sda <- list$sda[index,]
                graphics::plot(
                        range(list$sda$rt1, list$sda$rt2),
                        range(list$sda$ms1, list$sda$ms2),
                        type = "n",
                        xlab = "retention time(s)",
                        ylab = "m/z",
                        main = paste(sda$diff2[1],
                                     "group"),
                        ...
                )
                graphics::segments(sda$rt1,
                                   sda$ms1,
                                   sda$rt2,
                                   sda$ms2,
                                   lwd = 1.5,
                                   pch = 19)
                graphics::points(
                        sda$rt1,
                        sda$ms1,
                        pch = 19,
                        col = grDevices::rgb(0,
                                             0, 1, alpha = 0.318)
                )
                graphics::points(
                        sda$rt2,
                        sda$ms2,
                        pch = 19,
                        col = grDevices::rgb(0,
                                             0, 1, alpha = 0.318)
                )
        }

}
