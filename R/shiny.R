#' Shiny application for PMD analysis
#' @export
runPMD <- function() {
    if (!requireNamespace("rmarkdown", quietly = TRUE)) {
        stop("Package 'rmarkdown' is required. Install with install.packages('rmarkdown').",
             call. = FALSE)
    }
    file <- system.file("shinyapp", "PMD.Rmd", package = "pmd")
    if (file == "") {
        stop("Could not find directory. Try re-installing `pmd`.",
            call. = FALSE)
    }
    rmarkdown::run(file)
}
#' Shiny application for PMD network analysis
#' @export
runPMDnet <- function() {
    if (!requireNamespace("rmarkdown", quietly = TRUE)) {
        stop("Package 'rmarkdown' is required. Install with install.packages('rmarkdown').",
             call. = FALSE)
    }
    file <- system.file("shinyapp", "pmdnet.Rmd", package = "pmd")
    if (file == "") {
        stop("Could not find directory. Try re-installing `pmd`.",
             call. = FALSE)
    }
    rmarkdown::run(file)
}
