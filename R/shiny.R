#' @export
runPMD <- function() {
    file <- system.file("shinyapp", "PMD", "PMD.Rmd",
        package = "pmd")
    if (file == "") {
        stop("Could not find directory. Try re-installing `pmd`.",
            call. = FALSE)
    }
    rmarkdown::run(file)
}
