% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pmdvis.R
\name{plotsda}
\alias{plotsda}
\title{Plot the specific structure directed analysis(SDA) groups}
\usage{
plotsda(list, ...)
}
\arguments{
\item{list}{a list from getpmd function}

\item{...}{other parameters for plot function}
}
\description{
Plot the specific structure directed analysis(SDA) groups
}
\examples{
data(spmeinvivo)
re <- getpmd(spmeinvivo,pmd=78.9)
plotsda(re)
}
\seealso{
\code{\link{getstd}}, \code{\link{globalstd}},\code{\link{plotstd}},\code{\link{plotpaired}},\code{\link{plotstdrt}}
}
