% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sda.R
\name{getcda}
\alias{getcda}
\title{Perform correlation directed analysis for peaks list.}
\usage{
getcda(list, corcutoff = 0.9, rtcutoff = 10, accuracy = 4)
}
\arguments{
\item{list}{a list with mzrt profile}

\item{corcutoff}{cutoff of the correlation coefficient, default NULL}

\item{rtcutoff}{cutoff of the distances in retention time hierarchical clustering analysis, default 10}

\item{accuracy}{measured mass or mass to charge ratio in digits, default 4}
}
\value{
list with correlation directed analysis results
}
\description{
Perform correlation directed analysis for peaks list.
}
\examples{
data(spmeinvivo)
cluster <- getpseudospectrum(spmeinvivo)
cbp <- enviGCMS::getfilter(cluster,rowindex = cluster$stdmassindex2)
cda <- getcda(cbp)
}
\seealso{
\code{\link{getsda}},\code{\link{getrda}}
}
