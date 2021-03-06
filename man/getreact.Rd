% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sda.R
\name{getreact}
\alias{getreact}
\title{Get quantitative paired peaks list for specific reaction/pmd}
\usage{
getreact(
  list,
  pmd,
  rtcutoff = 10,
  digits = 2,
  accuracy = 4,
  ratiocv = 30,
  outlier = FALSE,
  ...
)
}
\arguments{
\item{list}{a list with mzrt profile and data}

\item{pmd}{a specific paired mass distances}

\item{rtcutoff}{cutoff of the distances in retention time hierarchical clustering analysis, default 10}

\item{digits}{mass or mass to charge ratio accuracy for pmd, default 2}

\item{accuracy}{measured mass or mass to charge ratio in digits, default 4}

\item{ratiocv}{ratio cv cutoff for quantitative paired peaks, default 30}

\item{outlier}{logical, if true, outlier of ratio will be removed, default False.}

\item{...}{other parameters for getpmd}
}
\value{
list with quantitative paired peaks.
}
\description{
Get quantitative paired peaks list for specific reaction/pmd
}
\examples{
data(spmeinvivo)
pmd <- getreact(spmeinvivo,pmd=15.99)
}
\seealso{
\code{\link{getpaired}},\code{\link{getstd}},\code{\link{getsda}},\code{\link{getrda}},\code{\link{getpmd}},
}
