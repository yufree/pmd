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
  cvcutoff = 30,
  outlier = FALSE,
  method = "static",
  ...
)
}
\arguments{
\item{list}{a list with mzrt profile and data}

\item{pmd}{a specific paired mass distances}

\item{rtcutoff}{cutoff of the distances in retention time hierarchical clustering analysis, default 10}

\item{digits}{mass or mass to charge ratio accuracy for pmd, default 2}

\item{accuracy}{measured mass or mass to charge ratio in digits, default 4}

\item{cvcutoff}{ratio or intensity cv cutoff for quantitative paired peaks, default 30}

\item{outlier}{logical, if true, outlier of ratio will be removed, default False.}

\item{method}{quantification method can be 'static' or 'dynamic'. See details.}

\item{...}{other parameters for getpmd}
}
\value{
list with quantitative paired peaks.
}
\description{
Get quantitative paired peaks list for specific reaction/pmd
}
\details{
PMD based reaction quantification methods have two options: 'static' will only consider the stable mass pairs across samples and such reactions will be limited by the enzyme or other factors than substrates. 'dynamic' will consider the unstable paired masses by normalization the relatively unstable peak with stable peak between paired masses and such reactions will be limited by one or both peaks in the paired masses.
}
\examples{
data(spmeinvivo)
pmd <- getreact(spmeinvivo,pmd=15.99)
}
\seealso{
\code{\link{getpaired}},\code{\link{getstd}},\code{\link{getsda}},\code{\link{getrda}},\code{\link{getpmd}},
}
