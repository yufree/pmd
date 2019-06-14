#' A dataset containing common Paired mass distances of substructure, ions replacements, and reaction
#' @docType data
#' @usage data(sda)
#' @format A data frame with 94 rows and 4 variables:
#' \describe{
#'   \item{PMD}{Paired mass distances}
#'   \item{origin}{potentical sources}
#'   \item{Ref.}{references}
#'   \item{mode}{positive, negative or both mode to find corresponding PMDs}
#'   }
"sda"

#' A peaks list dataset containing 9 samples from 3 fish with triplicates samples for each fish from LC-MS.
#' @docType data
#' @usage data(spmeinvivo)
#' @format A list with 4 variables from 1459 LC-MS peaks:
#' \describe{
#'   \item{mz}{mass to charge ratios}
#'   \item{rt}{retention time}
#'   \item{data}{intensity matrix}
#'   \item{group}{group information}
#'   }
"spmeinvivo"

#' A list dataset containing HMDB unique accurate mass pmd analysis results
#' @docType data
#' @usage data(hmdbp)
#' @format A list with two vectors
#' \describe{
#'   \item{massp}{all unique hmdb mass probability across all pmds}
#'   \item{pmdp}{pmds probability across all unique hmdb mass}
#'   }
"hmdbp"

#' A dataframe containing HMDB top 10000 unique accurate mass pmd and related reactions
#' @docType data
#' @usage data(hmdb)
#' @format A dataframe with atoms numbers of C, H, O, N, P, S
#' \describe{
#'   \item{percentage}{accuracy of atom numbers prediction}
#'   \item{pmd}{pmd with two digits}
#'   }
"hmdb"
