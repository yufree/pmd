#' A dataset containing common Paired mass differeces of substructure, ions replacements, and reaction
#' @docType data
#' @usage data(sda)
#' @format A data frame with 94 rows and 4 variables:
#' \describe{
#'   \item{PMD}{Paired mass differeces}
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
