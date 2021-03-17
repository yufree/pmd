#' A dataset containing common Paired mass distances of substructure, ions replacements, and reaction
#' @docType data
#' @usage data(sda)
#' @format A data frame with 94 rows and 4 variables:
#' \describe{
#'   \item{PMD}{Paired mass distances}
#'   \item{origin}{potential sources}
#'   \item{Ref.}{references}
#'   \item{mode}{positive, negative or both mode to find corresponding PMDs}
#'   }
"sda"

#' mass spectrometry contaminants database for PMD check
#' @docType data
#' @usage data(sda)
#' @format A data frame with 308 rows and 5 variables:
#' \describe{
#'   \item{id}{MaConDa ID}
#'   \item{name}{contaminants}
#'   \item{formula}{contaminants fomula}
#'   \item{exact_mass}{exact mass of contaminants}
#'   \item{type_of_contaminant}{type of contaminant}
#'   }
#' @source \url{https://academic.oup.com/bioinformatics/article/28/21/2856/236679}
"MaConDa"

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

#' A dataframe containing HMDB with unique accurate mass pmd with three digits frequency larger than 1 and accuracy percentage larger than 0.9.
#' @docType data
#' @usage data(hmdb)
#' @format A dataframe with atoms numbers of C, H, O, N, P, S
#' \describe{
#'   \item{percentage}{accuracy of atom numbers prediction}
#'   \item{pmd2}{pmd with two digits}
#'   \item{pmd}{pmd with three digits}
#'   }
"hmdb"

#' A dataframe containing multiple reaction database ID and their related accurate mass pmd and related reactions
#' @docType data
#' @usage data(omics)
#' @format A dataframe with reaction and their realted pmd
#' \describe{
#'   \item{KEGG}{KEGG reaction ID}
#'   \item{RHEA_ID}{RHEA_ID}
#'   \item{DIRECTION}{reaction direction}
#'   \item{MASTER_ID}{master reaction RHEA ID}
#'   \item{ec}{ec reaction ID}
#'   \item{ecocyc}{ecocyc reaction ID}
#'   \item{macie}{macie reaction ID}
#'   \item{metacyc}{metacyc reaction ID}
#'   \item{reactome}{reactome reaction ID}
#'   \item{compounds}{reaction related compounds}
#'   \item{pmd}{pmd with two digits}
#'   \item{pmd2}{pmd with three digits}
#'   }
"omics"

#' A dataframe containing reaction related accurate mass pmd and related reaction formula with KEGG ID
#' @docType data
#' @usage data(keggrall)
#' @format A dataframe with KEGG reaction, their realted pmd and atoms numbers of C, H, O, N, P, S
#' \describe{
#'   \item{ID}{KEGG reaction ID}
#'   \item{pmd}{pmd with three digits}
#'
#'   }
"keggrall"
