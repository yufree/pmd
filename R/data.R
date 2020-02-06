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

#' A dataframe containing HMDB with unique accurate mass pmd frequence larger than 100.
#' @docType data
#' @usage data(hmdb)
#' @format A dataframe with atoms numbers of C, H, O, N, P, S
#' \describe{
#'   \item{percentage}{accuracy of atom numbers prediction}
#'   \item{pmd}{pmd with two digits}
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
#'   }
"keggrall"

#' A list containing HMDB qqq MS/MS data with peaks larger than 10 percentage for PMD annatation
#' @docType data
#' @usage data(qqq)
#' @format A list containing HMDB qqq MS/MS data with peaks larger than 10 percentage for PMD annatation
#' \describe{
#'   \item{name}{HMDB ID}
#'   \item{mz}{mass to charge ratio}
#'   \item{msms}{msms pmd}
#'   \item{msmsraw}{raw msms data}
#'   }
"qqq"

#' A list containing HMDB qtof MS/MS data with peaks larger than 10 percentage for PMD annatation
#' @docType data
#' @usage data(qtof)
#' @format A list containing HMDB qqq MS/MS data with peaks larger than 10 percentage for PMD annatation
#' \describe{
#'   \item{name}{HMDB ID}
#'   \item{mz}{mass to charge ratio}
#'   \item{msms}{msms pmd}
#'   \item{msmsraw}{raw msms data}
#'   }
"qtof"

#' A list containing HMDB orbitrap MS/MS data with peaks larger than 10 percentage for PMD annatation
#' @docType data
#' @usage data(orb)
#' @format A list containing HMDB orbitrap MS/MS data with peaks larger than 10 percentage for PMD annatation
#' \describe{
#'   \item{name}{HMDB ID}
#'   \item{mz}{mass to charge ratio}
#'   \item{msms}{msms pmd}
#'   \item{msmsraw}{raw msms data}
#'   }
"orb"
