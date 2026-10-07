#' A dataset containing common Paired mass distances of substructure, ions replacements, and reaction
#' @docType data
#' @usage data(sda)
#' @format A data frame with 146 rows and 4 variables:
#' \describe{
#'   \item{PMD}{Paired mass distances}
#'   \item{origin}{potential sources}
#'   \item{Ref.}{references}
#'   \item{mode}{b for biological reaction and e for environmental reaction}
#'   }
"sda"

#' Curated database of known multi-step reaction chains
#'
#' A hand-curated set of named biotransformation chains, each an ordered sequence
#' of paired mass distances (PMDs) for known-pathway screening of untargeted
#' feature lists via \code{getchainseq(list, db = pmdchain)}. Each chain spans
#' several rows (one per step). Every step PMD is the exact monoisotopic mass of
#' its CHNOPS element delta. \code{keggcount} is an observable specificity score:
#' the number of compound-connected KEGG reaction paths (from \code{\link{keggrall}})
#' whose PMD signature matches the whole chain -- lower means a rarer, more
#' specific signature (xenobiotic phase-II conjugations are most specific; common
#' edits such as methylation/phosphorylation are highly degenerate). A chain match
#' is a relational annotation, not compound identification.
#' @docType data
#' @usage data(pmdchain)
#' @format A data frame with one row per chain step:
#' \describe{
#'   \item{chain_id}{short identifier of the chain}
#'   \item{name}{human-readable chain name}
#'   \item{class}{natural metabolism or xenobiotic (phase-I/II) chemistry}
#'   \item{nstep}{number of steps in the chain}
#'   \item{step}{1-based step index}
#'   \item{transformation}{the single-step transformation block name}
#'   \item{pmd}{exact monoisotopic mass of the step's element delta (Da)}
#'   \item{dC, dH, dN, dO, dP, dS}{signed CHNOPS element deltas of the step}
#'   \item{keggcount}{number of KEGG reaction paths sharing the chain's PMD signature (specificity proxy; lower = rarer)}
#'   \item{description}{biochemical description of the chain}
#'   }
"pmdchain"

#' mass spectrometry contaminants database for PMD check
#' @docType data
#' @usage data(MaConDa)
#' @format A data frame from \doi{doi:10.1093/bioinformatics/bts527} with 308 rows and 5 variables:
#' \describe{
#'   \item{id}{MaConDa ID}
#'   \item{name}{contaminants}
#'   \item{formula}{contaminants fomula}
#'   \item{exact_mass}{exact mass of contaminants}
#'   \item{type_of_contaminant}{type of contaminant}
#'   }
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

