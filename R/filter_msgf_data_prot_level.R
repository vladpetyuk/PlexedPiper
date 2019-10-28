#' Filtering MSGF Data
#'
#' Filtering MSGF data.
#'
#' @param msnid (MSnID object) collated MSGF output
#' @param protein_FDR_threshold (numeric) Maximum acceptable FDR rate. Default is 0.01.
#' @param n.iter.grid (numeric) number of grid-distributed evaluation points. Default 500.
#' @param n.iter.nm (numeric) number of iterations for Nelder-Mead optimization algorithm. Default 100.
#' @return (MSnID object) filtered MSGF output
#' @importFrom MSnID MSnIDFilter MSnIDFilter optimize_filter apply_filter
#' @export filter_msgf_data_protein_level
#' @examples
#' path_to_MSGF_results <- system.file("extdata/global/msgf_output", package = "PlexedPiperTestData")
#' msnid <- read_msgf_data(path_to_MSGF_results)
#' msnid <- MSnID::correct_peak_selection(msnid)
#' show(msnid)
#' msnid <- filter_msgf_data_peptide_level(msnid)
#' show(msnid)
#' path_to_FASTA <- system.file("extdata/Rattus_norvegicus_NCBI_RefSeq_2018-04-10.fasta.gz", package = "PlexedPiperTestData")
#' msnid <- compute_num_peptides_per_1000aa(msnid, path_to_FASTA)
#' msnid <- filter_msgf_data_protein_level(msnid)
#' show(msnid)
#'

filter_msgf_data_protein_level <- function(msnid,
                                           protein_FDR_threshold=0.01,
                                           n.iter.grid=500,
                                           n.iter.nm=100){
   # setup
   filtObj <- MSnIDFilter(msnid)
   filtObj$peptides_per_1000aa <- list(comparison=">", threshold=1.0)
   # step 1
   filtObj.grid <- optimize_filter(filtObj,
                                   msnid,
                                   fdr.max=0.01,
                                   method="Grid",
                                   level="accession",
                                   n.iter=n.iter.grid)
   # step 2
   filtObj.nm <- optimize_filter(filtObj.grid,
                                 msnid,
                                 fdr.max=0.01,
                                 # method="Nelder-Mead",
                                 method="SANN",
                                 level="accession",
                                 n.iter=n.iter.nm)
   apply_filter(msnid, filtObj.nm)
}

