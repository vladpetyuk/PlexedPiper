#' Filtering MSGF Data
#'
#' Filtering MSGF data. In this implementation the function optmized both ppm and
#' PepQValue thresholds to achieve maximum number of peptide identifications within
#' given FDR constrain.
#'
#' @param msnid (MSnID object) collated MSGF output
#' @param peptide_FDR_threshold (numeric) Maximum acceptable FDR rate. Default is 0.01.
#' @param n.iter.grid (numeric) number of grid-distributed evaluation points. Default 500.
#' @param n.iter.nm (numeric) number of iterations for Nelder-Mead optimization algorithm. Default 100.
#' @return (MSnID object) filtered MSGF output
#' @importFrom MSnID MSnIDFilter MSnIDFilter optimize_filter mass_measurement_error apply_filter
#' @export filter_msgf_data_peptide_level
#' @examples
#' path_to_MSGF_results <- system.file("extdata/global/msgf_output", package = "PlexedPiperTestData")
#' msnid <- read_msgf_data(path_to_MSGF_results)
#' msnid <- MSnID::correct_peak_selection(msnid)
#' show(msnid)
#' msnid <- filter_msgf_data_peptide_level(msnid, 0.01) # 1% FDR
#' show(msnid)

filter_msgf_data_peptide_level <- function(msnid,
                                           peptide_FDR_threshold=0.01,
                                           n.iter.grid=500,
                                           n.iter.nm=100){
   # setup
   msnid$msmsScore <- -log10(msnid$PepQValue)
   msnid$absParentMassErrorPPM <- abs(mass_measurement_error(msnid))
   filtObj <- MSnIDFilter(msnid)
   filtObj$absParentMassErrorPPM <- list(comparison="<", threshold=10.0)
   filtObj$msmsScore <- list(comparison=">", threshold=2.0)
   # step 1
   filtObj.grid <- optimize_filter(filtObj,
                                   msnid,
                                   fdr.max=0.01,
                                   method="Grid",
                                   level="peptide",
                                   n.iter=n.iter.grid)
   # step 2
   filtObj.nm <- optimize_filter(filtObj.grid,
                                 msnid,
                                 fdr.max=0.01,
                                 method="Nelder-Mead",
                                 level="peptide",
                                 n.iter=n.iter.nm)
   apply_filter(msnid, filtObj.nm)
}

