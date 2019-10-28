#' Reading MASIC Results. Generic.
#'
#' Reading MASIC output from a single directory
#'
#' @param path_to_MASIC_results (path string) to directory with MASIC results for all datasets
#' @param extra_metrics (logical) fetch extra metrics that MASIC extracts from dataset or not. Default is FALSE.
#' @return (data.frame) with reporter ion intensities and other metrics
#' @importFrom dplyr inner_join
#' @examples
#' path_to_MASIC_results <- system.file("extdata/global/masic_output", package = "PlexedPiperTestData")
#' x <- read_masic_data(path_to_MASIC_results, extra_metrics=TRUE)
#' head(x)
#'
#' @export
read_masic_data <- function(path_to_MASIC_results, extra_metrics=FALSE){
   out <- collate_files(path_to_MASIC_results, "_ReporterIons.txt")
   if(extra_metrics){
      extra <- collate_files(path_to_MASIC_results, "_SICstats.txt")
      out <- inner_join(out, extra,
                        by=c("Dataset", "ScanNumber" = "FragScanNumber"))
   }
   return(out)
}
