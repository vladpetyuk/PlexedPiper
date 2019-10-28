#' Converter for TMT Reporter Channels
#'
#' Converter from human-readable notation of TMT reporter channels to
#' m/z-type (reported by MASIC).
#'
#' @format A data frame with 10 rows and 2 variables:
#' \describe{
#'   \item{ReporterName}{human-readable notation}
#'   \item{ReporterIon}{reporter channel names extracted by MASIC}
#'   ...
#' }
#' @source \url{https://www.thermofisher.com/us/en/home/life-science/protein-biology/protein-mass-spectrometry-analysis/protein-quantitation-mass-spectrometry/tandem-mass-tag-systems.html}
"reporter_converter"
