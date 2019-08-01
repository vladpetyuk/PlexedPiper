#' Makes Cross-Tab With Quantitative Data.
#'
#' Links filtered MS/MS IDs with reporter intensities. Divides reporter ion
#' intensities by corresponding reference and then returns cross-tab. Rows are
#' species (e.g. proteins in global, phosphopeptides in phosphoproteomic experiment),
#' columns are sample names.
#'
#' @param msnid (MSnID object) filtered MS/MS identifications
#' @param reporter_intensities (data.frame) collated table with filtered reporter
#' intensities.
#' @param aggregation_level (string vector) defines what intensities needs to be aggregated.
#' At this point the only aggregation function is `sum`.
#' Typically intensities from different fractions of the same plex are aggregated.
#' Also e.g. in global proteomics intensities from different scans identifiying peptides
#' from the same protein aggregated togeher too.
#' @param converter (data.frame) A conversion table to covert TMT channels to 
#' reporter ions calculated by MASIC. If not provided, will default to TMT10
#' table provided in `PlexedPiperTestData` package called `reporter_converter`.
#' @return (matrix) with log2-transformed relative reporter ion intensities.
#' Row names are the names of the measured species.
#' Column names are the names of the samples.
#' @importFrom data.table data.table setkey setkeyv melt dcast
#' @export create_crosstab
#' @examples
#'
#' # Prepare MS/MS IDs
#' path_to_MSGF_results <- system.file("extdata/global/msgf_output", package = "PlexedPiperTestData")
#' msnid <- read_msgf_data(path_to_MSGF_results)
#' msnid <- MSnID::correct_peak_selection(msnid)
#' show(msnid)
#' msnid <- filter_msgf_data_peptide_level(msnid, 0.01)
#' show(msnid)
#' path_to_FASTA <- system.file("extdata/Rattus_norvegicus_NCBI_RefSeq_2018-04-10.fasta.gz", package = "PlexedPiperTestData")
#' msnid <- compute_num_peptides_per_1000aa(msnid, path_to_FASTA)
#' msnid <- filter_msgf_data_protein_level(msnid, 0.01)
#' show(msnid)
#'
#' # Prepare table with reporter ion intensities
#' path_to_MASIC_results <- system.file("extdata/global/masic_output", package = "PlexedPiperTestData")
#' masic_data <- read_masic_data(path_to_MASIC_results, extra_metrics=TRUE)
#' masic_data <- filter_masic_data(masic_data, 0.5, 0)
#'
#' # Creating cross-tab
#' library(readr)
#' fractions <- read_tsv(system.file("extdata/study_design/fractions.txt", package = "PlexedPiperTestData"))
#' samples <- read_tsv(system.file("extdata/study_design/samples.txt", package = "PlexedPiperTestData"))
#' references <- read_tsv(system.file("extdata/study_design/references.txt", package = "PlexedPiperTestData"))
#' aggregation_level <- c("PlexID","accession")
#'
#' out <- create_crosstab(msnid, masic_data, aggregation_level, fractions, samples, references)
#' dim(out)
#' head(out)

create_crosstab <- function(msnid,
                            reporter_intensities,
                            aggregation_level,
                            fractions,
                            samples,
                            references,
                            converter){

   # merges MS/MS IDs with reporter intensities
   quant_data <- link_msms_and_reporter_intensities(msnid, reporter_intensities)
   # aggregates reporter intensities to a given level
   quant_data <- aggregate_reporters(quant_data, fractions, aggregation_level)
   # taking ratios of reporter ion intensities to whatever the reference is
   out <- converting_to_relative_to_reference(quant_data, 
                                              samples, 
                                              references, 
                                              converter)
   return(out)
}




# no export
link_msms_and_reporter_intensities <- function(msnid, reporter_intensities)
{
   # prepare MS/MS data
   msms <- data.table(psms(msnid))
   msms <- msms[,.(Dataset, Scan, peptide, accession)]
   setkey(msms, Dataset, Scan)

   # prepare reporter intensities data
   # check that there are no extra columns
   stopifnot(all(c("Dataset","ScanNumber") %in% colnames(reporter_intensities)))
   other_columns <- setdiff(colnames(reporter_intensities), c("Dataset","ScanNumber"))
   stopifnot(all(grepl("^Ion.*\\d$",other_columns)))
   reporter_intensities <- data.table(reporter_intensities)
   reporter_intensities[,Scan := ScanNumber][,ScanNumber := NULL]

   # main link
   quant_data <- merge(msms, reporter_intensities)

   return(quant_data)
}

# no export
aggregate_reporters <- function(quant_data, fractions, aggregation_level)
{
   aggregation_level <- c("PlexID", aggregation_level)

   fractions <- data.table(fractions, key="Dataset")
   quant_data <- merge(quant_data, fractions)

   setkeyv(quant_data, aggregation_level)
   quant_data <- quant_data[,lapply(.SD,sum,na.rm=TRUE),
                            by=aggregation_level,
                            .SDcols=grep("^Ion_1.*\\d$", colnames(quant_data), value=T)]

   #
   specie_cols <- setdiff(aggregation_level, "PlexID")
   quant_data[, Specie := do.call(paste, c(.SD, sep = "@")), .SDcols=specie_cols]
   quant_data[,(specie_cols) := NULL]

   return(quant_data)
}

# no export
converting_to_relative_to_reference <- function(quant_data, 
                                                samples, 
                                                references,
                                                converter)
{

   # transforming from wide to long form
   quant_data <- melt(quant_data,
                     id.vars=c("PlexID","Specie"),
                     measure.vars=grep("Ion_", colnames(quant_data), value=TRUE),
                     variable.name="ReporterIon",
                     value.name="Abundance")
   setkey(quant_data, PlexID, ReporterIon)

   # preparing sample info
   samples <- data.table(samples)
   if (is.null(converter)) {
     converter <- data.table(reporter_converter)
   }
   else {
     converter <- data.table(converter)
   }
   samples <- merge(samples, converter)
   setkey(samples, PlexID, ReporterIon)

   # merging reporter intensities with sample info
   quant_data <- merge(quant_data, samples)
   setkey(quant_data, PlexID, QuantBlock)

   out <- list()
   #~~~
   for(i in seq_len(nrow(references))){

      # unique PlexID/QuantBlock combo
      ref_i <- data.table(references[i,])
      setkey(ref_i, PlexID, QuantBlock)

      # subset a piece that is unique to that reference
      quant_data_i <- quant_data[ref_i]

      # now make wide by ReporterAlias
      quant_data_i_w <- dcast(quant_data_i,
                              Specie ~ ReporterAlias,
                              value.var='Abundance')
      species <- quant_data_i_w[,"Specie"]
      quant_data_i_w <- quant_data_i_w[,!"Specie"]

      # compute reference values for this particular PlexID/QuantBlock combo
      # TODO. Maybe there is a more elegant, more data.table-esque way.
      if (ref_i$Reference == "ref") {
        ref_values <- quant_data_i_w$ref
      }
      else {
        ref_values <- with(quant_data_i_w, eval(parse(text=ref_i$Reference)))
      }
      
      # take the ratios over the reference
      quant_data_i_w <- quant_data_i_w/ref_values

      # switching from ReporterAlias to MeasurementName
      quant_data_i_w <- data.table(species, quant_data_i_w)
      quant_data_i_l <- melt(quant_data_i_w, id.vars="Specie",
                             variable.name="ReporterAlias",
                             value.name="Ratio")

      setkey(samples, PlexID, QuantBlock)
      sample_naming <- samples[ref_i
                               ][,.(ReporterAlias,MeasurementName)
                                 ][!is.na(MeasurementName)]

      # converting sample names
      setkey(quant_data_i_l, ReporterAlias)
      setkey(sample_naming, ReporterAlias)
      quant_data_i_l <- merge(quant_data_i_l, sample_naming)
      quant_data_i_l[,ReporterAlias := NULL]
      out <- c(out, list(quant_data_i_l))
   }
   #~~~
   out <- rbindlist(out)
   out <- dcast(out, Specie ~ MeasurementName, value.var="Ratio")
   species <- out[,Specie]
   out <- as.matrix(out[,!"Specie"])
   rownames(out) <- species

   # Inf and 0 values turn into NA
   out[is.infinite(out)] <- NA
   out[out == 0] <- NA

   # log2-transform
   out <- log2(out)
   return(out)
}


