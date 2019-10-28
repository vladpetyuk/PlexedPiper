#' Reading MASIC Results. PNNL system only.
#'
#' Reading MASIC output from PNNL's DMS
#'
#' @param path_to_MASIC_results (path string) to directory with MASIC results for all datasets
#' @param extra_metrics (logical) fetch extra metrics that MASIC extracts from dataset or not. Default is FALSE.
#' @return (data.frame) with reporter ion intensities and other metrics
#' @importFrom plyr llply
#' @importFrom data.table rbindlist
#' @export fetch_masic_data

fetch_masic_data <- function(path_to_MASIC_results, extra_metrics=FALSE){
   library("plyr")
   library("data.table")
   results = llply( path_to_MASIC_results,
                    fetch_masic_data_for_single_datset,
                    fileNamePttrn=tool2suffix[["MASIC_Finnigan"]],
                    .progress = "text")
   results.dt <- rbindlist(results)
   return(as.data.frame(results.dt))
}


fetch_masic_data_for_single_datset <- function(pathToFile, fileNamePttrn ){
   pathToFile = list.files( path=as.character(pathToFile),
                            pattern=fileNamePttrn,
                            full.names=T)
   if(length(pathToFile) == 0){
      stop("can't find the results file")
   }
   if(length(pathToFile) > 1){
      stop("ambiguous results files")
   }
   results = read.delim( pathToFile, header=T, stringsAsFactors = FALSE)
   dataset = strsplit( basename(pathToFile), split=fileNamePttrn)[[1]]
   out = data.table(Dataset=dataset, results)
   return(out)
}





#
#
#
#
# library("devtools")
# source_url("https://raw.githubusercontent.com/vladpetyuk/PNNL_misc/master/PNNL_DMS_utils.R")
#
#
# # 1
# jobRecords <- lapply(c(2930,2988,3087), get_job_records_by_dataset_package)
# jobRecords <- Reduce(rbind, jobRecords)
# jobRecords <- subset(jobRecords, Tool == "MASIC_Finnigan")
#
#
#
# system.time({
#     masicData <- get_results_for_multiple_jobs.dt(jobRecords)
# })
#
#
# library(plyr)
# library(data.table)
# results = llply( jobRecords[["Folder"]],
#                  get_results_for_single_job.dt,
#                  fileNamePattern="_SICstats.txt",
#                  .progress = "text")
# results.dt <- rbindlist(results)
# masicStats <- as.data.frame(results.dt)
#
# # hack to remove redundant Dataset column
# masicData  <- masicData[,-2] # use make.names and then remove Dataset.1
# masicStats <- masicStats[,-2]
# colnames(masicStats)[colnames(masicStats) == "FragScanNumber"] <- "ScanNumber" # dplyr::rename
#
# library(dplyr)
#
# x <- select(masicData, Dataset, ScanNumber, starts_with("Ion"), -contains("Resolution"))
# y <- select(masicStats, Dataset, ScanNumber, InterferenceScore)
# z <- inner_join(x, y)
# masicData <- z
# save(masicData, file="masicData_original.RData")
#
#
#
#
#
#
