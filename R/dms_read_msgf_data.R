#' Reading MSGF Results. Generic.
#'
#' Reading MSGF output from PNNL's DMS
#'
#' @param path_to_MSGF_results (path string) to directory with MSGF results for all datasets
#' @param DataPkgNumber (Numeric or Character vector) containing Data Package ID(s) located in DMS
#' @return (MSnID) MSnID object
#' @importFrom dplyr mutate
#' @importFrom MSnID MSnID
#' @importMethodsFrom MSnID psms<-
#' @examples
#' path_to_MSGF_results <- system.file("extdata/global/msgf_output", package = "PlexedPiperTestData")
#' msnid <- read_msgf_data(path_to_MSGF_results)
#' print(msnid)
#' head(MSnID::psms(msnid))

#' @export
dms_read_msms_data <- function(path_to_MSGF_results, DataPkgNumber = NULL){
  msnid <- MSnID(".")
  # accession -> Protein
  # calculatedMassToCharge -> f(MH, Charge) MSnID:::.PROTON_MASS
  # (MH + (Charge-1)*MSnID:::.PROTON_MASS)/Charge
  # chargeState -> Charge
  # experimentalMassToCharge -> PrecursorMZ
  # isDecoy -> grep("^XXX", Protein)
  # peptide -> Peptide
  # spectrumFile -> Dataset
  # spectrumID -> Scan
  
  if (!is.null(DataPkgNumber)) {
    
    # Fetch job records for data package(s)
    if (length(DataPkgNumber) > 1) {
      job_rec_ls <- lapply(DataPkgNumber, get_job_records_by_dataset_package)
      jobRecords <- Reduce(rbind, job_rec_ls)
    }
    
    else {
      jobRecords <- get_job_records_by_dataset_package(DataPkgNumber)
    }
    
    jobRecords <- jobRecords[grepl("MSGFPlus", jobRecords$Tool),]
    
    x <- get_results_for_multiple_jobs.dt(jobRecords)
    x <- x %>% 
      mutate(accession = Protein,
             calculatedMassToCharge = (MH + (Charge-1)*MSnID:::.PROTON_MASS)/Charge,
             chargeState = Charge,
             experimentalMassToCharge = PrecursorMZ,
             isDecoy = grepl("^XXX", Protein),
             peptide = Peptide,
             spectrumFile = Dataset,
             spectrumID = Scan)
    
    x <- mutate(x, pepSeq = MSnID:::.get_clean_peptide_sequence(peptide))
    
    psms(msnid) <- x
    return(msnid)
  }
  
  else{
    x <- collate_files(path_to_MSGF_results, "_msgfplus_syn.txt") %>%
      mutate(accession = Protein,
             calculatedMassToCharge = (MH + (Charge-1)*MSnID:::.PROTON_MASS)/Charge,
             chargeState = Charge,
             experimentalMassToCharge = PrecursorMZ,
             isDecoy = grepl("^XXX", Protein),
             peptide = Peptide,
             spectrumFile = Dataset,
             spectrumID = Scan)
    # extra piece that may be removed in the future.
    # decoy accessions stripped of XXX_ prefix
    # x <- x %>%
    #    mutate(accession = sub("^XXX_","",accession))
    
    # clean peptide sequences
    x <- mutate(x, pepSeq = MSnID:::.get_clean_peptide_sequence(peptide))
    
    psms(msnid) <- x
    
    return(msnid)
  }
}


#' # (yet) non-exported helper function
#' #' @importFrom Biostrings AA_STANDARD
#' #' @importFrom purrr map_chr
#' #' @importFrom stringr str_replace_all
#' peptide_to_sequence <- function(pep){
#'    pep_no_flank <- sub(".\\.(.*)\\..", "\\1", pep)
#'    present_chars <- paste0(pep_no_flank, collapse = '') %>%
#'       strsplit(split='') %>%
#'       `[[`(1) %>%
#'       unique()
#'    other_chars <- setdiff(present_chars, AA_STANDARD)
#'    if(length(other_chars) > 0){
#'       message("Detected extra chararacters in the peptide sequences!")
#'       # erase other chars in TrimmedPeptide
#'       other_chars_pttrn <- other_chars %>%
#'          map_chr(~paste0("\\",.x)) %>%
#'          paste0(collapse='') %>%
#'          paste0("[",.,"]")
#'       pep_clean <- str_replace_all(pep_no_flank, other_chars_pttrn, "")
#'    }
#'    return(pep_clean)
#' }
