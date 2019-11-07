#' Reading MSGF output from PNNL's DMS
#'
#' @param DataPkgNumber (Numeric or Character vector) containing Data Package ID(s) located in DMS
#' @return (MSnID) MSnID object
#' @importFrom dplyr mutate
#' @importFrom MSnID MSnID
#' @importMethodsFrom MSnID psms<-
#' @examples
#' msnid <- read_msgf_data_from_DMS(3442)
#' print(msnid)
#' head(MSnID::psms(msnid))

#' @export
read_msms_data_from_DMS <- function(DataPkgNumber){
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

}

