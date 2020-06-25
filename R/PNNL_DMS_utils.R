#' Utilities for interacting with PNNL DMS
#'
#' @description
#' Fetches results of different data analysis tools.
#' Not useful outside of PNNL unless connected through VPN.
#' Works on windows out of the box. To make it work on Mac/Linux FreeTDS needs
#' to be installed. Mac/Linux functionality is not well tested.
#'
#' @description
#' * `get_server_name_for_mtdb()`: returns the server name it is stored at given mass tag database
#' * `get_dms_job_records()`: returns job records given variety keyword patterns
#' * `get_tool_output_files_for_job_number()`: returns the output of the tool given tool name and file pattern
#' * `get_output_folder_for_job_and_tool()`: returns the path given job id and tool name
#'
#' * `get_AScore_results()`: returns Ascore results given data packge number
#' * `get_job_records_by_dataset_package()`: returns job records given data packge number
#' * `get_results_for_multiple_jobs()`: returns concatenated results given job numbers
#' * `get_results_for_multiple_jobs.dt()`: returns results as concatenated data.table given job numbers
#' * `get_results_for_single_job()`: returns results given job number
#' * `get_results_for_single_job.dt()`: returns results as data.table given job number
#' @md
#'
#'
#' @param mtdbName (string) name of mass tag database
#' @param dataPkgNumber (interger) data package id
#' @param jobs (integer) DMS job ID
#' @param xPttrn (string) substrings for names of dataset, experiment, tool, parameter file, settings file, fasta file, protein options, instrument name
#' @param mostRecent (logical) only most recent or all output files
#'
#' @importFrom RODBC odbcDriverConnect sqlQuery
# '@importFrom data.table
#'
#' @name pnnl_dms_utils
#'
#' @examples
#' get_output_folder_for_job_and_tool(863951, "DTA_Refinery")
NULL

# dictionary that defines the suffix of the files given the analysis tool
tool2suffix = list("MSGFPlus"="_msgfplus_syn.txt",
                   "MSGFPlus_MzML"="_msgfplus_syn.txt",
                   "MSGFPlus_DTARefinery"="_msgfplus_syn.txt",
                   "MSGFDB_DTARefinery"="_msgfdb_syn.txt",
                   "MASIC_Finnigan"="_ReporterIons.txt")

get_driver <- function(){
   if(.Platform$OS.type == "unix"){
      return("FreeTDS")
   }else if(.Platform$OS.type == "windows"){
      return("SQL Server")
   }else{
      stop("Unknown OS type.")
   }
}

get_auth <- function(){
   if(.Platform$OS.type == "unix"){
      return("PORT=1433;UID=dmsreader;PWD=dms4fun;")
   }else if(.Platform$OS.type == "windows"){
      return("")
   }else{
      stop("Unknown OS type.")
   }
}



#' @export
#' @rdname pnnl_dms_utils
get_server_name_for_mtdb = function(mtdbName)
   # a way to find which server MT DB is located
{
   con_str <- sprintf("DRIVER={%s};SERVER=Pogo;DATABASE=MTS_Master;%s",
                      get_driver(),
                      get_auth())
   con <- odbcDriverConnect(con_str)
   strSQL = sprintf("SELECT Server_Name
                    FROM V_MTS_MT_DBs
                    WHERE MT_DB_Name = '%s'", mtdbName)
   dbServer = sqlQuery(con,strSQL)
   close(con)
   return(as.character(dbServer[1,1]))
}

#' @export
#' @rdname pnnl_dms_utils
get_dms_job_records = function(
   jobs = NULL,
   datasetPttrn = "",
   experimentPttrn = "",
   toolPttrn = "",
   parPttrn = "",
   settingsPttrn = "",
   fastaPttrn = "",
   proteinOptionsPttrn = "",
   intrumentPttrn = ""){

   # first check if the input is valid
   x = as.list(environment())
   x[["jobs"]] = NULL
   if( all(x == "") & is.null(jobs) ){
      stop("insufficients arguments provided")
   }
   if( any(x != "") & !is.null(jobs) ){
      stop("can't provide both: job list and search terms")
   }

   # initialize connection
   con_str <- sprintf("DRIVER={%s};SERVER=gigasax;DATABASE=dms5;%s",
                      get_driver(),
                      get_auth())
   con <- odbcDriverConnect(con_str)

   # set-up query based on job list
   if(!is.null(jobs)){
      strSQL = sprintf("SELECT *
                       FROM V_Mage_Analysis_Jobs
                       WHERE [Job] IN ('%s')
                       ",
                       paste(jobs,sep="",collapse="',\n'"))
   }else{
      strSQL = sprintf("SELECT *
                       FROM V_Mage_Analysis_Jobs
                       WHERE [Dataset] LIKE '%%%s%%'
                       AND [Experiment] LIKE '%%%s%%'
                       AND [Tool] LIKE '%%%s%%'
                       AND [Parameter_File] LIKE '%%%s%%'
                       AND [Settings_File] LIKE '%%%s%%'
                       AND [Protein Collection List] LIKE '%%%s%%'
                       AND [Protein Options] LIKE '%%%s%%'
                       AND [Instrument] LIKE '%%%s%%'
                       ",
                       datasetPttrn,
                       experimentPttrn,
                       toolPttrn,
                       parPttrn,
                       settingsPttrn,
                       fastaPttrn,
                       proteinOptionsPttrn,
                       intrumentPttrn)
   }
   locationPointersToMSMSjobs = sqlQuery(con, strSQL, stringsAsFactors=FALSE)
   close(con)
   return(locationPointersToMSMSjobs)
}


#' @export
#' @rdname pnnl_dms_utils
get_tool_output_files_for_job_number <- function(jobNumber, toolName,
                                                 filePattern, mostRecent=TRUE)
{
   # get job records first. This will be useful to get dataset folder
   jobRecord <- get_dms_job_records(jobNumber)
   datasetFolder <- dirname( as.character(jobRecord$Folder))

   # get tool's subfolder
   if( is.missing(toolName) ){
      toolFolder = ''
   }else{
      # return stuff from the main dataset folder
      toolFolder = get_output_folder_for_job_and_tool(jobNumber, toolName, mostRecent)
   }
   #
   candidateFiles = list.files(file.path(datasetFolder, toolFolder),
                               pattern=filePattern, full.names=TRUE,
                               ignore.case=TRUE)
   #
   if(length(candidateFiles) == 1){
      return(candidateFiles)
   }else{
      return(NA)
   }
}


#' @export
#' @rdname pnnl_dms_utils
get_output_folder_for_job_and_tool <- function(jobNumber, toolName, mostRecent=TRUE)
{
   con_str <- sprintf("DRIVER={%s};SERVER=gigasax;DATABASE=DMS_Pipeline;%s",
                      get_driver(),
                      get_auth())
   con <- odbcDriverConnect(con_str)

   strSQLPattern = "SELECT Output_Folder
   FROM V_Job_Steps_History
   WHERE (Job = %s) AND (Tool = '%s') AND (Most_Recent_Entry = 1)"
   strSQL = sprintf( strSQLPattern, jobNumber, toolName)
   qry = sqlQuery(con, strSQL)
   close(con)
   return(as.character(qry[1,1]))
}


#' @export
#' @rdname pnnl_dms_utils
# Get AScore results for a given data package
get_AScore_results <- function(dataPkgNumber)
{
   #
   con_str <- sprintf("DRIVER={%s};SERVER=gigasax;DATABASE=dms5;%s",
                      get_driver(),
                      get_auth())
   con <- odbcDriverConnect(con_str)
   strSQL <- sprintf("SELECT *
                     FROM V_Mage_Analysis_Jobs
                     WHERE (Dataset LIKE 'DataPackage_%s%%')", dataPkgNumber)
   jobs <- sqlQuery(con, strSQL, stringsAsFactors=FALSE)
   close(con)
   #
   if(nrow(jobs) == 1){
      # library("RSQLite")
      dlist <- dir(as.character(jobs["Folder"]))
      idx <- which.max(as.numeric(sub("Step_(\\d+)_.*", "\\1", dlist))) # the Results supposed to be in the last folder
      ascoreResultDB <- file.path( jobs["Folder"], dlist[idx], "Results.db3")
      db <- dbConnect(SQLite(), dbname = ascoreResultDB)
      AScores <- dbGetQuery(db, "SELECT * FROM t_results_ascore")
      dbDisconnect(db)
      return(AScores)
   }else{
      return(NULL)
   }
}


#' @export
#' @rdname pnnl_dms_utils
get_job_records_by_dataset_package <- function(dataPkgNumber)
{
   con_str <- sprintf("DRIVER={%s};SERVER=gigasax;DATABASE=dms5;%s",
                      get_driver(),
                      get_auth())
   con <- odbcDriverConnect(con_str)
   strSQL = sprintf("
                    SELECT *
                    FROM V_Mage_Data_Package_Analysis_Jobs
                    WHERE Data_Package_ID = %s",
                    dataPkgNumber)
   jr <- sqlQuery(con, strSQL, stringsAsFactors=FALSE)
   close(con)
   return(jr)
}


#' @export
#' @rdname pnnl_dms_utils
get_results_for_multiple_jobs = function( jobRecords){
   toolName = unique(jobRecords[["Tool"]])
   if (length(toolName) > 1){
      stop("Contains results of more then one tool.")
   }
   results = ldply( jobRecords[["Folder"]],
                    get_results_for_single_job,
                    fileNamePttrn=tool2suffix[[toolName]],
                    .progress = "text")
   return( results )
}


#' @export
#' @rdname pnnl_dms_utils
get_results_for_multiple_jobs.dt <- function( jobRecords){
   toolName = unique(jobRecords[["Tool"]])
   if (length(toolName) > 1){
      stop("Contains results of more then one tool.")
   }
   results = llply( jobRecords[["Folder"]],
                    get_results_for_single_job.dt,
                    fileNamePttrn=tool2suffix[[toolName]],
                    .progress = "text")
   results.dt <- rbindlist(results)
   return( as.data.frame(results.dt) ) # in the future I may keep it as data.table
}


#' @export
#' @rdname pnnl_dms_utils
get_results_for_single_job = function(pathToFile, fileNamePttrn ){
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
   out = data.frame(Dataset=dataset, results, stringsAsFactors = FALSE)
   return(out)
}


#' @export
#' @rdname pnnl_dms_utils
get_results_for_single_job.dt = function(pathToFile, fileNamePttrn ){
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

#' @export
#' @rdname pnnl_dms_utils
get_study_design_by_dataset_package <- function(dataPkgNumber) {
   
   con_str <- sprintf("DRIVER={%s};SERVER=gigasax;DATABASE=DMS_Data_Package;%s",
                      get_driver(),
                      get_auth())
   con <- dbConnect(odbc(), .connection_string=con_str)
   
   ## fetch table with path to DataPackage
   strSQL <- sprintf("
                    SELECT *
                    FROM V_Data_Package_Detail_Report
                    WHERE ID = %s",
                     dataPkgNumber)
   qry <- dbSendQuery(con, strSQL)
   dataPkgReport <- dbFetch(qry)
   dbClearResult(qry)
   
   if(.Platform$OS.type == "unix"){
      local_folder <- "~/temp_msms_results"
      if(file.exists(local_folder)){
         unlink(local_folder, recursive = T)
      }
      dir.create(local_folder)
      
      remote_folder <- gsub("\\\\","/", dataPkgReport$`Share Path`)
      mount_cmd <- sprintf("mount -t smbfs %s %s", remote_folder, local_folder)
      system(mount_cmd)
   }else if(.Platform$OS.type == "windows"){
      local_folder <- pathToFile
   }else{
      stop("Unknown OS type.")
   }
   
   ## fetch samples.txt
   pathToFile <- list.files(path=local_folder,
                            pattern="^samples\\.txt$",
                            full.names=T)
   if(length(pathToFile) == 0){
      stop("Could not find 'samples.txt' file in DMS Data Package folder.")
   }
   
   samples <- readr::read_tsv(pathToFile, col_types=readr::cols(), progress=FALSE)
   if (!setequal(colnames(samples), c("PlexID",
                                      "QuantBlock",
                                      "ReporterAlias",
                                      "ReporterName",
                                      "MeasurementName"))) {
      stop(paste0("Please check the required columns in the 'samples.txt' file.\n",
                  paste0(c("PlexID",
                           "QuantBlock",
                           "ReporterAlias",
                           "ReporterName",
                           "MeasurementName"), collapse = "\n")))
   }
   
   ## fetch fractions.txt
   pathToFile <- list.files(path=local_folder,
                            pattern="^fractions\\.txt$",
                            full.names=T)
   if (length(pathToFile) == 0){
      stop("Could not find 'fractions.txt' file in DMS Data Package folder.")
   }
   
   
   fractions <- readr::read_tsv(pathToFile, col_types=readr::cols(), progress=FALSE)
   if (!setequal(colnames(fractions), c("PlexID",
                                        "Dataset"))) {
      stop(paste0("Please check the required columns in the 'fractions.txt' file.\n",
                  paste0(c("PlexID",
                           "Dataset"), collapse = "\n")))
   }
   
   ## fetch references.txt
   pathToFile <- list.files(path=local_folder,
                            pattern="^references\\.txt$",
                            full.names=T)
   if(length(pathToFile) == 0){
      stop("Could not find 'references.txt' file in DMS Data Package folder.")
   }
   
   references <- readr::read_tsv(pathToFile, col_types=readr::cols(), progress=FALSE)
   if (!setequal(colnames(references), c("PlexID",
                                         "QuantBlock",
                                         "Reference"))) {
      stop(paste0("Please check the required columns in the 'references.txt' file.\n",
                  paste0(c("PlexID",
                           "QuantBlock",
                           "Reference"), collapse = "\n")))
   }
   
   
   
   if(.Platform$OS.type == "unix"){
      umount_cmd <- sprintf("umount %s", local_folder)
      system(umount_cmd)
      unlink(local_folder, recursive = T)
   }
   
   study_des <- list(samples = samples, 
                     fractions = fractions, 
                     references = references)
   return(study_des)
}
