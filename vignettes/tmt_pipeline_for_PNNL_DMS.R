## ----setup, echo=FALSE-----------------------------------------------------
knitr::opts_chunk$set(message=FALSE, warning=FALSE)
# knitr::opts_chunk$set(echo=T, message=F, warning=F, fig.align='center', out.width='10cm')

## ----libraries-------------------------------------------------------------
library(PlexedPiper)
library(MSnID)
library(tidyverse)

## ----msms, eval=FALSE------------------------------------------------------
#  msms <- dms_read_msms_data(3442)
#  
#  con_str <- sprintf("DRIVER={%s};SERVER=gigasax;DATABASE=dms5;%s",
#          "FreeTDS", "PORT=1433;UID=dmsreader;PWD=dms4fun;")
#  con_str <- "DRIVER={FreeTDS};SERVER=gigasax;DATABASE=dms5;PORT=1433;UID=dmsreader;PWD=dms4fun"
#  con_str <- "DRIVER={};SERVER=gigasax;DATABASE=dms5;PORT=1433;UID=dmsreader;PWD=dms4fun"
#  library(RODBC)
#  con <- odbcDriverConnect(con_str)
#  
#  library(odbc)
#  con_str <- "DRIVER={FreeTDS};SERVER=gigasax;DATABASE=dms5;PORT=1433;UID=dmsreader;PWD=dms4fun"
#  con <- dbConnect(odbc::odbc(), .connection_string=con_str)
#  con <- dbConnect(odbc::odbc(), .connection_string=con_str)
#  
#  
#  strSQL = sprintf("\n                    SELECT *\n                    FROM V_Mage_Data_Package_Analysis_Jobs\n                    WHERE Data_Package_ID = %s",
#          3442)
#  jr <- dbSendQuery(con, strSQL, stringsAsFactors = FALSE)
#  stuff <- dbFetch(jr)
#  
#  
#  
#  pth <- "\\\\proto-3\\QExactHF03\\2017_4\\MoTrPAC_Pilot_TMT_W_S1_01_12Oct17_Elm_AQ-17-09-02\\MSG201804111041_Auto1572014"
#  x <- get_results_for_single_job.dt(pth, "_msgfplus_syn.txt")
#  
#  
#  # mount smbfs -t //proto-3/QExactHF03/2017_4/MoTrPAC_Pilot_TMT_W_S1_01_12Oct17_Elm_AQ-17-09-02/MSG201804111041_Auto1572014 ~/temp2
#  
#  
#  

