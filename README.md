# PlexedPiper

R package used at PNNL for processing isobaric labeling (e.g. TMT) proteomics
data. The key inputs are:

* MS/MS identifications from MS-GF+ search engine
* Reporter ion intensities extracted usign MASIC tool
* Tables outlying study design
   * table linking dataset to plexes
   * table linking reporter channels with sample names
   * table identifying reference within each plex

```{r}
if(!require("remotes", quietly = T)) install.packages("remotes")
remotes::install_github("vladpetyuk/PlexedPiper", build_vignettes = TRUE)
library(PlexedPiper)
vignette("tmt_pipeline_v1")
```

Compantion R package with test data based on MoTrPAC pilot study is available
here [PlexedPiperTestData](https://github.com/vladpetyuk/PlexedPiperTestData)

Note, if decided to run using docker, you need a few system libraries 
to make package run. An example `Dockerfile` looks like this:
```{yaml}
FROM rocker/tidyverse:3.6.1
RUN apt-get update \
 && apt-get install -y --no-install-recommends \
 unixodbc \
 unixodbc-dev \
 freetds-dev \
 freetds-bin \
 tdsodbc \
 libcurl4-openssl-dev \
 libxml2-dev \
 libnetcdf-dev \
 libssl-dev
RUN R -e 'remotes::install_github("vladpetyuk/PlexedPiper", build_vignettes = TRUE")'
```

