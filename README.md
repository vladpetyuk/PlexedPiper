# PlexedPiper

R package used at PNNL for processing isobaric labeling (e.g. TMT) proteomics
data. The key inputs are:
* MS/MS identifications from the [MS-GF+ search engine](https://github.com/MSGFPlus/msgfplus)
* Reporter ion intensities extracted using [MASIC](https://github.com/pnnl-comp-mass-spec/MASIC)
* Tables outlying study design
   * table linking dataset to plexes
   * table linking reporter channels with sample names
   * table identifying reference within each plex

## R Installation and Usage

```{r}
if(!require("remotes", quietly = T)) install.packages("remotes")
remotes::install_github("vladpetyuk/PlexedPiper", build_vignettes = TRUE)
library(PlexedPiper)
vignette("tmt_pipeline_v1")
```

### Example Data

A companion R package with test data based on the MoTrPAC pilot study is available
here
* [PlexedPiperTestData](https://github.com/vladpetyuk/PlexedPiperTestData)

## Docker/Linux installation

PlexedPiper can be run within a [Docker Container](https://www.docker.com/resources/what-container)
* This example `Dockerfile` shows the required system libraries, starting with the base [rocker/TidyVerse](https://hub.docker.com/r/rocker/tidyverse/dockerfile) image

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
RUN R -e 'remotes::install_github("vladpetyuk/PlexedPiper", build_vignettes = TRUE)'
```

## MacOS installation

On MacOS, install [Homebrew](https://brew.sh/), then use

```{shell}
brew install unixodbc
brew install freetds
```
Note, the `--with-unixodbc` option in freetds installation is deprecated.

Create `~/.odbcinst.ini` file and add
```{yaml}
[FreeTDS]
Driver = /usr/local/lib/libtdsodbc.so
```
If your location of `libtdsodbc.so` differs, use the proper location.
