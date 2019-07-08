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
remotes::install_github("vladpetyuk/PlexedPiper")
library(PlexedPiper)
vignette("tmt_pipeline_v1")
```

Compantion R package with test data based on MoTrPAC pilot study is available
here [PlexedPiperTestData](https://github.com/vladpetyuk/PlexedPiperTestData)