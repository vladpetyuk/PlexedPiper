## ----setup, echo=FALSE-----------------------------------------------------
knitr::opts_chunk$set(message=FALSE, warning=FALSE)
# knitr::opts_chunk$set(echo=T, message=F, warning=F, fig.align='center', out.width='10cm')

## ----libraries-------------------------------------------------------------
library(PlexedPiper)
library(MSnID)

## ----msms------------------------------------------------------------------
path_to_MSGF_results <- system.file(
   "extdata/global/msgf_output",
   package = "PlexedPiperTestData")
msnid <- read_msgf_data(path_to_MSGF_results)
show(msnid)

## ----correct_isotopes------------------------------------------------------
msnid <- correct_peak_selection(msnid)

## ----peptide_level_filter--------------------------------------------------
show(msnid)
msnid <- filter_msgf_data_peptide_level(msnid, 0.01)
show(msnid)

## ----refseq_2_gene---------------------------------------------------------
msnid <- remap_accessions_refseq_to_gene(msnid, 
                                         organism_name="Rattus norvegicus")

## ----fasta_refseq_2_gene---------------------------------------------------
path_to_FASTA <- system.file(
   "extdata/Rattus_norvegicus_NCBI_RefSeq_2018-04-10.fasta.gz",
   package = "PlexedPiperTestData")
temp_dir <- tempdir()
file.copy(path_to_FASTA, temp_dir)
path_to_FASTA <- file.path(temp_dir, basename(path_to_FASTA))
path_to_FASTA_gene <- remap_accessions_refseq_to_gene_fasta(
   path_to_FASTA, organism_name="Rattus norvegicus")

## ----protein_level_filter--------------------------------------------------
msnid <- compute_num_peptides_per_1000aa(msnid, path_to_FASTA_gene)
show(msnid)
msnid <- filter_msgf_data_protein_level(msnid, 0.01)
show(msnid)

## ----inference-------------------------------------------------------------
msnid <- infer_parsimonious_accessions(msnid, unique_only=TRUE)
show(msnid)

## ----remove_decoy----------------------------------------------------------
msnid <- apply_filter(msnid, "!isDecoy")

## ----read_masic------------------------------------------------------------
path_to_MASIC_results <- system.file("extdata/global/masic_output", package = "PlexedPiperTestData")
masic_data <- read_masic_data(path_to_MASIC_results, interference_score=TRUE)
head(masic_data)

## ----filter_masic----------------------------------------------------------
nrow(masic_data)
masic_data <- filter_masic_data(masic_data, 0.5, 0)
nrow(masic_data)

## ----read_fractions--------------------------------------------------------
library(readr)
fractions <- read_tsv(system.file("extdata/study_design/fractions.txt", package = "PlexedPiperTestData"))
head(fractions)

## ----read_samples----------------------------------------------------------
samples <- read_tsv(system.file("extdata/study_design/samples.txt", package = "PlexedPiperTestData"))
head(samples,10)

## ----read_references-------------------------------------------------------
references <- read_tsv(system.file("extdata/study_design/references.txt", package = "PlexedPiperTestData"))
head(references)

## ----crosstab_compillation-------------------------------------------------
aggregation_level <- c("accession")
quant_cross_tab <- create_crosstab(msnid, 
                                   masic_data, 
                                   aggregation_level, 
                                   fractions, samples, references)
dim(quant_cross_tab)
head(quant_cross_tab)

## ----cleanup, echo=FALSE---------------------------------------------------
unlink(".Rcache", recursive=TRUE)

