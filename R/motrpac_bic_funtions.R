#' Format Tables for BIC
#'
#' @description
#' Assembles data in format compliant with BIC requirements. 
#'
#' @description
#' * `make_rii_peptide_gl()`: returns 'RII_peptide.txt' table formatted for BIC (global)
#' * `make_results_ratio_gl()`: returns 'results_ratio.txt' table (global)
#' * `make_rii_peptide_ph()`: returns 'RII_peptide.txt' table (phospho)
#' * `make_results_ratio_ph()`: returns 'results_ratio.txt' table (phospho)
#' * `map_flanking_sequence()`: returns MSnID object with +/- 7 amino acids (by default) surrounding each PTM
#' * `assess_redundant_proteins()`: appends proteins matched to multiple peptides
#' * `assess_noninferable_proteins()`: appends proteins with identical peptide sets
#' * `compute_protein_coverage()`: determines what percent of protein is covered by given peptides.
#' 
#' @md
#'
#' @param msnid (MSnID-object) final filtered version of MSnID object
#' @param masic_data (data.frame) final filtered version of MASIC table
#' @param fractions (data.frame) Study design table linking Dataset with PlexID
#' @param samples (data.frame) Study design table linking sample names with TMT channels and PlexID
#' @param references (data.frame) Study design table describing reference value calculation
#' @param org_name (character) Organism name. Default is 'Rattus norvegicus'
#' @param sep (character) Single character used to concatenate protein, SiteID, and peptide
#'
#' @importFrom dplyr select inner_join left_join mutate %>% case_when rename group_by summarize
#' @importFrom tibble rownames_to_column
#' @importFrom purrr map map2
#' @importFrom IRanges IRanges IRangesList reduce
#'
#'
#' @name motrpac_bic_output
#'
#' @examples
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
#' library(readr)
#' fractions <- read_tsv(system.file("extdata/study_design/fractions.txt", package = "PlexedPiperTestData"))
#' samples <- read_tsv(system.file("extdata/study_design/samples.txt", package = "PlexedPiperTestData"))
#' references <- read_tsv(system.file("extdata/study_design/references.txt", package = "PlexedPiperTestData"))
#' 
#' results_ratio <- make_results_ratio_gl(msnid, masic_data, fractions, samples, references, org_name = "Rattus norvegicus")
#' rii_peptide <- make_rii_peptide_gl(msnid, masic_data, fractions, samples, references, org_name = "Rattus norvegicus")


#' @export
#' @rdname motrpac_bic_output
make_rii_peptide_gl <- function(msnid, masic_data, fractions, samples, 
                                references, org_name = "Rattus norvegicus") {
  ## Make RII study design tables
  samples_rii <- samples %>%
    mutate(MeasurementName = case_when(is.na(MeasurementName) ~ paste0("Ref", "_", PlexID),
                                       TRUE ~ MeasurementName))
  references_rii <- references %>%
    mutate(Reference = 1)
  
  ## Create crosstab
  aggregation_level <- c("accession", "peptide")
  crosstab <- create_crosstab(msnid, 
                              masic_data, 
                              aggregation_level, 
                              fractions, samples_rii, references_rii)
  crosstab <- 2^crosstab # undo log2
  
  crosstab <- as.data.frame(crosstab) %>%
    rownames_to_column("Specie")
  
  ## Create RII peptide table
  rii_peptide <- crosstab %>%
    select(Specie) %>%
    mutate(protein_id = sub("(^.*\\.\\d+)@.*", "\\1", Specie),
           sequence = sub("(^.*)@(.*)", "\\2", Specie),
           organism_name = org_name) %>%
    mutate(REFSEQ = sub("(^.*)\\.\\d+", "\\1", protein_id))
  
  ## Attach Gene symbol and Entrez ID
  conv <- fetch_conversion_table(org_name, from = "REFSEQ", "SYMBOL")
  
  conv <- fetch_conversion_table(org_name, from = "REFSEQ", "ENTREZID") %>%
    inner_join(., conv)
  
  rii_peptide <- rii_peptide %>%
    inner_join(conv) %>%
    rename(gene_symbol = SYMBOL,
           entrez_id = ENTREZID) %>%
    select(-REFSEQ)
  
  ## Additional info from MS/MS
  ids <- psms(msnid) %>%
    select(accession, peptide, redundantAccessions, MSGFDB_SpecEValue) %>%
    rename(protein_id = accession,
           sequence = peptide,
           redundant_accessions = redundantAccessions) %>%
    group_by(protein_id, sequence, redundant_accessions) %>%
    summarize(peptide_score = min(MSGFDB_SpecEValue))
  
  rii_peptide <- inner_join(rii_peptide, ids)
  
  ## Join with crosstab
  rii_peptide <- inner_join(rii_peptide, crosstab) %>%
    select(-Specie)
  
  return(rii_peptide)
}


#' @export
#' @rdname motrpac_bic_output
make_results_ratio_gl <- function(msnid, masic_data, fractions, samples, 
                                  references, org_name = "Rattus norvegicus") {
  ## Create crosstab
  aggregation_level <- c("accession")
  crosstab <- create_crosstab(msnid, masic_data, aggregation_level, fractions,
                              samples, references)
  
  crosstab <- as.data.frame(crosstab) %>%
    rownames_to_column("protein_id")
  
  ## Create results ratio table
  results_ratio <- crosstab %>%
    select(protein_id) %>%
    mutate(REFSEQ = sub("(^.*)\\.\\d+", "\\1", protein_id))
  
  ## Attach Gene symbol and Entrez ID
  conv <- fetch_conversion_table(org_name, from = "REFSEQ", "SYMBOL")
  
  conv <- fetch_conversion_table(org_name, from = "REFSEQ", "ENTREZID") %>%
    inner_join(., conv)
  
  results_ratio <- results_ratio %>%
    inner_join(conv) %>%
    rename(gene_symbol = SYMBOL,
           entrez_id = ENTREZID) %>%
    select(-REFSEQ)
    
  
  ## Additional info from MS/MS
  ids <- psms(msnid) %>%
    select(accession, peptide, noninferableProteins, percentAACoverage,
           MSGFDB_SpecEValue) %>%
    rename(protein_id = accession,
           sequence = peptide,
           noninferable_proteins = noninferableProteins,
           percent_coverage = percentAACoverage) %>%
    group_by(protein_id, sequence, noninferable_proteins, percent_coverage) %>%
    summarize(peptide_score = min(MSGFDB_SpecEValue)) %>%
    group_by(protein_id, noninferable_proteins, percent_coverage) %>%
    summarize(protein_score = min(peptide_score),
              num_peptides = n()) 
  
  results_ratio <- inner_join(results_ratio, ids)
  
  ## Join with crosstab
  results_ratio <- inner_join(results_ratio, crosstab)
  
  return(results_ratio)
}


#' @export
#' @rdname motrpac_bic_output
make_rii_peptide_ph <- function(msnid, masic_data, fractions, samples, references,
                                org_name = "Rattus norvegicus", sep="_") {
  ## Make RII study design tables
  samples_rii <- samples %>%
    mutate(MeasurementName = case_when(is.na(MeasurementName) ~ paste0("Ref", "_", PlexID),
                                       TRUE ~ MeasurementName))
  
  
  references_rii <- references %>%
    mutate(Reference = 1)
  
  ## Create crosstab
  aggregation_level <- c("accession", "peptide", "SiteID")
  crosstab <- create_crosstab(msnid, 
                              masic_data, 
                              aggregation_level, 
                              fractions, samples_rii, references_rii)
  
  crosstab <- 2^crosstab # undo log2
  
  crosstab <- as.data.frame(crosstab) %>%
    rownames_to_column("Specie")
  
  ## Create RII peptide table
  rii_peptide <- crosstab %>%
    select(Specie) %>%
    mutate(protein_id = sub("(^.*)@(.*)@(.*)", "\\1", Specie),
           sequence = sub("(^.*)@(.*)@(.*)", "\\2", Specie),
           ptm_id = sub("(^.*)@(.*)@(.*)", "\\3", Specie),
           organism_name = org_name) %>%
    mutate(REFSEQ = sub("(^.*)\\.\\d+", "\\1", protein_id))
  
  
  ## Add Genes + EntrezID
  conv <- fetch_conversion_table(org_name, from = "REFSEQ", "SYMBOL")
  
  conv <- fetch_conversion_table(org_name, from = "REFSEQ", "ENTREZID") %>%
    inner_join(., conv)
  
  rii_peptide <- inner_join(rii_peptide, conv) %>%
    rename(gene_symbol = SYMBOL,
           entrez_id = ENTREZID) %>%
    select(-REFSEQ)
  
  
  ## Additional info from MS/MS
  ids <- psms(msnid) %>%
    select(accession, peptide, SiteID,
           flankingSequence, redundantAccessions, MSGFDB_SpecEValue, maxAScore) %>%
    rename(protein_id = accession,
           sequence = peptide,
           ptm_id = SiteID,
           flanking_sequence = flankingSequence,
           redundant_accessions = redundantAccessions) %>%
    group_by(protein_id, sequence, ptm_id, flanking_sequence, redundant_accessions) %>%
    summarize(peptide_score = min(MSGFDB_SpecEValue),
              confident_score = max(maxAScore)) %>%
    mutate(confident_site = case_when(confident_score >= 17 ~ TRUE,
                                      confident_score < 17 ~ FALSE))
  
  rii_peptide <- inner_join(rii_peptide, ids) %>%
    mutate(ptm_id = gsub("-", sep, ptm_id))
  
  ## Join with crosstab
  rii_peptide <- inner_join(rii_peptide, crosstab) %>%
    select(-Specie)
  
  return(rii_peptide)
}

#' @export
#' @rdname motrpac_bic_output
make_results_ratio_ph <- function(msnid, masic_data, fractions, samples, 
                                  references, org_name = "Rattus norvegicus", sep="_") {
  
  aggregation_level <- c("accession", "SiteID")
  crosstab <- create_crosstab(msnid, masic_data, aggregation_level, fractions,
                              samples, references)
  crosstab <- as.data.frame(crosstab) %>% 
    rownames_to_column('Specie')
  
  ## Create RII peptide table
  results_ratio <- crosstab %>%
    select(Specie) %>%
    mutate(protein_id = sub("(^.*)@(.*)", "\\1", Specie),
           ptm_id = sub("(^.*)@(.*)", "\\2", Specie),
           organism_name = org_name) %>%
    mutate(REFSEQ = sub("(^.*)\\.\\d+", "\\1", protein_id))
  
  ## Add Genes + EntrezID
  conv <- fetch_conversion_table(org_name, from = "REFSEQ", "SYMBOL")
  
  conv <- fetch_conversion_table(org_name, from = "REFSEQ", "ENTREZID") %>%
    inner_join(., conv)
  
  results_ratio <- inner_join(results_ratio, conv) %>%
    rename(gene_symbol = SYMBOL,
           entrez_id = ENTREZID) %>%
    select(-REFSEQ)
  
  ## Additional info from MS/MS
  ids <- psms(msnid) %>%
    select(accession, peptide, SiteID,
           flankingSequence, MSGFDB_SpecEValue, maxAScore) %>%
    rename(protein_id = accession,
           sequence = peptide,
           ptm_id = SiteID,
           flanking_sequence = flankingSequence) %>%
    # group at peptide level to calculate peptide score, confident score
    group_by(protein_id, sequence, ptm_id, flanking_sequence) %>%
    summarize(peptide_score = min(MSGFDB_SpecEValue),
              confident_score = max(maxAScore)) %>%
    # regroup at siteID level and recalculate ptm score
    group_by(protein_id, ptm_id, flanking_sequence) %>%
    summarize(ptm_score = min(peptide_score),
              confident_score = max(confident_score)) %>%
    mutate(confident_site = case_when(confident_score >= 17 ~ TRUE,
                                      confident_score < 17 ~ FALSE))
  
  results_ratio <- inner_join(results_ratio, ids) %>%
    mutate(ptm_id = gsub("-", sep, ptm_id))
  
  ## Join with crosstab
  results_ratio <- inner_join(results_ratio, crosstab) %>%
    select(-Specie)
  
  return(results_ratio)
}



#' @export
#' @rdname motrpac_bic_output
map_flanking_sequence <- function (msnid, fasta, radius=7L, collapse="|") {
  
  # This function takes every phosphosite
  # and appends the amino acids within +/- 7 neighborhood
  # along the corresponding protein.
  # If there are multiple phosphosites on a single peptide then
  # they are pasted together.
  
  if (!("SiteID" %in% names(msnid))) {
    stop("No SiteID found. Call map_mod_sites.")
  }
  
  x <- psms(msnid) %>%
    dplyr::select(accession, Peptide, SiteLoc) %>%
    distinct()
  
  x <- fasta %>%
    as.data.frame() %>%
    rownames_to_column("accession") %>%
    dplyr::mutate(accession = sub("^(.P_\\d+\\.\\d+)?\\s.*", "\\1", accession)) %>%
    dplyr::rename(ProtSeq = x) %>%
    dplyr::mutate(ProtSeqWidth = nchar(ProtSeq)) %>%
    inner_join(x, .)
  
  
  f <- function(ProtSeq_i, SiteLoc_i) {
    flankingSequences <- c()
    for (k in unlist(SiteLoc_i)) {
      
      site_left <- substr(ProtSeq_i, max(k-radius,1), k-1)
      
      if (k-radius < 1) {
        site_left <- paste0(paste(rep("-", 1-(k-radius)),collapse=""), site_left)
      }
      
      site_right <- substr(ProtSeq_i, k+1, min(k+radius,nchar(ProtSeq_i)))
      
      if (k+radius > nchar(ProtSeq_i)) {
        site_right <- site_right <- paste0(site_right, paste(rep("-", k+radius-nchar(ProtSeq_i)),collapse=""))
      }
      
      mod_aa <- tolower(substr(ProtSeq_i, k, k))
      
      flank <- paste0(site_left, tolower(mod_aa), site_right)
      
      flankingSequences <- c(flankingSequences, flank)
    }
    flankingSequences <- paste(flankingSequences, collapse=collapse)
    return(flankingSequences)
  }
  
  
  x$flankingSequence <-  map2(x$ProtSeq, x$SiteLoc, f)
  x$flankingSequence <- as.character(x$flankingSequence)
  
  msnid@psms <- psms(msnid) %>% mutate(flankingSequence=NULL) %>%
    left_join(x) %>% data.table()
  
  return(msnid)
}



#' @export
#' @rdname motrpac_bic_output
assess_redundant_protein_matches <- function(msnid, collapse="|") {
  msnid@psms <- psms(msnid) %>%
    dplyr::select(accession, peptide) %>%
    distinct() %>%
    group_by(peptide) %>%
    summarize(redundantAccessions = paste(accession, collapse=collapse)) %>%
    left_join(psms(msnid), .) %>%
    data.table()
  return(msnid)
}

#' @export
#' @rdname motrpac_bic_output
assess_noninferable_proteins <- function(msnid, collapse="|") {
  
  # assign each accession to its peptide set
  x <- psms(msnid) %>%
    dplyr::select(accession, Peptide) %>%
    group_by(accession) %>%
    arrange(Peptide) %>%
    summarize(peptide_set = paste(Peptide, collapse=collapse))
  
  # group together accessions with identical peptide sets
  x <- x %>%
    group_by(peptide_set) %>%
    summarize(noninferableProteins = paste(accession, collapse=collapse)) %>%
    left_join(x) %>%
    dplyr::select(-peptide_set)
    
  msnid@psms <- psms(msnid) %>% mutate(noninferableProteins=NULL) %>%
    left_join(x) %>% data.table()
  
  return(msnid)
}


#' @export
#' @rdname motrpac_bic_output
compute_protein_coverage <- function(msnid, path_to_FASTA) {
  
  fasta <- readAAStringSet(path_to_FASTA)
  names(fasta) <- sub("^(\\S+)\\s.*", "\\1", names(fasta))
                                             
  if(any(duplicated(names(fasta)))) {
        stop("FASTA entry names are not unique!")
  }
  
  if(length(intersect(msnid$Protein, names(fasta))) == 0) {
    stop("There is zero overlap in protein IDs and FASTA entry names!")
  }
  
  get_coverage_for_single_protein <- function(protein_i, ids, fasta) {
    
    
    
    x <- ids %>% 
      filter(Protein == protein_i) %>%
      distinct()
    
    protAAstring <- fasta[[protein_i]]
    
    # main loop
    irl <- IRangesList()
    for(i in 1:nrow(x)) {
      mtch <- regexpr(x$pepSeq[i], protAAstring)
      start <- as.numeric(mtch)
      width <- attr(mtch, "match.length")
      tgt <- IRanges(start=start, width=width, names=x$pepSeq[i])
      irl[[i]] <- tgt
    }
    
    fullAACoverage <- sum(width(reduce(unlist(irl))))
    percentAACoverage <- 100*fullAACoverage/length(protAAstring)
    return(percentAACoverage)
  }
  
  ids <- psms(msnid) %>% dplyr::select(Protein, pepSeq) %>%
    distinct()
  
  proteins <- ids %>% dplyr::select(Protein) %>% distinct()
  
  proteins$percentAACoverage <-  map(proteins$Protein, get_coverage_for_single_protein, ids, fasta)
  proteins$percentAACoverage <- as.numeric(proteins$percentAACoverage)
  
  msnid@psms <- psms(msnid) %>% mutate(percentAACoverage=NULL) %>%
    left_join(proteins) %>% data.table()
  
  return(msnid)
}




