### PEGRA Example ###

# Written by Viviane Reber
# 10.01.2024
# Please let me know if you find any issues =)

# Load packages

library(tidyverse)
library(protti)
library(multimode)
library(stats)
library(outliers)

######### Functions needed #############

source("pegRa.R")

####### Analysis #######

DIA_raw <- read_protti("directDIA_Report_rapamycin_vs_DMSO_1.tsv ") 

# Filter based on q-value and quantity 
# (filtering on RT-differnce is optional, not included here)

DIA_clean <- DIA_raw %>%
  filter(r_condition %in% c("DMSO", "RAPA")) %>% 
  filter(pep_is_proteotypic) %>% 
  filter(!grepl(";", pg_protein_accessions)) %>% 
  filter(eg_qvalue < 1e-4) %>%
  filter(fg_quantity > 1024) %>% 
  mutate(intensity_log2 = log2(fg_quantity)) %>%
  normalise(r_file_name, intensity_log2) %>%
  mutate(uniprot_id = str_extract(pg_protein_accessions, "\\w+\\-?\\w?"))

# Fetch uniprot data of all proteins in dataset for later annotation

unis <- DIA_clean %>%
  pull(pg_protein_accessions) %>%
  strsplit(";") %>%
  unlist() %>%
  unique()

uniprot <-
  fetch_uniprot(
    unis,
    columns = c(
      "protein_name",
      "length",
      "sequence",
      "gene_names",
      "gene_primary"
    )
  ) %>%
  dplyr::rename(
    protein_sequence = sequence,
    length_protein = length,
    uniprot_id = accession
  )

# Impute only MNAR peptides

set.seed(24)
DIA_clean_imputed <- DIA_clean %>% 
  impute_vivi(
    sample = r_file_name,
    replicate = r_replicate,
    protein = pg_protein_accessions,
    grouping = eg_precursor_id,
    condition = r_condition,
    treatment_condition = "RAPA",
    reference_condition = "DMSO",
    intensity = normalised_intensity_log2,
    cutoff_MAR = 3,
    cutoff_MNAR = 1
  ) %>% 
  left_join(DIA_clean %>% distinct(pep_stripped_sequence, eg_precursor_id), 
            by = c("eg_precursor_id"), relationship = "many-to-many")

# Find peptides with a low CV to use. Since no significance testing is done on the
# peptide between conditions, this is how noisy peptides are removed from the data

low_cv_peptides <- DIA_clean_imputed %>% 
  mutate(normalised_intensity = 2 ^ imputed_intensity) %>% 
  mutate(peptide = paste(eg_precursor_id, r_condition, sep = "-")) %>% 
  qc_cvs(grouping = eg_precursor_id,
         condition = peptide,
         intensity = normalised_intensity,
         plot = FALSE) %>% 
  mutate(peptide = str_replace(peptide, "[-].+", "")) %>% 
  group_by(peptide) %>% 
  mutate(max_cv = max(median_cv)) %>% 
  ungroup() %>% 
  filter(max_cv < 25) %>% 
  pull(peptide)


DIA_clean_uniprot <- DIA_clean_imputed %>%
  filter(r_condition %in% c("RAPA", "DMSO")) %>% 
  left_join(uniprot, by = c("pg_protein_accessions" = "uniprot_id")) %>%
  find_peptide(protein_sequence, pep_stripped_sequence) %>%
  assign_peptide_type(aa_before, last_aa)  %>% 
  filter(eg_precursor_id %in% low_cv_peptides) %>% 
  group_by(eg_precursor_id, r_condition) %>% 
  filter(n() >= 3) %>% 
  ungroup() %>% 
  group_by(eg_precursor_id) %>% 
  filter(n() >= 6) %>% 
  ungroup()

# Run PEGRA 

PEGRA_result <- DIA_clean_uniprot %>% 
  group_peptides_by_similarity_to_median_differential_abundance_in_largest_cluster(
    protein = pg_protein_accessions,
    peptide = eg_precursor_id,
    peptide_intensity = imputed_intensity,
    condition = r_condition,
    treatment_condition = "RAPA",
    reference_condition = "DMSO"
  ) %>%
  mutate(pep_corr_to_prot = ifelse(predicted_type_0_peptide == TRUE, 1, 0)) %>% 
  calculate_peptide_profile_significance(
    protein = pg_protein_accessions,
    peptide = eg_precursor_id,
    peptide_intensity = imputed_intensity,
    pep_corr_to_prot = pep_corr_to_prot,
    condition = r_condition,
    treatment_condition = "RAPA",
    reference_condition = "DMSO",
    structural_group = eg_precursor_id
  ) 

# Make a list of all significant proteins and save

significant_proteins <- PEGRA_result %>% 
  filter(significant == TRUE) %>% 
  filter(largest_cluster_selected_randomly == FALSE) %>% 
  pull(pg_protein_accessions) %>% 
  unique()


significant_proteins_uniprot <- uniprot %>% 
  filter(uniprot_id %in% significant_proteins)

significant_proteins_uniprot %>% write.csv("significant_proteins_uniprot.csv")

# Load biogrid interactors to help interpret structural changes

biogrid_interactors <- read_delim("BIOGRID-GENE-108757-4.4.229.tab3.txt", 
                                                                      delim = "\t") %>% 
  janitor::clean_names() %>% 
  mutate(official_symbol_interactor_a = ifelse(official_symbol_interactor_a == "MTOR", official_symbol_interactor_b, official_symbol_interactor_a)) %>% 
  pull(official_symbol_interactor_a) %>% 
  unique()

# Levels for column "new_sample_id", check column "new_sample_id" after imputation
# Can be replaces with any unique sample name, just change it in the peptide profile
# plot function as well!

levels = c("RAPA_1", "RAPA_2", "RAPA_3", "RAPA_4", "DMSO_1", "DMSO_2", "DMSO_3", "DMSO_4") 


# Plot the peptide profiles of all significant proteins (here only 10 sampled proteins
# to save time)

colorset = c('Significant'='#d9363c','Not significant'='#5680C1')

PEGRA_result_peptide_profile_plots <- PEGRA_result %>% 
  mutate(significant = ifelse(is.na(significant), 0, significant)) %>%
  mutate(significant = ifelse(significant ==1, "Significant", "Not significant")) %>% 
  mutate(significant = factor(significant, levels = c("Significant", "Not significant"))) %>% 
  mutate(new_sample_id = factor(new_sample_id, levels = levels)) %>% 
  peptide_profile_plot_2(
    sample = new_sample_id,
    peptide = eg_precursor_id,
    intensity_log2 = imputed_intensity,
    pep_corr_to_prot = significant,
    grouping = gene_primary,
    targets = significant_proteins_uniprot %>% filter(gene_primary %in% biogrid_interactors) %>% pull(gene_primary), 
    export = TRUE,
    export_name = "biogrid_interactors_significant_proteins_"
  )

