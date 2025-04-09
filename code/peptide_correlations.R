
library(tidyverse)
library(future.apply)

source('code/utility_functions.R')
source('code/plotting_functions.R')

##### Load data #####

# Protein-centric MS data
ms_data <- read_MS_data()

# Olink data
olink_data <- read_OlinkExplore_data(withControlSamples = FALSE,
                                     withReplicateAssays = FALSE,
                                     withBelowLOD = TRUE,
                                     withQCWarning = TRUE,
                                     removeAllLOD = TRUE)

olink_data_wide <- olink_long_to_wide(olink_data, 'UniProt')

# Protein annotations
protein_metadata <- read.csv('data/metadata/protein_metadata.csv') |> 
  filter(UniProt %in% ms_data$UniProt | OlinkID %in% unique(olink_data$OlinkID)) |> 
  mutate(Gene.Name = ifelse(is.na(Gene.Name), Assay, Gene.Name))

protein_metadata_MS <- read.csv('data/metadata/protein_metadata_MS.csv')

# Peptide-centric data
peptides <- read.delim('data/processed_data/peptide_table_noReplicates.txt')
peptide_metadata <- read.csv('data/metadata/peptide_metadata_MS.csv', row.names = 1)

path <- 'results/peptide_correlations'

if (!dir.exists(path)) {
  dir.create(path, recursive = TRUE)
}

# Number of peptides, proteins, and genes
nrow(peptide_metadata) # Peptides
peptide_metadata %>% 
  separate_longer_delim(UniProt.ID, ';') %>%
  pull(UniProt.ID) %>%
  unique() %>%
  length() # Proteins
peptide_metadata %>% 
  separate_longer_delim(Gene.Name, ';') %>%
  pull(Gene.Name) %>%
  unique() %>%
  length() # Genes

# 30,569 peptides from 3,240 proteins and 3,193 genes



##### All peptide correlations #####

## Filter out peptides mapping to multiple gene names
peptide_metadata_filt1 <- peptide_metadata |> 
  filter(!str_detect(Gene.Name, ';'))

peptides_filt1 <- peptides |>
  filter(Peptide.label %in% peptide_metadata_filt1$Peptide.label)

# Check that all Protein.IDs are for the same UniProt entry
peptide_metadata_filt1 |> 
  rowwise() |> 
  mutate(N.UniProt = length(unique(str_remove(unlist(str_split(UniProt.ID, ';')), '-.+')))) |> 
  pull(N.UniProt) |> 
  table() # All values should be 1

## Filter out peptides from proteins that are not found in the protein table
rows_to_keep1 <- apply(peptide_metadata_filt1, 1, function(x) {
  # Check whether at least one protein ID is found in the protein table
  any(unlist(str_split(x['UniProt.ID'], ';')) %in% ms_data$UniProt)
})
  
peptide_metadata_filt2 <- peptide_metadata_filt1[rows_to_keep1,]
peptides_filt2 <- peptides_filt1[rows_to_keep1,]

# Number of peptides, proteins, and genes
nrow(peptide_metadata_filt2) # Peptides
peptide_metadata_filt2 %>% 
  separate_longer_delim(UniProt.ID, ';') %>%
  pull(UniProt.ID) %>%
  unique() %>%
  length() # Proteins
peptide_metadata_filt2 %>% 
  pull(Gene.Name) %>%
  unique() %>%
  length() # Genes

# 28,785 peptides from 2,583 proteins and 2,543 genes

## Check how many Gene names have peptides uniquely mapping to several isoforms
isoform_genes <- peptide_metadata_filt2 |> 
  group_by(Gene.Name) |>
  filter(!str_detect(UniProt.ID, ';')) |> # unique peptides
  summarise(N.isoforms = n_distinct(UniProt.ID)) |> 
  filter(N.isoforms > 1) |> 
  pull(Gene.Name)

length(isoform_genes) # 32 genes with multiple isoforms in MS data

# Add to peptide annotation table
peptide_metadata_filt2 <- peptide_metadata_filt2 |> 
  mutate(MS.isoforms = Gene.Name %in% isoform_genes)

# Order columns so they match between peptides and Olink data
peptides_filt2 <- peptides_filt2 %>% 
  dplyr::select(Peptide.label, all_of(colnames(olink_data_wide[-1])))

# Calculate correlations
# Parallelize with chunking

num_chunks <- 100
chunk_size <- ceiling(nrow(peptides_filt2) / num_chunks)
chunks <- split(peptides_filt2, rep(1:num_chunks, each = chunk_size, length.out = nrow(peptides_filt2)))


plan(multisession)
t1 <- Sys.time()
peptide_cors <- future_lapply(chunks, function(chunk) {
  as.data.frame(cor(t(chunk[-1]), t(olink_data_wide[-1]),
                    method = 'spearman', use = 'pairwise.complete.obs'))
})
t2 <- Sys.time()
t2 - t1

peptide_cors <- bind_rows(peptide_cors)

# Change missing Peptide label (unknown protein) to UniProt ID
rownames(peptide_cors) <- peptides_filt2$Peptide.label
colnames(peptide_cors) <- olink_data_wide$UniProt

write.csv(peptide_cors, file.path(path, 'peptide-olink_cor_matrix.csv'))


# Number of samples used for calculating correlations
N_df <- as.data.frame(crossprod(!is.na(t(peptides_filt2[-1])), !is.na(t(olink_data_wide[-1]))))
rownames(N_df) <- peptides_filt2$Peptide.label
colnames(N_df) <- olink_data_wide$UniProt

write.csv(N_df, file.path(path, 'peptide-olink_cors_N_matrix.csv'))


## Correlations in long format

peptide_cors_long <- peptide_cors %>% 
  rownames_to_column('Peptide.label') %>% 
  # Add MS protein IDs
  left_join(peptide_metadata_filt2, by = 'Peptide.label') |> 
  # Long format
  pivot_longer(-all_of(colnames(peptide_metadata_filt2)), 
               names_to = 'UniProt.Olink', 
               values_to = 'Correlation') |> 
  mutate(UniProt.Olink = str_replace(UniProt.Olink, '\\.', '-')) |>
  # Rename columns
  dplyr::rename(UniProt.MS = UniProt.ID,
                Gene.Name.MS = Gene.Name) |> 
  # Add Olink assay name and gene name
  left_join(protein_metadata[c('UniProt', 'Assay', 'Gene.Name', 'OlinkID')], by = c('UniProt.Olink' = 'UniProt')) |>
  dplyr::rename(Gene.Name.Olink = Gene.Name) |>
  # Add N observations
  left_join(N_df |> 
              rownames_to_column('Peptide.label') |> 
              pivot_longer(-Peptide.label,
                           names_to = 'UniProt.Olink', 
                           values_to = 'N') |> 
              mutate(UniProt.Olink = str_replace(UniProt.Olink, '\\.', '-')),
            by = join_by(Peptide.label, UniProt.Olink))

write.csv(peptide_cors_long, file.path(path, 'peptide-olink_cors_long.csv'), row.names = FALSE)


## Correlations for matching proteins

# Filter to correlations for matching peptides/Olink proteins
# match by gene name to include correlations for different isoforms in the Shiny app
peptide_cors_matches <- peptide_cors_long |> 
  filter(Gene.Name.MS == Gene.Name.Olink) # matching proteins

write.csv(peptide_cors_matches, file.path(path, 'peptide_cors_overlappingProteins.csv'))

## Filter peptides for input to Shiny app

# Filter to peptides based on at least 15 observations
peptide_cors_matches_filt <- peptide_cors_matches |> 
  filter(N >= 15)

# Number of peptides per gene
n_peptides <- peptide_cors_matches_filt |> 
  count(Gene.Name.MS)

peptide_cors_matches_filt <- peptide_cors_matches_filt |> 
  # Filter out peptides from genes with only one peptide
  filter(Gene.Name.MS %in% n_peptides$Gene.Name.MS[n_peptides$n > 1]) |> 
  # Calculate min, median, and max number of PSMs per peptide
  rowwise() |> 
  mutate(min.PSMs = min(c_across(contains('PSM.count')), na.rm = TRUE),
         median.PSMs = median(c_across(contains('PSM.count')), na.rm = TRUE),
         max.PSMs = max(c_across(contains('PSM.count')), na.rm = TRUE)) |>
  ungroup() |> 
  # reorder columns
  dplyr::select(Peptide.label, Peptide.sequence, UniProt.MS, Gene.Name.MS, 
                Description, Correlation, N, UniProt.Olink, Gene.Name.Olink, OlinkID, Assay,
                everything(), -contains('set')) |> 
  dplyr::rename(MS.Total.number.proteins.in.group = Total.number.proteins.in.group,
                MS.Coverage = Coverage)

# Save peptide-olink correlation table
# used as input for Shiny app
write.csv(peptide_cors_matches_filt, file.path(path, 'peptide_cors_overlappingProteins_filt.csv'))

nrow(peptide_cors_matches_filt) # 13,856 peptides
length(unique(unlist(str_split(peptide_cors_matches_filt$UniProt.MS, ';')))) # 847 MS proteins
length(unique(peptide_cors_matches_filt$Gene.Name.MS)) # 822 genes
length(unique(peptide_cors_matches_filt$Assay)) # 822 Olink assays/proteins
