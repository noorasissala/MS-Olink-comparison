

library(tidyverse)
library(OlinkAnalyze)

#### MS sample metadata ####

metadata <- read.delim('data/metadata/sample_metadata.txt')

# Randomly pick one of each duplicate sample to keep
# With seed 1, this keeps  LCP27a,  LCP106b, LCP65a, LCP76a, LCP109b, LCP90a
duplicates <- metadata$Sample.ID[str_detect(metadata$Sample.ID, 'a|b')]
set.seed(1)
keep <- c()
for (i in unique(str_remove(duplicates, 'a|b'))) {
  keep <- c(keep, sample(duplicates[str_detect(duplicates, i)], size = 1))
}
remove <- duplicates[!duplicates %in% keep]

# Metadata without replicates and internal standards
metadata_noReplicates <- metadata %>% 
  filter(!Sample.ID %in% remove,
         !str_detect(Sample.ID, 'IS')) |> 
  mutate(Sample.ID = str_remove(Sample.ID, 'a|b'))

write.table(metadata_noReplicates, 'data/metadata/sample_metadata_noReplicates.txt', 
            sep = '\t',
            row.names = FALSE)


#### MS data ####

##### Protein-centric #####

proteins <- read.delim('data/raw_data/proteins_table.txt')

# Extract UniProt IDs
proteins$UniProt <- str_match(proteins$Protein.ID, '\\|(.*?)\\|')[,2]

# Extract protein metadata
protein_metadata_MS <- proteins |> 
  dplyr::select(UniProt, 1:7, contains('Fully.quanted.PSM.count'), 
                contains('Peptide.count', ignore.case = F),
                contains('peptide.count', ignore.case = F)) |> 
  dplyr::select(-Gene.ID, -Protein.ID)

# Select quant columns
proteins <- dplyr::select(proteins, UniProt, starts_with(c('S', 'IS'), ignore.case = FALSE))

# Remove proteins with NA in all channels
proteins <- filter(proteins, !if_all(-UniProt, is.na))
protein_metadata_MS <- filter(protein_metadata_MS, UniProt %in% proteins$UniProt)

write.csv(protein_metadata_MS, 'data/metadata/protein_metadata_MS.csv')
write.table(proteins, 'data/processed_data/protein_table_all.txt',
            row.names = F, sep = '\t')

# Remove duplicates and internal standards
proteins_noReplicates <- proteins %>% 
  dplyr::select(-all_of(remove), -contains('IS')) |> 
  rename_with(.fn = function(x) str_remove(x, 'a$|b$'))

write.table(proteins_noReplicates, 'data/processed_data/protein_table_noReplicates.txt',
            row.names = FALSE, sep = '\t')


##### Peptide-centric #####

peptides <- read.delim('data/raw_data/peptides_table.txt')

# Extract UniProt IDs
ids <- str_split(peptides$Protein.s., ';')
ids <- unlist(lapply(ids, function(x) paste0(str_match(x, '\\|(.*?)\\|')[,2], collapse = ';')))
peptides$Protein.s. <- ids

# Create unique readable peptide IDs
make.unique.2 = function(x, sep='.'){
  ave(x, x, FUN=function(a){if(length(a) > 1){paste(a, 1:length(a), sep=sep)} else {a}})
}
peptides$Peptide.label <- make.unique.2(peptides$Gene.name.s., sep = '_')

# Extract peptide metadata
peptide_metadata <- peptides |> 
  dplyr::select(Peptide.label, 1:7, contains('Fully.quanted.PSM.count'), -Gene.ID.s.) |> 
  dplyr::rename(Gene.Name = Gene.name.s.,
                UniProt.ID = Protein.s.,
                Description = Description.s.,
                Coverage = Coverage.s.)

# Select quant columns
peptides <- dplyr::select(peptides, Peptide.label, starts_with(c('S', 'IS'), ignore.case = FALSE))

# Remove peptides with NA in all channels
peptides <- filter(peptides, !if_all(-Peptide.label, is.na))
peptide_metadata <- filter(peptide_metadata, Peptide.label %in% peptides$Peptide.label)

write.csv(peptide_metadata, 'data/metadata/peptide_metadata_MS.csv')
write.table(peptides, 'data/processed_data/peptide_table_all.txt',
            row.names = F, sep = '\t')

# Remove duplicates and internal standards
peptides_noReplicates <- peptides %>% 
  dplyr::select(-all_of(remove), -contains('IS')) |> 
  rename_with(.fn = function(x) str_remove(x, 'a$|b$'))

write.table(peptides_noReplicates, 'data/processed_data/peptide_table_noReplicates.txt',
            row.names = FALSE, sep = '\t')


##### PSM table #####

psm_table <- read.delim('data/raw_data/target_psmtable.txt', 
                        sep = '\t') %>% 
  rename_with(.fn = function(x) str_remove(x, '\\.$')) %>%
  rename_with(.fn = function(x) paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)))) %>% 
  dplyr::rename(Master.protein = Master.protein.s,
         Protein.group.content = Protein.group.s..content,
         N.proteins.in.group = Total.number.of.matching.proteins.in.group.s,
         Missed.cleavages = Missed_cleavage,
         TMT.set = Biological.set) %>% 
  # Add UniProt ID column
  rowwise() %>%
  mutate(UniProt = paste0(unlist(
    str_match_all(Master.protein, '\\|(.*?)\\|')[[1]][,2]), 
    collapse = ';'), .before = 1) |> 
  ungroup() |> 
  mutate(TMT.set = replace(TMT.set, TMT.set == 'setE_setErerun', 'setE'))

write.table(psm_table, 'data/processed_data/target_psmtable.txt', 
            sep = '\t',
            row.names = FALSE)


#### Olink Explore data ####

# Read NPX data
NPX_data <- read_NPX('data/raw_data/VB-3207.npx.csv') %>% 
  mutate(Panel = str_replace(Panel, '_', ' ')) |> 
  dplyr::rename(Sample.ID = SampleID)

write.table(NPX_data, 'data/processed_data/OlinkExplore_data.txt', sep = '\t')

# Randomly pick one of each replicate assay to keep
# With seed 1, this keeps OID31014, OID30225, OID30563, OID20074, OID20911, OID20153
olink_proteins <- distinct(NPX_data, OlinkID, UniProt)
duplicate_proteins <- olink_proteins[duplicated(olink_proteins$UniProt),]

duplicates <- list()

for (i in unique(duplicate_proteins$UniProt)) {
  duplicates[[i]] <- unique(olink_proteins$OlinkID[olink_proteins$UniProt == i])
}

set.seed(1)
keep <- c()
for (i in duplicates) {
  keep <- c(keep, sample(i, size = 1))
} 
remove <- unlist(duplicates)[!unlist(duplicates) %in% keep]

# Remove replicate assays
NPX_data_noReplicates <- NPX_data %>% 
  filter(!OlinkID %in% remove)

write.table(NPX_data_noReplicates, 'data/processed_data/OlinkExplore_data_noReplicateAssays.txt', sep = '\t')


#### Metadata for overlapping samples #####

metadata_overlap <- metadata_noReplicates %>% 
  filter(Sample.ID %in% unique(NPX_data$Sample.ID))

write.table(metadata_overlap, 'data/metadata/sample_metadata_overlap.txt', sep = '\t')


#### Protein metadata ####

protein_metadata <- full_join(protein_metadata_MS[1:3], 
                             distinct(NPX_data[c('Panel', 'UniProt', 'OlinkID', 'Assay')]),
                         by = 'UniProt')

write.csv(protein_metadata, 'data/metadata/protein_metadata.csv', row.names = F)


#### Human Protein Atlas data ####

# Download data
if (!file.exists('data/external_data/HPA_v24.tsv')) {
  download.file('https://v24.proteinatlas.org/search?format=tsv&download=yes',
                'data/external_data/HPA_v24.tsv')
}

hpa_data <- read_tsv('data/external_data/HPA_v24.tsv', name_repair = 'universal_quiet', 
                     show_col_types = F) %>% 
  # Clean up column names
  rename_with(function(x) {str_replace_all(string = x, pattern = '\\.\\.\\.', replacement = '\\.')}) %>% 
  rename_with(function(x) {str_replace_all(string = x, pattern = '\\.\\.', replacement = '\\.')}) %>% 
  # Rename columns
  dplyr::rename(Blood.conc.MS.pgL = Blood.concentration.Conc.blood.MS.pg.L.,
                Blood.conc.IM.pgL = Blood.concentration.Conc.blood.IM.pg.L.) %>% 
  # Format Subcellular.location column to match formatting of the others
  mutate(Subcellular.location = str_replace_all(Subcellular.location, ',', ', ')) %>% 
  # Extract which tissue(s) a protein is enriched / enhanced / group enriched in
  mutate(Tissue = gsub(':.*(?=;)|:.*(?=$)', '', RNA.tissue.specific.nTPM, perl = TRUE),
         Tissue = gsub(';', ', ', Tissue),
         Tissue = gsub(' 1', '', Tissue),
         # Column that shows which tissue a protein is enriched in
         Enriched.tissue = case_when(
           RNA.tissue.specificity == 'Tissue enriched' ~ Tissue,
           RNA.tissue.specificity != 'Tissue enriched' ~ 'No enrichment')) %>% 
  rowwise() %>% 
  # Primarily use MS concentration data
  # For proteins with no MS concentration data, use the immunoassay data
  mutate(Blood.conc.pgL = ifelse(is.na(Blood.conc.MS.pgL), Blood.conc.IM.pgL, Blood.conc.MS.pgL),
         .after = Blood.conc.MS.pgL) %>% 
  ungroup() %>% 
  # Convert concentration from pg/L to ng/mL
  mutate(Blood.conc.ngmL = Blood.conc.pgL/10^6,
         .after = Blood.conc.pgL)

# Add UniProt IDs for proteins that do not have them in the HPA data
# through additional matching to MS and Olink data by gene name
# as a result, proteins will be matched between data sets by UniProt ID primarily
# and by gene name secondarily
gene_matches1 <- inner_join(
  hpa_data, protein_metadata_MS, by = join_by(Gene == Gene.Name),
  relationship = 'many-to-many') |> 
  filter(is.na(Uniprot)) |> 
  dplyr::select(UniProt, Uniprot, Gene, Gene.description, Description)

gene_matches2 <- inner_join(
  hpa_data, protein_metadata, by = join_by(Gene == Assay),
  relationship = 'many-to-many') |>
  filter(is.na(Uniprot)) |> 
  dplyr::select(UniProt, Uniprot, Gene, Gene.description, Description)

# Add UniProt IDs to HPA data
hpa_data <- hpa_data %>% 
  mutate(Uniprot = case_when(
    is.na(Uniprot) & Gene %in% gene_matches1$Gene ~ gene_matches1$UniProt[match(Gene, gene_matches1$Gene)],
    is.na(Uniprot) & Gene %in% gene_matches2$Gene ~ gene_matches2$UniProt[match(Gene, gene_matches2$Gene)],
    TRUE ~ Uniprot))

write.table(hpa_data, 'data/processed_data/HPA_v24_clean.txt', sep = '\t')

