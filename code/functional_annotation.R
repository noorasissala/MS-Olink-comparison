

##### Setup #####

###### Load functions and packages ######
 
source('code/utility_functions.R')
source('code/plotting_functions.R')

library(tidyverse)
library(patchwork)
library(scales)
library(ggpp)
library(ggvenn)
library(eulerr)
library(rlang)
library(clusterProfiler)
library(org.Hs.eg.db)


###### Read in data ######

# Output directory
path <- 'results/functional_annotation'

if (!dir.exists(path)) {
  dir.create(path, recursive = TRUE)
}

# Read in MS data
ms_data <- read_MS_data()

# Olink Explore data
olink_data <- read_OlinkExplore_data()

# Olink Explore data in wide format
olink_data_wide <- olink_long_to_wide(olink_data, id_cols = 'UniProt')

# Metadata
metadata <- read_metadata()

# Protein annotations
protein_metadata <- read.csv('data/metadata/protein_metadata.csv') |> 
  filter(UniProt %in% ms_data$UniProt | OlinkID %in% olink_data$OlinkID)

# Vectors of protein IDs
ms_proteins <- unique(ms_data$UniProt)
olink_proteins <- unique(olink_data$UniProt)
all_proteins <- union(ms_proteins, olink_proteins)
ms_unique <- setdiff(ms_proteins, olink_proteins)
olink_unique <- setdiff(olink_proteins, ms_proteins)

# Load HPA data
hpa_data <- read.delim('data/processed_data/HPA_v24_clean.txt')

# HPA data for proteins detected with MS
hpa_ms <- filter(hpa_data, Uniprot %in% ms_proteins) |> 
  left_join(protein_metadata, by = join_by(Uniprot == UniProt))  # match by UniProt ID
  
# Identify duplicates and do additional matching based on gene name
ms_dups <- unique(hpa_ms$Uniprot[duplicated(hpa_ms$Uniprot)])

# Apply filtering to HPA data
hpa_ms <- hpa_ms %>% 
  filter(!Uniprot %in% ms_dups |
           (Uniprot %in% ms_dups & Gene.Name == Gene)) |> 
  distinct(Uniprot, Gene.Name, .keep_all = TRUE) # if there are still duplicates, keep the first one

# Check that no matches were lost entirely through additional filtering
stopifnot(ms_dups %in% unique(hpa_ms$Uniprot))
sum(duplicated(hpa_ms$Uniprot))

# HPA data for proteins detected with Olink Explore
hpa_olink <- filter(hpa_data, Uniprot %in% olink_proteins) |> 
  left_join(protein_metadata, by = join_by(Uniprot == UniProt))  # match by UniProt ID

# Identify duplicates and do additional matching based on gene name
olink_dups <- unique(hpa_olink$Uniprot[duplicated(hpa_olink$Uniprot)])

# Apply filtering to HPA data
hpa_olink <- hpa_olink %>% 
  mutate(Assay2 = Assay) |>  # keep original assay name
  separate_longer_delim(Assay2, '_') |> 
  filter(!Uniprot %in% olink_dups |
           (Uniprot %in% olink_dups & Gene == Assay2)) |> 
  distinct(Uniprot, Assay, .keep_all = TRUE) |>  # if there are still duplicates, keep the first one
  dplyr::select(-Assay2)  # remove temporary column

# Check that no matches were lost entirely through additional filtering
stopifnot(olink_dups %in% unique(hpa_olink$Uniprot))

# All detected proteins found in HPA
hpa_detected <- bind_rows(hpa_ms, hpa_olink) |> 
  distinct()

# Number and proportion of MS / Olink proteins not found in HPA
sum(!ms_proteins %in% hpa_data$Uniprot)
sum(!ms_proteins %in% hpa_data$Uniprot)/length(ms_proteins)*100
sum(!olink_proteins %in% hpa_data$Uniprot)
sum(!olink_proteins %in% hpa_data$Uniprot)/length(olink_proteins)*100


#### Functions ####

# Calculate frequencies for annotations:

calc_freq <- function(df, col) {
  
  freq <- table(unlist(str_split(df[[col]], ', ')), useNA = 'always')
  names(freq)[is.na(names(freq))] <- 'Unknown'
  names(freq) <- str_to_sentence(names(freq))
  
  freq <- freq %>% 
    enframe(name = col, value = 'N') %>% 
    mutate(N = as.numeric(N),
           Freq = N/nrow(df)*100)
  
  freq
  
}


# Test for difference in proportions:

fisher_pvalue <- function(a, b, c, d) {
  data <- matrix(c(a, b, c, d), ncol = 2, nrow = 2, byrow = T)
  
  if (any(is.na(data))) {
    return(NA_real_)
  }
  else {
    return(fisher.test(data)$p.value)
  }
}


# Calculate enrichment (overrepresentation analysis):

hyper.test <- function(overlap, protein_set, size_target, size_bg) {
  # Hypergeometric test for overrepresentation analysis
  # overlap: number of elements of protein set observed in target (i.e overlap between protein set and target)
  # protein_set: number of elements of protein set observed in background (i.e size of protein set)
  # size_target: total number of proteins in target (i.e sample size)
  # size_bg: total number of proteins in background (i.e database size)
  
  # Test for overrepresentation
  phyper(q = overlap-1, m = protein_set, n = size_bg - protein_set, k = size_target, lower.tail = FALSE)
  
}

# Merge results and perform tests for differences in proportions
# and enrichment compared to background gene list:

ann_freq_summary <- function(bg_ann, ms_ann, olink_ann, col) {
  # bg_ann: Protein annotations for desired background protein list
  # ms_ann: Annotations filtered for proteins detected with MS
  # olink_ann: Annotations filtered for proteins detected with Olink
  # col: Column in the data to calculate frequencies on
  
  # Calculate frequencies of different categories
  freq_bg <- calc_freq(bg_ann, col)
  freq_ms <- calc_freq(ms_ann, col)
  freq_olink <- calc_freq(olink_ann, col)
  
  # Test for difference in proportions and enrichment compared to background
  ann_summary <- freq_bg %>% 
    full_join(freq_ms, by = col, suffix = c('_BG', '_MS')) %>% 
    full_join(freq_olink, by = col) %>% 
    dplyr::rename(N_Olink = N,
                  Freq_Olink = Freq,
                  Annotation = all_of(col)) %>% 
    add_column(HPA.category = col, .before = 1) %>%
    # Set NA frequencies to 0
    mutate(across(N_MS:Freq_Olink, ~ifelse(is.na(.), 0, .))) %>% 
    rowwise() %>% 
    mutate('N.diff' = N_MS - N_Olink,
           'Freq.diff' = Freq_MS - Freq_Olink,
           'Fisher.pvalue' = fisher_pvalue(N_MS, N_Olink, nrow(ms_ann)-N_MS, nrow(olink_ann)-N_Olink),
           'Higher' = case_when(Freq.diff < 0 ~ 'Olink',
                                Freq.diff > 0 ~ 'MS',
                                Freq.diff == 0 ~ 'No difference'),
           'Enr.pvalue_MS' = hyper.test(N_MS, N_BG, nrow(ms_ann), nrow(bg_ann)),
           'Enr.pvalue_Olink' = hyper.test(N_Olink, N_BG, nrow(olink_ann), nrow(bg_ann))) %>% 
    ungroup() %>% 
    mutate('Fisher.adj.p' = p.adjust(Fisher.pvalue, method = 'fdr'),
           'Enr.adj.p_MS' = p.adjust(Enr.pvalue_MS, method = 'fdr'),
           'Enr.adj.p_Olink' = p.adjust(Enr.pvalue_Olink, method = 'fdr'),
           'Enriched' = case_when(Enr.adj.p_MS < 0.05 & Enr.adj.p_Olink >= 0.05 ~ 'MS',
                                  Enr.adj.p_Olink < 0.05 & Enr.adj.p_MS >= 0.05 ~ 'Olink',
                                  Enr.adj.p_Olink < 0.05 & Enr.adj.p_MS < 0.05 ~ 'Both',
                                  TRUE ~ 'None')) |> 
    dplyr::relocate(Fisher.adj.p, .after = Fisher.pvalue) |> 
    dplyr::relocate(Enr.adj.p_MS, .after = Enr.pvalue_MS) |>
    dplyr::relocate(Enr.adj.p_Olink, .after = Enr.pvalue_Olink)
    
  
  return(ann_summary)
  
}


# Frequency barplot:

freq_barplot <- function(ann_summary, show.terms = c('most frequent', 'most different', 'enriched')) {
  # ann_summary: result from ann_freq_summary
  
  show.terms <- match.arg(show.terms)
  
  col <- unique(ann_summary$HPA.category)
  
  plot_data <- ann_summary %>% 
    pivot_longer(cols = c(Freq_MS, Enr.pvalue_MS, Enr.adj.p_MS, 
                          Freq_Olink, Enr.pvalue_Olink, Enr.adj.p_Olink),
                 names_to = c('.value', 'Platform'),
                 names_sep = '_') |> 
    filter(!Annotation == 'Unknown', !Annotation == 'No enrichment')
  
  if (show.terms == 'most frequent') {
    
    plot_data <- plot_data %>% 
      group_by(Platform) %>% 
      arrange(-Freq) %>% 
      mutate(Rank = row_number()) %>%
      group_by(Annotation) %>%
      mutate(Keep = ifelse(any(Rank <= 15), TRUE, FALSE)) %>% 
      filter(Keep)
    
  }
  
  else if (show.terms == 'most different') {
    
    plot_data <- plot_data %>% 
      group_by(Higher) %>% 
      arrange(-abs(Freq.diff)) %>% 
      mutate(Rank = row_number()) %>%
      group_by(Annotation) %>%
      mutate(Keep = ifelse(any(Rank <= 15), TRUE, FALSE)) %>% 
      filter(Keep)
    
  }
  
  else if (show.terms == 'enriched') {
    
    plot_data <- plot_data %>% 
      filter(Enriched %in% c('MS', 'Olink', 'Both'))
    
  }
  
  plot_data %>% 
    mutate(Label = ifelse(Fisher.adj.p < 0.05 & Platform == Higher, '*', NA)) %>% 
    ggplot(aes(y = reorder(Annotation, Freq), x = Freq, fill = Platform)) +
    geom_col(position = position_dodge(width = 0.9, preserve = 'single')) +
    geom_text(aes(label = Label),
              position = position_dodge(width = 0.9), 
              vjust = 0.75, hjust = -0.5, size = 6/.pt) +
    labs(title = str_replace_all(col, '\\.', ' '), y = '', x = 'Frequency (%)') +
    scale_fill_manual(values = TDP_palette[c(1,4)]) +
    theme_publ()
  
}


#### HPA annotations ####

cols <- c('Chromosome', 'Protein.class', 'Biological.process', 'Molecular.function',
          'Disease.involvement', 'Evidence', 'HPA.evidence', 'UniProt.evidence',
          'NeXtProt.evidence', 'RNA.tissue.specificity', 'RNA.tissue.distribution',
          'RNA.single.cell.type.specificity', 'RNA.single.cell.type.distribution',
          'RNA.cancer.specificity', 'RNA.cancer.distribution', 
          'RNA.brain.regional.specificity', 'RNA.brain.regional.distribution',
          'RNA.blood.cell.specificity', 'RNA.blood.cell.distribution',
          'RNA.blood.lineage.specificity', 'RNA.blood.lineage.distribution',
          'Blood.expression.cluster', 'Tissue.expression.cluster',
          'Subcellular.location', 'Secretome.location', 'Subcellular.main.location',
          'Subcellular.additional.location', 'Enriched.tissue')

# Run analysis
hpa_ann_freq <- setNames(
  lapply(cols, 
         ann_freq_summary, 
         bg_ann = hpa_detected, 
         ms_ann = hpa_ms, 
         olink_ann = hpa_olink), 
  cols)

hpa_ann_freq_df <- do.call(rbind.data.frame, hpa_ann_freq)

write.csv(hpa_ann_freq_df, paste0(path, '/HPA_ann_frequencies.csv'),
          row.names = FALSE)

# Plot results
hpa_freq_plots <- setNames(lapply(hpa_ann_freq, freq_barplot), cols)

# Selected a subset of interesting annotations to show in the paper
p <- wrap_plots(hpa_freq_plots[c('Protein.class', 'Secretome.location', 'Enriched.tissue', 'Subcellular.location')]) +
  plot_layout(ncol = 2, axis_titles = 'collect')

freq_plot <- guide_area() + p +
  plot_layout(guides = 'collect', nrow = 2, heights = c(0.4,10), design = 'A\nB')
  
freq_plot

ggsave(freq_plot, filename = 'HPA_freq_barplots.pdf',
       path = path,
       width = 18.3,
       height = 12,
       units = 'cm')


##### GO terms #####
 
# Target protein lists
input_list <- list('MS' = ms_proteins,
                   'Olink' = olink_proteins)

# Background protein list
background <- all_proteins


GO.analysis <- function(input_list, background, ont = c('ALL', 'BP', 'MF', 'CC')) {
  
  ont <- match.arg(ont)
  
  # Perform GO enrichment analysis
  GO_res <- compareCluster(input_list, 
                           fun = 'enrichGO',
                           OrgDb = org.Hs.eg.db, 
                           keyType = 'UNIPROT', 
                           universe = background,
                           ont = ont)
  
  # Calculate fold enrichment and add column with indication of which 
  # platform each term is enriched in
  res_ms <- GO_res@compareClusterResult %>% filter(Cluster == 'MS')
  res_olink <- GO_res@compareClusterResult %>% filter(Cluster == 'Olink')
  overlapping_IDs <- intersect(res_ms$ID, res_olink$ID)
  
  GO_res@compareClusterResult <- GO_res@compareClusterResult %>% 
    rowwise() %>% 
    mutate(GeneRatio2 = eval(parse_expr(GeneRatio)), # calculate GeneRatio
           BgRatio2 = eval(parse_expr(BgRatio)), # calculate BgRatio
           foldEnr = GeneRatio2/BgRatio2, # calculate fold enrichment
           Enriched = case_when(ID %in% overlapping_IDs ~ 'Both',
                                !ID %in% overlapping_IDs & Cluster == 'MS' ~ 'MS',
                                !ID %in% overlapping_IDs & Cluster == 'Olink' ~ 'Olink'),
           Description = str_to_sentence(Description)) %>% # Description in sentence case
    ungroup()
  
  return(GO_res)
  
}

plot.enr.dotplot <- function(enr_res, showCategory = 10) {
  
  library(RColorBrewer)
  
  if ('ONTOLOGY' %in% colnames(enr_res@compareClusterResult)) {
      facet <- facet_wrap(vars(Cluster, ONTOLOGY), scales = 'free_y')
      format_facet <- theme(
        strip.background = element_rect(fill = 'grey90', color = 'black'),
        strip.text = element_text(hjust = 0.5))
      
      # Select terms to plot
      top <- enr_res@compareClusterResult %>% 
        group_by(Cluster, ONTOLOGY) %>% 
        slice_min(p.adjust, n = showCategory) %>% 
        pull(ID) %>% 
        unique()
      
    } else {
      facet <- facet_wrap(vars(Cluster), scales = 'free_y')
      format_facet <- NULL
      
      # Select terms to plot
      top <- enr_res@compareClusterResult %>% 
        group_by(Cluster) %>% 
        slice_min(p.adjust, n = showCategory) %>% 
        pull(ID) %>% 
        unique()
      
    }
    
    # Filter results
    dat <- enr_res@compareClusterResult %>% 
      filter(ID %in% top) %>% 
      arrange(desc(Cluster), foldEnr) %>% 
      mutate(Order = row_number()) 
    
    # Plot results
    p <- ggplot(dat, aes(x = foldEnr, y = reorder(Description, Order), 
                         fill = Count)) +
      geom_point(shape = 21, size = 1.8, color = 'grey30', stroke = 0.2) +
      scale_fill_gradientn(colors = rev(brewer.pal(5, 'Spectral'))) +
      labs(x = 'Fold enrichment', y = '') +
      theme_publ() +
      theme(legend.key.width = unit(1, 'cm')) +
      facet +
      format_facet
    
}


GO_res_BP <- GO.analysis(input_list = input_list, 
                         background = background,
                         ont = 'BP')

# Result data frame
GO_res_BP_df <- as.data.frame(GO_res_BP)
write.csv(GO_res_BP_df, file.path(path, 'GO_BP_res.csv'))

# Plot venn diagram of enriched GO terms
res_ms <- GO_res_BP@compareClusterResult %>% filter(Cluster == 'MS')
res_olink <- GO_res_BP@compareClusterResult %>% filter(Cluster == 'Olink')
ggvenn(list('MS' = res_ms$ID, 'Olink' = res_olink$ID),
             fill_color = c("#21a2ac","#f38895"),
             stroke_color = rep(c("#21a2ac","#f38895"), each = 100),
             stroke_size = 1,
             text_size = 6/.pt,
             set_name_size = 7/.pt) + 
  coord_cartesian(xlim = c(-2, 2), ylim = c(-2, 2)) +
  labs(title = paste('ORA of GO BP terms')) +
  theme(title = element_text(size = 8))

ggsave(filename = 'GO_BP_venn.pdf',
       path = path,
       width = 5.5,
       height = 6.5,
       units = 'cm')

# Plot dotplots of enriched GO terms
GO_BP_dotplot <- plot.enr.dotplot(GO_res_BP)

GO_BP_dotplot

ggsave(GO_BP_dotplot, 
       filename = 'GO_BP_dotplot.pdf',
       path = path,
       width = 18.3,
       height = 6.5,
       units = 'cm')

##### Coverage of clinically used biomarkers #####

# Load FDA biomarker data from Anderson (Clinical Chemistry, 2010)
biomarkers <- readxl::read_excel('data/external_data/126706_supplemental_table_1.xls',
                                 range = 'B1:R110', .name_repair = 'universal_quiet') |> 
  dplyr::rename(Protein.Name = Protein.name.from.FDA,
                UniProt = SwissProt.Accession.s.,
                Blood.conc.ugmL = Normal.value..Hortin..ug.ml.) |> 
  mutate(UniProt = ifelse(UniProt %in% c('?', '~'), NA, UniProt),
         UniProt = str_replace_all(UniProt, ' \\+ ', ','),
         UniProt = str_replace_all(UniProt, ' or ', ','))

# All FDA biomarkers with a UniProt ID
fda_biomarkers <- biomarkers |> 
  filter(!is.na(UniProt)) |> 
  pull(Protein.Name)

# FDA biomarkers found with MS
ms_biomarkers_uniprot <- biomarkers |> 
  filter(!is.na(UniProt)) |> 
  separate_longer_delim(UniProt, ',') |> 
  inner_join(filter(protein_metadata, UniProt %in% ms_proteins), by = 'UniProt') |> 
  group_by(Protein.Name) |>
  summarize(UniProt_MS = paste(UniProt, collapse = ';'),
            Gene.Name_MS = paste(Gene.Name, collapse = ';'))

ms_biomarkers <- ms_biomarkers_uniprot |>
  distinct(Protein.Name) |> 
  pull(Protein.Name)

# FDA biomarkers found with Olink
olink_biomarkers_uniprot <- biomarkers |>
  filter(!is.na(UniProt)) |>
  separate_longer_delim(UniProt, ',') |>
  inner_join(filter(protein_metadata, UniProt %in% olink_proteins), by = 'UniProt') |> 
  group_by(Protein.Name) |>
  summarize(UniProt_Olink = paste(UniProt, collapse = ';'),
            Assay_Olink = paste(Assay, collapse = ';'),
            OlinkID = paste(OlinkID, collapse = ';'))

olink_biomarkers <- olink_biomarkers_uniprot |>
  distinct(Protein.Name) |> 
  pull(Protein.Name)

# All biomarkers detected with MS or Olink
all_detected_biomarkers <- union(ms_biomarkers, olink_biomarkers)

# Calculate proportions
ms_prop <- length(intersect(ms_biomarkers, fda_biomarkers)) / length(fda_biomarkers)*100
olink_prop <- length(intersect(olink_biomarkers, fda_biomarkers)) / length(fda_biomarkers)*100
ms_olink_prop <- length(intersect(all_detected_biomarkers, fda_biomarkers)) / length(fda_biomarkers)*100

# Data for bar plot
plot_data <- data.frame(
  Protein_Set = factor(c(rep("MS", 2), rep("Olink", 2), rep("MS + Olink", 2)),
                       levels = c('MS', 'Olink', 'MS + Olink')),
  Source = factor(c('MS', 'FDA', 'Olink', 'FDA', 'MS + Olink', 'FDA'),
                  levels = c('MS', 'Olink', 'MS + Olink', 'FDA')),
  N = c(length(ms_biomarkers), length(fda_biomarkers) - length(ms_biomarkers),
        length(olink_biomarkers), length(fda_biomarkers) - length(olink_biomarkers),
        length(all_detected_biomarkers), length(fda_biomarkers) - length(all_detected_biomarkers)),
  Proportion = c(ms_prop, 100 - ms_prop, 
                 olink_prop, 100 - olink_prop, 
                 ms_olink_prop, 100 - ms_olink_prop)) |> 
  mutate(Label = paste0(round(Proportion), '%'),
         Label = ifelse(Source == 'FDA', NA, Label))

# Bar plot
plot_data |> 
  ggplot(aes(fill = Source, alpha = Source, y = N, x = Protein_Set)) +
  geom_col(position = position_stack(reverse = T)) +
  geom_text(aes(label = Label),
            position = position_stack(vjust = 0.5, reverse = T), size = 5/.pt,
            show.legend = F) +
  scale_y_continuous(expand = expansion(mult = c(0,0.05))) +
  theme_publ() +
  labs(x = 'Platform', y = 'Number of proteins',
       title = 'Coverage of FDA-approved plasma protein biomarkers') +
  theme(text = element_text(size = 5)) +
  scale_fill_manual(values = as.character(c(TDP_palette[c(1,4,13)], 'grey80'))) +
  scale_alpha_manual(values = c(0.8,1,1,1))

ggsave(filename = 'FDA_biomarker_coverage.pdf',
       path = path,
       width = 6,
       height = 6.5,
       units = 'cm')


# Proportional venn
n_ms <- length(setdiff(ms_biomarkers, olink_biomarkers))
n_olink <- length(setdiff(olink_biomarkers, ms_biomarkers))
n_ms_olink <- length(intersect(ms_biomarkers, olink_biomarkers))

pdf(file = file.path(path, 'FDA_biomarker_venn.pdf'),
    width = 4*0.393701, height = 5*0.393701)

p <- plot(euler(
  combinations = c('MS' = n_ms,
                   'Olink' = n_olink,
                   'MS&Olink' = n_ms_olink)),
  fills = list(fill = TDP_palette[c(1,4)], alpha = 0.5), 
  edges = list(lwd = 2.5, col = TDP_palette[c(1,4)]),
  quantities = list(fontsize = 7, font = 'regular', fontfamily = 'Helvetica'), 
  labels = list(fontsize = 8, font = 'regular', fontfamily = 'Helvetica'))

# Adjust plot to smaller size
p$vp$width <- unit(0.95, "npc")
p$vp$height <- unit(0.95, "npc")

p

dev.off()

# Result table
biomarkers_table <- biomarkers |> 
  dplyr::select(Protein.Name, UniProt) |>
  dplyr::rename(UniProt_FDA = UniProt) |>
  mutate(Detected_MS = Protein.Name %in% ms_biomarkers) |>
  left_join(ms_biomarkers_uniprot, by = 'Protein.Name') |>
  mutate(Detected_Olink = Protein.Name %in% olink_biomarkers) |>
  left_join(olink_biomarkers_uniprot, by = 'Protein.Name') |>
  mutate(Detected.Platform = case_when(
           Detected_MS & Detected_Olink ~ 'Both',
           Detected_MS & !Detected_Olink ~ 'MS',
           Detected_Olink & !Detected_MS ~ 'Olink',
           !Detected_Olink & !Detected_MS ~ 'None'),
         Detected.Platform = ifelse(
           is.na(UniProt_FDA), NA_character_, Detected.Platform))

write.csv(biomarkers_table, file.path(path, 'FDA_biomarkers_table.csv'))


##### Analysis on proteins with max 50% missing values #####

###### Prepare data ######

# MS data
ms_data_all <- read_MS_data(withReplicateSamples = TRUE)
ms_data_max50NA <- na_omit_prop(ms_data_all, 0.5, ignore_cols = 'UniProt')

# Olink Explore data
olink_data_LOD <- read_OlinkExplore_data(withBelowLOD = FALSE,
                                         removeAllLOD = TRUE)
olink_data_wide_LOD <- olink_long_to_wide(olink_data_LOD, id_cols = 'UniProt')
olink_data_wide_max50NA <- na_omit_prop(
  olink_data_wide_LOD, 0.5, ignore_cols = 'UniProt')

# Vectors of protein IDs
ms_proteins_max50NA <- unique(ms_data_max50NA$UniProt)
olink_proteins_max50NA <- unique(olink_data_wide_max50NA$UniProt)
all_proteins_max50NA <- union(ms_proteins_max50NA, olink_proteins_max50NA)
ms_unique_max50NA <- setdiff(ms_proteins_max50NA, olink_proteins_max50NA)
olink_unique_max50NA <- setdiff(olink_proteins_max50NA, ms_proteins_max50NA)

# HPA data for proteins detected with MS
hpa_ms_max50NA <- filter(hpa_ms, Uniprot %in% ms_proteins_max50NA)

# HPA data for proteins detected with Olink Explore
hpa_olink_max50NA <- filter(hpa_olink, Uniprot %in% olink_proteins_max50NA)

# All detected proteins found in HPA
hpa_detected_max50NA <- bind_rows(hpa_ms_max50NA, hpa_olink_max50NA) |> 
  distinct()

# Number and proportion of MS / Olink proteins not found in HPA
sum(!ms_proteins_max50NA %in% hpa_data$Uniprot)
sum(!ms_proteins_max50NA %in% hpa_data$Uniprot)/length(ms_proteins_max50NA)*100
sum(!olink_proteins_max50NA %in% hpa_data$Uniprot)
sum(!olink_proteins_max50NA %in% hpa_data$Uniprot)/length(olink_proteins_max50NA)*100


###### HPA annotations #######

cols <- c('Chromosome', 'Protein.class', 'Biological.process', 'Molecular.function',
          'Disease.involvement', 'Evidence', 'HPA.evidence', 'UniProt.evidence',
          'NeXtProt.evidence', 'RNA.tissue.specificity', 'RNA.tissue.distribution',
          'RNA.single.cell.type.specificity', 'RNA.single.cell.type.distribution',
          'RNA.cancer.specificity', 'RNA.cancer.distribution', 
          'RNA.brain.regional.specificity', 'RNA.brain.regional.distribution',
          'RNA.blood.cell.specificity', 'RNA.blood.cell.distribution',
          'RNA.blood.lineage.specificity', 'RNA.blood.lineage.distribution',
          'Blood.expression.cluster', 'Tissue.expression.cluster',
          'Subcellular.location', 'Secretome.location', 'Subcellular.main.location',
          'Subcellular.additional.location', 'Enriched.tissue')

# Run analysis
hpa_ann_freq_max50NA <- setNames(
  lapply(cols, 
         ann_freq_summary, 
         bg_ann = hpa_detected_max50NA, 
         ms_ann = hpa_ms_max50NA, 
         olink_ann = hpa_olink_max50NA), 
  cols)

hpa_ann_freq_df_max50NA <- do.call(rbind.data.frame, hpa_ann_freq_max50NA)

write.csv(hpa_ann_freq_df_max50NA, file.path(path, 'HPA_ann_frequencies_max50NA.csv'),
          row.names = FALSE)

# Plot results

hpa_freq_plots_max50NA <- setNames(lapply(hpa_ann_freq_max50NA, freq_barplot), cols)

p <- wrap_plots(hpa_freq_plots_max50NA[c('Protein.class', 'Secretome.location', 'Enriched.tissue', 'Subcellular.location')]) +
  plot_layout(ncol = 2, axis_titles = 'collect')

freq_plot_max50NA <- guide_area() + p +
  plot_layout(guides = 'collect', nrow = 2, heights = c(0.4,10), design = 'A\nB')

freq_plot_max50NA

ggsave(freq_plot_max50NA, 
       filename = 'HPA_freq_barplots_max50NA.pdf',
       path = path,
       width = 18.3,
       height = 12,
       units = 'cm')

###### GO terms ######

# Target protein lists
input_list_max50NA <- list('MS' = ms_proteins_max50NA,
                           'Olink' = olink_proteins_max50NA)

# Background protein list
background_max50NA <- all_proteins_max50NA

GO_res_BP_max50NA <- GO.analysis(input_list = input_list_max50NA, 
                                 background = background_max50NA,
                                 ont = 'BP')


# Result data frame
GO_res_BP_max50NA_df <- as.data.frame(GO_res_BP_max50NA)
write.csv(GO_res_BP_max50NA_df, file.path(path, 'GO_BP_res_max50NA.csv'))

# Plot venn diagram of enriched GO terms
res_ms <- GO_res_BP_max50NA@compareClusterResult %>% filter(Cluster == 'MS')
res_olink <- GO_res_BP_max50NA@compareClusterResult %>% filter(Cluster == 'Olink')
ggvenn(list('MS' = res_ms$ID, 'Olink' = res_olink$ID),
       fill_color = c("#21a2ac","#f38895"),
       stroke_color = rep(c("#21a2ac","#f38895"), each = 100),
       stroke_size = 1,
       text_size = 6/.pt,
       set_name_size = 7/.pt) + 
  coord_cartesian(xlim = c(-2, 2), ylim = c(-2, 2)) +
  labs(title = paste('ORA of GO BP terms - max 50% missing values')) +
  theme(title = element_text(size = 8))

ggsave(filename = 'GO_BP_venn_max50NA.pdf',
       path = path,
       width = 5.5,
       height = 6.5,
       units = 'cm')


GO_BP_dotplot_max50NA <- plot.enr.dotplot(GO_res_BP_max50NA)

ggsave(GO_BP_dotplot_max50NA,
       filename = 'GO_BP_dotplot_max50NA.pdf',
       path = path,
       width = 18.3,
       height = 6.5,
       units = 'cm')


###### Coverage of clinically used biomarkers ######

ms_biomarkers_uniprot_max50NA <- biomarkers |> 
  filter(!is.na(UniProt)) |> 
  separate_longer_delim(UniProt, ',') |> 
  inner_join(filter(protein_metadata, UniProt %in% ms_proteins_max50NA), by = 'UniProt') |> 
  group_by(Protein.Name) |>
  summarize(UniProt_MS = paste(UniProt, collapse = ','),
            Gene.Name_MS = paste(Gene.Name, collapse = ','))

ms_biomarkers_max50NA <- ms_biomarkers_uniprot_max50NA |>
  distinct(Protein.Name) |> 
  pull(Protein.Name)

olink_biomarkers_uniprot_max50NA <- biomarkers |>
  filter(!is.na(UniProt)) |>
  separate_longer_delim(UniProt, ',') |>
  inner_join(filter(protein_metadata, UniProt %in% olink_proteins_max50NA), by = 'UniProt') |> 
  group_by(Protein.Name) |>
  summarize(UniProt_Olink = paste(UniProt, collapse = ';'),
            Assay_Olink = paste(Assay, collapse = ';'),
            OlinkID = paste(OlinkID, collapse = ';'))

olink_biomarkers_max50NA <- olink_biomarkers_uniprot_max50NA |>
  distinct(Protein.Name) |> 
  pull(Protein.Name)

all_detected_biomarkers_max50NA <- union(ms_biomarkers_max50NA, olink_biomarkers_max50NA)

# Calculate proportions
ms_prop_max50NA <- length(intersect(ms_biomarkers_max50NA, fda_biomarkers)) / length(fda_biomarkers)*100
olink_prop_max50NA <- length(intersect(olink_biomarkers_max50NA, fda_biomarkers)) / length(fda_biomarkers)*100
ms_olink_prop_max50NA <- length(intersect(all_detected_biomarkers_max50NA, fda_biomarkers)) / length(fda_biomarkers)*100

plot_data <- data.frame(
  Protein_Set = factor(c(rep("MS", 2), rep("Olink", 2), rep("MS + Olink", 2)),
                       levels = c('MS', 'Olink', 'MS + Olink')),
  Source = factor(c('MS', 'FDA', 'Olink', 'FDA', 'MS + Olink', 'FDA'),
                  levels = c('MS', 'Olink', 'MS + Olink', 'FDA')),
  N = c(length(ms_biomarkers_max50NA), length(fda_biomarkers) - length(ms_biomarkers_max50NA),
        length(olink_biomarkers_max50NA), length(fda_biomarkers) - length(olink_biomarkers_max50NA),
        length(all_detected_biomarkers_max50NA), length(fda_biomarkers) - length(all_detected_biomarkers_max50NA)),
  Proportion = c(ms_prop_max50NA, 100 - ms_prop_max50NA, 
                 olink_prop_max50NA, 100 - olink_prop_max50NA, 
                 ms_olink_prop_max50NA, 100 - ms_olink_prop_max50NA)) |> 
  mutate(Label = paste0(round(Proportion), '%'),
         Label = ifelse(Source == 'FDA', NA, Label))

plot_data |> 
  ggplot(aes(fill = Source, alpha = Source, y = N, x = Protein_Set)) +
  geom_col(position = position_stack(reverse = T)) +
  geom_text(aes(label = Label),
            position = position_stack(vjust = 0.5, reverse = T), size = 5/.pt,
            show.legend = F) +
  scale_y_continuous(expand = expansion(mult = c(0,0.05))) +
  theme_publ() +
  labs(x = 'Platform', y = 'Number of proteins',
       title = 'Coverage of FDA-approved plasma protein biomarkers') +
  theme(text = element_text(size = 5)) +
  scale_fill_manual(values = as.character(c(TDP_palette[c(1,4,13)], 'grey80'))) +
  scale_alpha_manual(values = c(0.8,1,1,1))

ggsave(filename = 'FDA_biomarker_coverage_max50NA.pdf',
       path = path,
       width = 6,
       height = 6.5,
       units = 'cm')

