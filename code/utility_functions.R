
read_metadata <- function(withReplicates = FALSE, samples = c('overlapping', 'MS')) {
  
  samples <- match.arg(samples)
  
  if (withReplicates) {
    
    if (samples == 'MS') {
      metadata <- read.delim('data/metadata/sample_metadata.txt') |> 
        filter(!str_detect(Sample.ID, 'IS'))
    } else {
      stop('Invalid combination of arguments')
    }
  } else {
    if (samples == 'MS') {
      metadata <- read.delim('data/metadata/sample_metadata_noReplicates.txt')
    } else if (samples == 'overlapping') {
      metadata <- read.delim('data/metadata/sample_metadata_overlap.txt')
    }
  }
  return(metadata)
}

read_MS_data <- function(withIS = FALSE, withReplicateSamples = FALSE) {
  # withIS: Include internal standards
  # withReplicateSamples: Include replicate samples
  # allProteins: Also include proteins with all NA values
  
  # MS data with internal standards and replicate samples
  if (!withReplicateSamples & !withIS) {
    ms_data <- read.delim('data/processed_data/protein_table_noReplicates.txt')
  } else {
    ms_data <- read.delim('data/processed_data/protein_table_all.txt')
    # Remove internal standards
    if (!withIS) {
      ms_data <- ms_data |> 
        dplyr::select(!contains('IS'))
    }
    if (!withReplicateSamples) {
      stop('Invalid combination of arguments')
    }
  }
  return(ms_data)
}

read_OlinkExplore_data <- function(withControlSamples = FALSE, withBelowLOD = TRUE, 
                                   withQCWarning = TRUE, withReplicateAssays = FALSE,
                                   removeAllLOD = TRUE) {
  
  if (withReplicateAssays) {
    olink_data <- read.delim('data/processed_data/OlinkExplore_data.txt')
  } else {
    olink_data <- read.delim('data/processed_data/OlinkExplore_data_noReplicateAssays.txt')
  }

  if (!withControlSamples) {
    olink_data <- olink_data |>  
      filter(!str_detect(Sample.ID, 'Control'))
  }
  
  if (removeAllLOD) {
    # Calculate number of <LOD values per protein
    to_remove <- olink_data |> 
      group_by(OlinkID) |> 
      # Calculate number of <LOF values across samples not containing "Control"
      filter(!str_detect(Sample.ID, 'Control')) |>
      summarise(Below.LOD = sum(NPX < LOD),
                n = n()) |> 
      filter(Below.LOD == n) |> 
      pull(OlinkID)
      
    olink_data <- olink_data |>
      filter(!OlinkID %in% to_remove)
  }
  
  if (!withBelowLOD) {
    olink_data <- olink_data |> 
      mutate(NPX = ifelse(NPX < LOD, NA_real_, NPX))
  }
  
  if (!withQCWarning) {
    olink_data <- olink_data |> 
      mutate(NPX = ifelse(QC_Warning=='WARN'|Assay_Warning=='WARN', NA_real_, NPX))
  }
  
  return(olink_data)
  
}

olink_long_to_wide <- function(NPX_data, id_cols = c()) {
  
  library(tidyverse)
  
  sample_col <- case_when('SampleID' %in% colnames(NPX_data) ~ 'SampleID',
                          'Sample.ID' %in% colnames(NPX_data) ~ 'Sample.ID',
                          TRUE ~ NA_character_)
  
  if (is.na(sample_col)) stop('Sample.ID/SampleID column not found')
  
  NPX_data |> pivot_wider(id_cols = all_of(id_cols),
                           names_from = all_of(sample_col),
                           values_from = 'NPX') |>
    as.data.frame()
  
}

olink_proteins_per_panel <- function(NPX_data) {
  
  proteins_per_panel <- NPX_data |>
    group_by(Panel) |>
    summarise(total.proteins = length(unique(OlinkID))) |>
    as.data.frame()
}

na_omit_prop <- function(x, max_na_prop = 0, ignore_cols = NULL) {
  
  if (is.null(ignore_cols)) {
    x[rowSums(is.na(x))/ncol(x) <= max_na_prop,] 
  } else {
    x[rowSums(is.na(dplyr::select(x, -all_of(ignore_cols))))/ncol(dplyr::select(x, -all_of(ignore_cols))) <= max_na_prop,]
    }
}


group_summary <- function(data, summary_var, ...) {
  
  library(tidyverse)
  
  data |> 
    group_by(...) |> 
    summarize(Min = min({{summary_var}}, na.rm = TRUE),
              Q1 = quantile({{summary_var}}, na.rm = TRUE)[[2]],
              Median = median({{summary_var}}, na.rm = TRUE),
              Mean = mean({{summary_var}}, na.rm = TRUE),
              Q3 = quantile({{summary_var}}, na.rm = TRUE)[[4]],
              Max = max({{summary_var}}, na.rm = TRUE),
              IQR = IQR({{summary_var}}, na.rm = TRUE),
              SD = sd({{summary_var}}, na.rm = TRUE),
              N = n())
  
}

separate_wider_mixed <- function(data, col_name, sep, include_name = FALSE) {
  
  col_name_str <- rlang::as_name(rlang::ensym(col_name))
  col_name_str <- str_replace_all(col_name_str, '\\.', ' ')
  
  # Split the specified column by the separator
  data |>
    mutate({{col_name}} := strsplit({{col_name}}, sep)) |>
    unnest({{col_name}}) |>
    mutate(value = {{col_name}}) |>
    mutate(value = ifelse(
      include_name & !is.na(value), 
      paste0(col_name_str, ': ', {{col_name}}), 
      {{col_name}})) |>
    pivot_wider(names_from = {{col_name}},
                values_from = value,
                values_fill = NA,
                names_repair = 'universal_quiet')
  
  
}