
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(data.table)
library(heatmaply)
library(RColorBrewer)
library(dplyr)

 
get_peptide_and_correlation_numeric <- function(fasta_list, peptide_seq_list, exclude = TRUE){
  peptide_indices <- list()
  for( i in 1:length(peptide_seq_list$peptides)){
    peptide_indices[[rownames(peptide_seq_list[i,])]] <- data.frame(nchar(gsub(paste0(peptide_seq_list$peptides[i], ".*"), "", fasta_list)),
                                                                    peptide_seq_list$peptide[i],
                                                                    peptide_seq_list$N[i],
                                                                    peptide_seq_list$target_correlation[i])
    names(peptide_indices[[rownames(peptide_seq_list[i,])]]) <- c("start_position", "peptide", "N", "target_correlation")
    
    if(exclude == TRUE){
      
      if(peptide_indices[[rownames(peptide_seq_list[i,])]][1] == nchar(fasta_list)){
        
        peptide_indices[[rownames(peptide_seq_list[i,])]] <- NULL
        
      }
      
    }
  }
  
  return(peptide_indices)
}

create_annotation <- function(ann_df, seq_length) {
  
  # transform into data frame with one column per annotation
  # number of rows equal to sequence length
  # and the annotation values filled in each row/column
  ann_mat <- matrix(NA, nrow = seq_length, ncol = length(unique(ann_df$feature)))
  unique_features <- unique(ann_df$feature)
  colnames(ann_mat) <- unique_features
  # fill in the annotation values
  for (i in 1:length(unique_features)) {
    ann_df_subset <- ann_df[ann_df$feature == unique_features[i],]
    for (j in 1:nrow(ann_df_subset)) {
      ann_mat[ann_df_subset$start[j]:ann_df_subset$end[j], i] <- unique_features[i]
    }
  }
  
  ann_mat
  
}


heatmap_plot <- function(gene_name, peptide_seq_list_file, fasta_data,
                         col_annotations = NULL, correlation_palette = c("red", "gray", "blue"),
                         hm_height = NULL, ann_height = NULL, ann_palette = NULL) {
  
  # Get fasta seq
  fasta_list <- unique(fasta_data$fasta[fasta_data$Gene.Name ==gene_name])
  
  peptide_seq_list <- readRDS(peptide_seq_list_file)
  peptide_seq_list <- peptide_seq_list[[gene_name]]
  peptide_indices <- tryCatch({get_peptide_and_correlation_numeric(fasta_list, peptide_seq_list)},
                              error = function(e){0})
  
  # Initialize meta_hmap with two columns: Amino acid residue, Olink target correlation
  meta_hmap <- matrix(nrow = nchar(fasta_list[[1]]), ncol = 2)
  colnames(meta_hmap) <- c("Amino acid residue", "Olink target correlation")
  meta_hmap[,'Amino acid residue'] <- as.numeric(1:nchar(fasta_list[[1]]))
  
  # Function to process peptides and handle overlaps (as before)
  process_peptides <- function(col) {
    overlap_peptides <- list()
    for (k in seq_along(peptide_indices)) {
      start <- as.numeric(peptide_indices[[k]][["start_position"]])
      end <- start + nchar(peptide_indices[[k]][["peptide"]]) - 1
      added_to_overlap <- FALSE
      
      if(all(is.na(meta_hmap[start:end, col]))){
        meta_hmap[start:end, col] <- rep(peptide_indices[[k]][["target_correlation"]], times = length(start:end))
      } else {
        added_to_overlap = TRUE
      }
      
      if (added_to_overlap) {
        overlap_peptides[[length(overlap_peptides) + 1]] <- peptide_indices[[k]]
      }
    }
    return(list(overlap_peptides=overlap_peptides, meta_hmap=meta_hmap))
  }
  
  # Manage overlaps, starting with the correlation column (2)
  col <- 2
  repeat {
    result <- process_peptides(col)
    current_overlaps <- result$overlap_peptides
    meta_hmap <- result$meta_hmap
    if (length(current_overlaps) == 0) {
      break
    }
    peptide_indices <- current_overlaps
    # If more overlaps remain, add another correlation column
    col <- ncol(meta_hmap) + 1
    meta_hmap <- cbind(meta_hmap, NA)
    colnames(meta_hmap)[col] <- paste("Correlation", col - 2)
  }
  
  meta_hmap <- as.data.frame(meta_hmap, stringsAsFactors = FALSE)
  
  # We now have: 
  # "Amino acid residue", "Olink target correlation", maybe more "Correlation X" columns
  
  # Convert this to numeric
  meta_hmap[,"Olink target correlation"] <- as.numeric(meta_hmap[,"Olink target correlation"])
  
  cor_data <- meta_hmap[,grepl("correlation", colnames(meta_hmap))|grepl("Correlation", colnames(meta_hmap)), drop=FALSE]
  # Convert cor_data into a numeric matrix for heatmaply
  cor_mat <- as.matrix(cor_data)
  rownames(cor_mat) <- meta_hmap[,"Amino acid residue"]
  
  if (is.null(hm_height)) {
    hm_height <- (ncol(cor_mat) + 1) * 30
  }
  
  if (is.null(ann_height)) {
    ann_height <- (ncol(col_annotations) + 1) * 20
  }
  
  total_height <- hm_height + ann_height

  heatmap <- heatmaply(
    t(data.matrix(cor_mat)),
    scale = "none",
    limits = c(-1,1),
    showticklabels = c(TRUE, FALSE),
    colors = correlation_palette,
    ColSideColors = col_annotations,
    col_side_palette = ann_palette,
    colv = F,
    Rowv = F,
    plot_method ='plotly',
    height = total_height,
    subplot_heights = c(ann_height/total_height, hm_height/total_height),
    dendrogram = "none"
  ) %>%
    colorbar(
      which = 2,
      title = "MS-Olink correlation", 
      titlefont = list(size = 9), 
      tickfont = list(size = 8),
      x = 1,
      y = 0.2,
      len = 100,
      lenmode = 'pixels',
      thickness = 15
    ) %>%
    layout(
      xaxis = list(
        title = 'Sequence Position',
        titlefont = list(size = 9),
        tickfont = list(size = 8),
        tickvals = floor(seq(1, nrow(cor_mat), length.out = 10)),  # Set tick positions every 100
        ticktext = floor(seq(1, nrow(cor_mat), length.out = 10)),   # Set corresponding labels
        tickangle = 0
      ),
      yaxis = list(
        tickfont = list(size = 8),
        side = 'right'
      ),
      yaxis2 = list(
        title = 'MS Peptides',
        titlefont = list(size = 9),
        tickfont = list(size = 8)
      )
    )
  
  heatmap$x$data[[1]]$showscale <- FALSE
  
  return(heatmap)
  
}

jenks_density_plot <- function(gene_name, peptide_seq_list_file, fasta_data) {
  library(ggplot2)
  library(dplyr)
  
  # Load peptide sequences for the gene of interest
  peptide_seq_list <- readRDS(peptide_seq_list_file)
  if (!gene_name %in% names(peptide_seq_list)) {
    stop("IOI not found in peptide_seq_list.")
  }
  
  gene_peptides <- peptide_seq_list[[gene_name]]
  
  # Get FASTA sequence
  fasta_seq <- unique(fasta_data$fasta[fasta_data$Gene.Name ==gene_name])
  
  jenks_df <- gene_peptides %>% 
    mutate(peptide_id = row_number()) %>% 
    mutate(jenks_class = cor_intervals(as.numeric(target_correlation)))
  
  # jenks_df now has a jenks_class column
  # Now, we need to map each peptide to its amino acid positions
  
  # For each peptide, we find its start position in the sequence
  # The start position is computed by removing the peptide and counting length
  # The first occurrence of the peptide in fasta_seq gives start position
  # But be aware of multiple occurrences. If multiple matches occur, handle accordingly!!!!! RECHECK THIS LOGIC
  
  # Assume the first occurrence is the correct one:
  create_position_data <- function(peptide, jenks_class, fasta_seq, used_positions) {
    # Find the start position of the peptide in fasta_seq
    # We can use gregexpr for this:
    match_pos <- gregexpr(peptide, fasta_seq, fixed = TRUE)[[1]]
    if (match_pos[1] == -1) {
      # Peptide not found, skip
      return(NULL)
    }
    
    # Find the first match that does not overlap with used_positions
    for (pos in match_pos) {
      start_pos <- pos
      end_pos <- start_pos + nchar(peptide) - 1
      if (!any(used_positions >= start_pos & used_positions <= end_pos)) {
        # Mark these positions as used
        used_positions[start_pos:end_pos] <- TRUE
        # Create a data frame of all positions covered by this peptide
        return(data.frame(position = start_pos:end_pos,
                          jenks_class = jenks_class,
                          stringsAsFactors = FALSE))
      }
    }
    # If all matches overlap, skip
    return(NULL)
  }
  
  # Initialize a vector to keep track of used positions
  used_positions <- rep(FALSE, nchar(fasta_seq))
  
  # Apply this to all peptides
  position_data_list <- lapply(seq_len(nrow(jenks_df)), function(i) {
    current_peptide <- jenks_df$peptides[i]
    current_class <- jenks_df$jenks_class[i]
    
    # Check if peptide is non-empty
    if (is.na(current_peptide) || current_peptide == "") return(NULL)
    create_position_data(current_peptide, current_class, fasta_seq, used_positions)
  })
  
  # Combine all
  position_data <- do.call(rbind, position_data_list)
  
  # If no data, return NULL or a simple message
  if (is.null(position_data) || nrow(position_data) == 0) {
    warning("No peptide position data available for plotting.")
    return(NULL)
  }
  
  # Now we have a data frame with columns: position, jenks_class
  # We can make a density plot by jenks_class
  p <- ggplot(position_data, aes(x = position, fill = jenks_class,
                                 text = paste("Class:", jenks_class))) +
    geom_density(alpha = 0.5, lwd = 0.25) +
    scale_fill_manual(values = c("No correlation" = "grey",
                                 "Weak correlation" = "#9ccfe7",  
                                 "Moderate correlation" = "#629db8", 
                                 "Strong correlation" = "#206e8c")) +
    theme_minimal() +
    labs(
      title = paste("Peptide Density by Class -",gene_name),
      x = "Sequence Position",
      y = "Density",
      fill = "Correlation Category"
    ) +
    theme(
      legend.key.size = unit(0.5, "lines"),
    )
  
  return(p)
}

# Combined Plot Function
heatmap_density_plot <- function(gene_name, peptide_seq_list_file, fasta_data, 
                                 annotations = NULL, correlation_palette = NULL,
                                 hm_height = NULL, ann_height = NULL, density_height = 150,
                                 path = NULL, ann_palette = NULL) {
  
  if (is.null(correlation_palette)) {
    correlation_palette <- colorRampPalette(rev(brewer.pal(11, 'Spectral')))(100)
  }
  
  # Generate the Heatmaply Plot
  hm <- heatmap_plot(
   gene_name = gene_name,
   peptide_seq_list_file = peptide_seq_list_file,
   fasta_data = fasta_data,
   col_annotations = annotations,
   correlation_palette = correlation_palette,
   hm_height = hm_height,
   ann_height = ann_height,
   ann_palette = ann_palette
  )
  
  # Generate the Density Plot using ggplot2
  density_plot_gg <- jenks_density_plot(
   gene_name = gene_name,
    peptide_seq_list_file = peptide_seq_list_file,
    fasta_data = fasta_data
  )
  
  # Convert ggplot2 Density Plot to Plotly Object
  density_plotly <- ggplotly(density_plot_gg, height = density_height) %>%
    layout(showlegend = TRUE,
           legend = list(
             itemwidth = 0.01,
             tracegroupgap = 0,
             font = list(size = 8), 
             y = 1,
             yanchor = 'container',
             title = list(font = list(size = 9))),
           xaxis = list(title = "", showticklabels = FALSE), # remove x-axis
           yaxis = list(showticklabels = FALSE)) # remove y-axis labels
  
  # Ensure the heatmap and density plot share the same x-axis range
  heatmap_xrange <- hm$x$layout$xaxis$range
  density_plotly <- layout(density_plotly, xaxis = list(range = heatmap_xrange))
  
  heatmap_height <- hm$height
  total_height <- heatmap_height + density_height
  
  # Combine Heatmap and Density Plot using subplot
  combined_plot <- subplot(
    density_plotly,
    hm,
    nrows = 2, 
    shareX = TRUE, 
    margin = 0,
    heights = c(density_height/total_height, heatmap_height/total_height)
  ) %>%
    layout(
      title = list(
        text = paste("MS Peptide Mapping & Correlation to Olink -", gene_name),
        font = list(size = 11)),
      xaxis = list(title = "Sequence Position"),
      margin = list(t = 30, b = 25, l = 50, r = 0),
      annotations = list(
        list(
          x = -0.01, y = 0.5,
          showarrow = FALSE,
          text = 'Density',
          xref = 'paper', yref = 'y1 domain',
          xanchor = 'right', yanchor = 'middle',
          font = list(size = 9),
          textangle = -90
        ),
        list(
          x = -0.01, y = 0.5,
          showarrow = FALSE,
          text = 'Sequence Features',
          xref = 'paper', yref = 'y3 domain',
          xanchor = 'right', yanchor = 'middle',
          font = list(size = 9),
          textangle = -90
        ),
        list(
          x = -0.01, y = 0.5,
          showarrow = FALSE,
          text = 'MS Peptides',
          xref = 'paper', yref = 'y2 domain',
          xanchor = 'right', yanchor = 'middle',
          font = list(size = 9),
          textangle = -90
        )
      )
    )
  
  # Save image
  if (!is.null(path)) {
    orca(
      combined_plot,
      file = file.path(path, paste0("peptide_plot_", gene_name, ".svg")),
      width = 600, 
      height = total_height
    )
  }
  
  return(combined_plot)
}


cor_intervals <- function(x) {
  cut(x, breaks = c(-1, 0.3, 0.5, 0.7, 1),
      labels = c('No correlation','Weak correlation', 
                 'Moderate correlation', 'Strong correlation'),
      include.lowest = TRUE, right = FALSE)
}
