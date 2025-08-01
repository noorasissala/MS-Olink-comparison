


library(ggplot2)
library(colorspace)
library(RColorBrewer)

#### Color palettes ####

TDP_palette <- c('#21A2AC', '#49BDC4', '#74CDD2', 
                 '#F38895', '#F5A5B5', '#F9C8D2', 
                 '#FDAB5E', '#FDBA7A', '#FDC691', 
                 '#9DD177', '#B5DC98', '#C8E5B2',
                 '#C3ABE3', '#D2C0EA', '#DDD0F0')

names(TDP_palette) <- c(rep('MS', 3), rep('Olink', 3), rep(NA_character_, 9))

panel_pal <- rep(brewer.pal(n = 4, name = 'Set2'), each = 2)
panel_pal <- desaturate(lighten(panel_pal, amount = 0.1), amount = 0.2)
names(panel_pal) <- c('Cardiometabolic', 'Cardiometabolic II', 
                      'Inflammation', 'Inflammation II', 
                      'Neurology', 'Neurology II',     
                      'Oncology', 'Oncology II')

blue_palette <- c('#206e8c', '#4485a2', '#629db8', '#7fb6cf', '#9ccfe7', '#b9e9ff')


#### Themes ####

theme_publ <- function() {
  # ggplot theme for publication figures
  
  theme_classic(base_size = 7, base_family = 'sans') +
  theme(legend.position = 'top',
        legend.key.size = unit(0.3, 'cm'),
        axis.text = element_text(size = 6),
        legend.text = element_text(size = 6),
        strip.text = element_text(hjust = 0, size = 7),
        strip.background = element_rect(fill = 'grey90', colour = NA))
}


#### Basic plots #####

ggplot.histogram <- function(data, X_var, Fill_var = NULL, bin_prop = NULL, 
                             binwidth = NULL, median_line = TRUE, median_text = TRUE, 
                             median_line_params = list(), median_text_params = list(),
                             title = NULL, xlab = NULL, ylab = NULL, text_hpos = 'left',
                             text_vpos = 'top', Facet_var = NULL, facet_scales = 'fixed', 
                             ...) {
  
  library(ggplot2)
  library(stringr)
  library(ggpp)
  
  X_var_str <- rlang::as_label(ensym(X_var))
  try(Fill_var <- ensym(Fill_var), silent = TRUE)
  if (is.null(Fill_var)) Fill_var_str <- NULL else Fill_var_str <- as_label(Fill_var)
  try(Facet_var <- ensym(Facet_var), silent = TRUE)
  if (is.null(Facet_var)) Facet_var_str <- NULL else Facet_var_str <- as_label(Facet_var)
  
  # Default parameters
  params <- list(...)
  if (is.null(params$boundary)) params$boundary <- 0
  if (is.null(params$alpha)) params$alpha <- ifelse(is.null(Fill_var), 1, 0.6)
  if (is.null(params$position)) params$position <- 'identity'
  
  if (is.null(median_line_params$linetype)) median_line_params$linetype <- 'dashed'
  if (is.null(median_text_params$size)) median_text_params$size <- 6/.pt
  
  # Calculate bin width
  if (!is.null(bin_prop)) binwidth <- (ceiling(range(data[[X_var_str]])[2])-floor(range(data[[X_var_str]])[1]))*bin_prop
  params$binwidth <- binwidth
  
  # Generate histogram
  p <- ggplot(data, aes(x = {{X_var}}, fill = {{Fill_var}}, color = {{Fill_var}})) +
    do.call(geom_histogram, params) +
    labs(title = title,
         x = ifelse(is.null(xlab), X_var_str, xlab),
         y = ifelse(is.null(ylab), 'Count', ylab)) +
    theme_publ() +
    scale_y_continuous(expand = expansion(c(0,0.05)))
  
  # Facet plot
  if (!is.null(Facet_var)) {
    p <- p + facet_wrap(vars({{Facet_var}}), scales = facet_scales)
  }
    
  # Add median line(s)
  if (median_line) {
    p <- p +
      do.call('geom_vline', modifyList(
        list(data = data %>% 
               group_by({{Facet_var}}, {{Fill_var}}) %>% 
               summarize(Median = median({{X_var}}, na.rm = T)),
             mapping = aes(xintercept = Median, color = {{Fill_var}})),
        median_line_params))
  }
  
  # Add median annotation(s)
  
  if (median_text) {
    
    if (!is.null(Fill_var)) {
      
      label_df <- data %>%
        group_by({{Fill_var}}, {{Facet_var}}) %>%
        summarize(Median = round(median({{X_var}}, na.rm = T), 2),
                  Label = str_glue('Median {group} = {Median}',
                                   group = cur_group()[[Fill_var_str]]))
      
      facet_n <- ifelse(is.null(Facet_var_str), 1, length(unique(data[[Facet_var_str]])))
      label_df$npcy <- calculate.label.npcy(y = text_vpos, 
                                            group_n = length(unique(data[[Fill_var_str]])),
                                            facet_n = facet_n)
      label_df$npcx <- text_hpos
    }
    
    else {
      
      label_df <- data %>%
        group_by({{Facet_var}}) %>%
        summarize(Median = round(median({{X_var}}, na.rm = T), 2),
                  Label = str_glue("Median = {Median}"),
                  npcy = text_vpos,
                  npcx = text_hpos)
      }
    
    p <- p + do.call('geom_text_npc', modifyList(
      list(data = label_df, 
           mapping = aes(npcx = npcx, npcy = npcy, label = Label, color = {{Fill_var}})), 
      median_text_params))
    }
  
  return(p)
    
}


ggplot.density <- function(data, X_var, Fill_var = NULL, median_line = TRUE, 
                           median_text = TRUE, median_line_params = list(), 
                           median_text_params = list(), title = NULL, xlab = NULL, 
                           ylab = NULL, Facet_var = NULL, facet_scales = 'fixed', 
                           text_hpos = 'left', ...) {
  
  library(tidyverse)
  library(ggplot2)
  library(stringr)
  library(ggpp)
  
  params <- list(...)
  
  X_var_str <- rlang::as_label(ensym(X_var))
  try(Fill_var <- ensym(Fill_var), silent = TRUE)
  if (is.null(Fill_var)) Fill_var_str <- NULL else Fill_var_str <- as_label(Fill_var)
  try(Facet_var <- ensym(Facet_var), silent = TRUE)
  if (is.null(Facet_var)) Facet_var_str <- NULL else Facet_var_str <- as_label(Facet_var)
  
  if (is.null(median_line_params$linetype)) median_line_params$linetype <- 'dashed'
  if (is.null(median_line_params$lwd)) median_line_params$lwd <- 0.3
  if (is.null(params$alpha)) params$alpha <- 0.5
  
  if (is.null(median_text_params$size)) median_text_params$size <- 6/.pt
  
  # Generate histogram
  p <- ggplot(data, aes(x = {{X_var}}, fill = {{Fill_var}}, color = {{Fill_var}})) +
    do.call('geom_density', params) +
    labs(title = title,
         x = ifelse(is.null(xlab), X_var_str, xlab),
         y = ifelse(is.null(ylab), 'Density', ylab)) +
    theme_publ()
  
  # Facet plot
  if (!is.null(Facet_var)) {
    p <- p + facet_wrap(vars({{Facet_var}}), scales = facet_scales)
  }
  
  # Add median line(s)
  if (median_line) {
    p <- p +
      do.call('geom_vline', modifyList(
        list(data = data %>% 
               group_by({{Facet_var}}, {{Fill_var}}) %>% 
               summarize(Median = median({{X_var}}, na.rm = T)),
             mapping = aes(xintercept = Median, color = {{Fill_var}})),
        median_line_params))
  }
  
  # Add median annotation(s)
  
  if (median_text) {
    
    if (!is.null(Fill_var)) {
      
      label_df <- data %>%
        group_by({{Fill_var}}, {{Facet_var}}) %>%
        summarize(Median = round(median({{X_var}}, na.rm = T), 2),
                  Label = str_glue('Median {group} = {Median}',
                                   group = cur_group()[[Fill_var_str]]))
      
      facet_n <- ifelse(is.null(Facet_var_str), 1, length(unique(data[[Facet_var_str]])))
      label_df$npcy <- calculate.label.npcy(y = 'top', 
                                            group_n = length(unique(data[[Fill_var_str]])),
                                            facet_n = facet_n)
      
      label_df$npcx <- text_hpos
      
    }
    
    else {
      
      label_df <- data %>%
        group_by({{Facet_var}}) %>%
        summarize(Median = round(median({{X_var}}, na.rm = T), 2),
                  Label = str_glue("Median = {Median}"),
                  npcy = 'top',
                  npcx = text_hpos)
    }
    
    p <- p + do.call('geom_text_npc', modifyList(
      list(data = label_df, 
           mapping = aes(npcx = npcx, npcy = npcy, label = Label, color = {{Fill_var}})), 
      median_text_params)) +
      scale_y_continuous(expand = expansion(c(0,0.05)))
  }
  
  return(p)
  
}

ggplot.barplot <- function(data, X_var, Y_var, Fill_var = NULL, 
                           title = NULL, xlab = NULL, ylab = NULL, 
                           Facet_var = NULL, facet_scales = 'fixed', ...) {
  
  library(ggplot2)
  library(stringr)
  
  X_var_str <- rlang::as_label(ensym(X_var))
  Y_var_str <- rlang::as_label(ensym(Y_var))
  try(Fill_var <- ensym(Fill_var), silent = TRUE)
  try(Facet_var <- ensym(Facet_var), silent = TRUE)
  
  # Generate plot
  p <- ggplot(data, aes(x = {{X_var}}, y = {{Y_var}}, fill = {{Fill_var}})) +
    geom_col(...) +
    labs(x = ifelse(is.null(xlab), X_var_str, xlab), 
         y = ifelse(is.null(ylab), Y_var_str, ylab),
         title = ifelse(is.null(title), str_glue('{Y_var_str} by {X_var_str}'), title)) +
    theme_bw()
  
  if (!is.null(Facet_var)) {
    p <- p + facet_wrap(vars({{Facet_var}}), scales = facet_scales)}
  
  return(p)

}


ggplot.scatterplot <- function(data, X_var, Y_var, Color_var = NULL, 
                               title = NULL, xlab = NULL, ylab = NULL, density = FALSE,
                               Facet_var = NULL, facet_scales = 'fixed', ...) {
  
  library(ggplot2)
  library(stringr)
  library(ggpointdensity)
  
  X_var_str <- rlang::as_label(ensym(X_var))
  Y_var_str <- rlang::as_label(ensym(Y_var))
  try(Facet_var <- ensym(Facet_var), silent = TRUE)
  
  # Set default parameters
  params <- list(...)
  point_params <- modifyList(list(alpha = 0.6), params)
  point_params$size <- coalesce(point_params$size, 0.5)

  # Add facet if applicable
  if (!is.null(Facet_var)) {
    facet <- facet_wrap(vars({{Facet_var}}), scales = facet_scales)
  } else {facet <- NULL}
  
  # Generate plot
  if (density == FALSE) {
    p <- ggplot(data, aes({{X_var}}, {{Y_var}}, color = {{Color_var}})) # for geom_point
  } else {p <- ggplot(data, aes({{X_var}}, {{Y_var}}))} # for geom_pointdensity

  p <- p +
    # Scatter plot
    do.call(ifelse(density == FALSE, 'geom_point', 'geom_pointdensity'), point_params) +
    theme_publ() +
    facet
  
  
  if (density) {
    p <- p +
      scale_color_gradientn(
        colors = desaturate(lighten(rev(brewer.pal(name = 'Spectral', n = 11)), 0.2), 0.05),
        guide = 'none') 
  }
  
  return(p)
  
  
}

ggplot.boxplot <- function(data, Group_var, Value_var, metadata = NULL, filename = NULL,
                           path = NULL, title = NULL, width = 5.5, height = 5.5, ...) {
  
  library(ggplot2)
  
  Group_var_str <- rlang::as_label(ensym(Group_var))
  Value_var_str <- rlang::as_label(ensym(Value_var))
  
  # Get group variable from metadata if not in data
  if (!Group_var_str %in% colnames(data) & !is.null(metadata)) {
    data$Group <- metadata[rownames(data), which(colnames(metadata) == Group_var_str)]
  }
  # Get group variable from data
  else if (Group_var_str %in% colnames(data)) {
    data$Group <- data[,which(colnames(data) == Group_var_str)]}
  else {
    stop('Group_var not found in data or metadata')
  }
  
  p <- ggplot(data, aes(x = {{Group_var}}, y = as.numeric({{Value_var}}), fill = {{Group_var}}))  + 
    geom_boxplot(...) +
    theme_publ() +
    labs(y = Value_var_str, title = title)
  
  if (!is.null(filename)) {
    ggsave(p, 
           filename = filename, 
           path = path, 
           width = width, 
           height = height, 
           units = 'cm')
  }
  return(p)
}

ggplot.jitterbox <- function(data, Group_var, Value_var, Color_var = NULL,
                             metadata = NULL, filename = NULL, path = NULL, 
                             title = NULL, width = 5.5, height = 5.5, ...) {
  
  library(ggplot2)
  
  Group_var_str <- rlang::as_name(ensym(Group_var))
  Value_var_str <- rlang::as_name(ensym(Value_var))
  
  # Get group variable from metadata if not in data
  if (!Group_var_str %in% colnames(data) & !is.null(metadata)) {
    data$Group <- metadata[rownames(data), which(colnames(metadata) == Group_var_str)]
  }
  # Get group variable from data
  else if (Group_var_str %in% colnames(data)) {
    data$Group <- data[,which(colnames(data) == Group_var_str)]}
  else {
    stop('Group_var not found in data or metadata')
  }
  
  if (is.null(substitute(Color_var))) Color_var <- substitute(Group_var)
  Color_var_str <- rlang::as_name(ensym(Color_var))
  
  # Set default parameters
  params <- list(...)
  if (Color_var_str == Group_var_str) {
    jitter_params <- modifyList(
      list(mapping = aes(color = {{Color_var}}), alpha = 0.8, size = 0.5, width = 0.3), 
      params)
  } else {
    jitter_params <- modifyList(
      list(mapping = aes(color = {{Color_var}}), alpha = 0.8, size = 0.5, 
           position = position_jitterdodge(jitter.width = 0.3, seed = 1)), 
      params)
  }
  boxplot_params <- modifyList(
    list(mapping = aes(fill = {{Color_var}}),
         outliers = F, alpha = 0.3, lwd = 0.3), params)
  
  p <- ggplot(data, aes(x = {{Group_var}}, y = as.numeric({{Value_var}}), fill = {{Color_var}}))  + 
    do.call('geom_jitter', jitter_params) +
    do.call('geom_boxplot', boxplot_params) +
    theme_publ()
  
  if (!is.null(filename)) {
    ggsave(p, 
           filename = filename, 
           path = path, 
           width = width, 
           height = height, 
           units = 'cm')
  }
  return(p)
  
}
  

ggplot.corrplot <- function(data, X_var, Y_var, Color_var = NULL, Facet_var = NULL, 
                            cor_method = c('pearson', 'spearman'), cor_coef_digits = 2, 
                            reg_line_params = list(), facet_scales = 'fixed', 
                            xlab = NULL, ylab = NULL, title = NULL, density = FALSE,
                            label_hpos = 'left', ...) {
  
  # Load packages
  library(tidyverse)
  library(ggpp)
  library(ggpointdensity)
  
  # Load functions
  source('code/correlation_analysis_functions.R')
  
  # Get input variables
  X_var <- ensym(X_var)
  Y_var <- ensym(Y_var)
  try(Color_var <- ensym(Color_var), silent = TRUE)
  try(Facet_var <- ensym(Facet_var), silent = TRUE)
  
  # Set default parameters
  params <- list(...)
  point_params <- modifyList(list(alpha = 0.6), params)
  point_params$size <- coalesce(point_params$size, 0.5)
  reg_line_params$color <- coalesce(reg_line_params$color, 'steelblue')
  reg_line_params$lwd <- coalesce(reg_line_params$lwd, 0.5)
  xlab <- coalesce(xlab, as_label(X_var))
  ylab <- coalesce(ylab, as_label(Y_var))
  title <- coalesce(title, str_glue("{X_var}-{Y_var} correlation"))
  
  # Create correlation labels to add to the plot
  n_methods <- length(cor_method)
  y_positions <- seq(0.95, 0.95 - 0.05 * (n_methods - 1), by = -0.05)
  
  correlation_labels <- NULL
  
  for (i in seq_along(cor_method)) {
    
    m <- cor_method[i]
    
    correlation_labels[[i]] <- data %>%
      group_by({{Facet_var}}) %>% # Group data if desired
      reframe(enframe(
        my.cor.test({{X_var}}, {{Y_var}}, # Calculate correlation
                    cor_method = m)[c('Correlation', 'P.value')],
        name = 'Variable', value = 'Value')) %>%
      mutate(Value = as.numeric(Value),
             Value = ifelse(Variable == 'Correlation', round(Value, cor_coef_digits), Value),
             Value = ifelse(Variable == 'P.value', round(Value, 3), Value)) %>%
      pivot_wider(id_cols = {{Facet_var}},
                  names_from = Variable,
                  values_from = Value) %>% 
      group_by({{Facet_var}}) %>% 
      mutate(Label = case_when( # Create correlation labels
        m == 'pearson' & P.value < 0.001 ~ str_glue("Pearson's R = {Correlation}, p < 0.001"),
        m == 'pearson' & P.value >= 0.001 ~ str_glue("Pearson's R = {Correlation}, p = {P.value}"),
        m == 'spearman' & P.value < 0.001 ~ str_glue("Spearman's \u03C1 = {Correlation}, p < 0.001"),
        m == 'spearman' & P.value >= 0.001 ~ str_glue("Spearman's \u03C1 = {Correlation}, p = {P.value}")),
        npcx = label_hpos,
        npcy = y_positions[i]) %>% 
      dplyr::select({{Facet_var}}, Label, npcx, npcy)
  }
  
  correlation_labels <- bind_rows(correlation_labels)
    
  # Add facet if applicable
  if (!is.null(Facet_var)) {
    facet <- facet_wrap(vars({{Facet_var}}), scales = facet_scales)
  } else {facet <- NULL}
  
  # Generate plot
  if (density == FALSE) {
    p <- ggplot(data, aes({{X_var}}, {{Y_var}}, color = {{Color_var}})) # for geom_point
  } else {p <- ggplot(data, aes({{X_var}}, {{Y_var}}))} # for geom_pointdensity

  p <- p +
    # Scatter plot
    do.call(ifelse(density == FALSE, 'geom_point', 'geom_pointdensity'), point_params) +
    # Add correlation coefficient
    geom_text_npc(data = correlation_labels, mapping = aes(npcx = npcx, npcy = npcy, label = Label),
                  lineheight = 1, size = 6/.pt) +
    # Add regression line
    do.call('geom_smooth', reg_line_params) +
    scale_y_continuous(expand = expansion(c(0.05, 0.25))) +
    labs(x = xlab, y = ylab, title = title) +
    theme_publ() +
    facet
  
  if (density) {
    p <- p +
      scale_color_gradientn(
        colors = desaturate(lighten(rev(brewer.pal(name = 'Spectral', n = 11)), 0.2), 0.05),
        guide = 'none') 
  }
  
  return(p)
  
}

ggplot.volcano <- function(res, diff_labels = c('Down', 'Up', 'No difference'),
                           title = 'Volcano plot', label_points = TRUE, point_size = 3,
                           label_size = 6/.pt, p_threshold = 0.05, path = NULL, 
                           filename = NULL, width = 8, height = 8) {
  # volcano plot function for differential abundance analysis results
  
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  
  res$sig <- case_when(res$adj.p < p_threshold & res$log2FC > 0 ~ diff_labels[2],
                       res$adj.p < p_threshold & res$log2FC < 0 ~ diff_labels[1],
                       TRUE ~ 'No difference')
  
  p <- ggplot() +
    geom_point(data = res, 
               aes(x = log2FC, y = -log10(adj.p), color = sig),
               alpha = 0.9, size = point_size) +
    geom_hline(yintercept = -log10(0.05), linetype = 'dashed', color = 'grey25',
               lwd = 0.3) +
    geom_vline(xintercept = 0, linetype = 'dashed', color = 'grey25',
               lwd = 0.3) +
    labs(x = 'log2-fold change',
         y = '-log10(FDR)',
         title = title,
         color = 'Differential abundance') +
    theme_publ() + 
    theme(legend.position = 'bottom') +
    scale_y_continuous(expand = expansion(c(0,0.01)))
  
  if (label_points == TRUE) {
    
    p <- p + geom_text_repel(data = filter(res, adj.p < p_threshold),
                             aes(x = log2FC, y = -log10(adj.p), label = variable),
                             max.overlaps = Inf,
                             size = label_size)
  }
  
  if (!is.null(filename)) {
    
    ggsave(plot = p, 
           path = path,
           filename = filename,
           width = width, 
           height = height,
           units = 'cm')
    
  }
  
  
  return(p)
}


#### Miscellaneous #####

calculate.label.npcy <- function(y, group_n, facet_n, gap = 0.05) {
  # Function for calculating y position of labels added to plots using geom_text_npc
  # when coloring by a grouping variable and/or faceting
  
  if (y == 'top') {
    y_pos <- seq(0.95, 0.95-(group_n-1)*gap, -gap)
  }
  else if (y == 'center') {
    y_pos <- seq(0.45, 0.45-(group_n-1)*gap, -gap) + (gap/2)*(group_n-1)
  }
  else if (y == 'bottom') {
    y_pos <- seq(0.05+(group_n-1)*gap, 0.05, -gap)
  }
  else if(is.numeric(y)) {
    y_pos <- seq(y, y-(group_n-1)*gap, -gap)
  }
  
  y_pos <- rep(y_pos, each = facet_n)
  
  return(y_pos)
}

small_logticks <- function(sides) {
  # Function for adding axis ticks to log10 axes
  
  annotation_logticks(sides = sides,
                      short = unit(0.075, "cm"),
                      mid = unit(0.125, "cm"),
                      long = unit(0.2, "cm"),
                      colour = "black",
                      size = 0.2,
                      linewidth = 0.3)
  
  
}

log10_format <- function(x) {
  # Function for formatting log10 axis labels
  
  library(scales)
  
  trans_format('log10', math_format(10^.x))
  
}


format_pval_label <- function(p) {
  if (p >= 0.001) {
    paste0("p == ", sprintf("%.3f", p))
  } else {
    paste0("p == ", gsub("e", " %*% 10^", sprintf("%.2e", p)))
  }
}

##### Project-specific plots #####

cor_hist <- function(data, table_y, table_x = -0.75, median_y = 0.55, 
                     perc_colname = 'Percentage')  {
  # MS-Olink correlation histogram colored by correlation category
  
  cor_pal <- brewer.pal(4, 'Set2')
  cor_pal <- c('grey', blue_palette[c(5,3,1)])
  names(cor_pal) <- c('No correlation', 'Weak correlation', 
                      'Moderate correlation', 'Strong correlation')
  
  plot_table <- data |> 
    filter(!is.na(Correlation)) |> 
    count(Cor.interval, name = 'N') |> 
    mutate({{perc_colname}} := round(N/sum(N)*100, 1)) |> 
    dplyr::rename(Correlation = Cor.interval)
  
  median_df <- data.frame(x = 'left', y = median_y,
                          Median = median(data$Correlation, na.rm = T),
                          IQR.lower = quantile(data$Correlation, na.rm = T)[[2]],
                          IQR.upper = quantile(data$Correlation, na.rm = T)[[4]])
  data |> 
    filter(!is.na(Correlation)) |> 
    ggplot(aes(x = Correlation, fill = Cor.interval)) +
    geom_histogram(binwidth = 0.05, alpha = 0.9, boundary = 0, size = 0.3, color = 'grey20',
                   closed = 'left', pad = TRUE) +
    geom_vline(aes(xintercept = median(Correlation)), 
               lwd = 0.3, color = 'grey20', linetype = 'dashed') +
    geom_text_npc(data = median_df,
                  aes(npcx = x, npcy = y, 
                      label = paste0('Median = ', round(Median, 2), ' (IQR ', 
                                     round(IQR.lower, 2), '-', round(IQR.upper, 2), ')')),
                  size = 5/.pt) +
    labs(x = 'MS-Olink correlation', y = 'Count', fill = 'Category', color = 'Category') +
    scale_y_continuous(expand = expansion(c(0,0.05))) +
    xlim(-0.75,1) +
    theme_publ() +
    ggpp::annotate('table', label = plot_table, x = table_x, y = table_y, 
                   size = 5/.pt) +
    scale_fill_manual(values = cor_pal) +
    guides(fill = guide_legend(nrow = 2))
  
}

DAA_volcano <- function(DAA_res, labels, filename = NULL, path = NULL) {
  # Function for creating volcano plots for DAA between females and males
  
  p <- ggplot.volcano(DAA_res, 
                      diff_labels = c('Down in female', 'Up in female', 'No difference'),
                      label_points = FALSE,
                      point_size = 1,
                      title = '') +
    scale_color_manual(
      values = setNames(c(TDP_palette[c(10, 13)], 'gray'), 
                        c('Down in female', 'Up in female', 'No difference'))
    ) +
    geom_text_repel(
      data = labels, aes(x = log2FC, y = -log10(adj.p), label = Gene.Name),
      size = 6/.pt, min.segment.length = 0.2, 
      segment.size = unit(0.3, 'points'), 
      segment.color = 'grey40') +
    theme_publ()
  
  if (!is.null(filename) & !is.null(path)) {
    ggsave(p,
           filename = filename,
           path = path,
           width = 8.7,
           height = 7,
           units = 'cm')
  }
  
  return(p)
  
}
