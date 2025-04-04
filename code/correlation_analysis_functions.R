
my.cor.test <- function(x, y, cor_method = 'pearson', verbose = TRUE) {
  
  # Correlation test that uses pairwise complete observations
  # and outputs the correlation, p-value, number of complete pairwise observations,
  # and the calculation method
  
  library(stringr)
  
  res <- tryCatch({cor.test(x, y, use = 'pairwise.complete.obs', method = cor_method)[c('estimate', 'p.value')]},
    error = function(cond) {
      if(verbose) message('Not enough observations, returning NA')
      return(list('estimate' = NA, 'p.value' = NA))})
  res$estimate <- ifelse(is.null(res$estimate), NA, as.vector(res$estimate))
  res$N <- sum(complete.cases(x, y))
  res$Cor.method <- str_to_sentence(cor_method)
    names(res)[1:2] <- c('Correlation', 'P.value')
  
  return(res)
  
}


cor.test.df <- function(x, y, match_by = 'rownames', cor_type = c('total', 'rowwise', 'colwise'),
                        cor_method = 'pearson', cols_to_keep = NULL) {
  
  # Function for calculating the correlation between two data frames in wide format.
  # Columns are matched by name
  # Rows are matched by row names or values of a specified column/columns.
  # All numeric columns are used to calculate the correlation.
  # match_by: How should the rows in the two data frames be matched and ordered? 
  #   Specify either the name of a column found in both data frames, 
  #   a vector of two column names where the first one is found in x and 
  #   the second one in y, or "rownames" to use the row names.
  
  library(dplyr) 
  
  cor_type <- match.arg(cor_type)
  
  x <- as.data.frame(x)
  y <- as.data.frame(y)
  
  
  if (length(match_by) == 1) {
    
    if (match_by == 'rownames') {
      
      message('Matching rows by row names')
      
      # Order rows
      x <- x[rownames(x) %in% rownames(y),]
      y <- y[rownames(y) %in% rownames(x),]
      y <- y[rownames(x),]
      
      if (cor_type == 'rowwise') id_col <- rownames(x)
      
    } else if (match_by != 'rownames') {
      
      stopifnot('Invalid match_by column' = match_by %in% colnames(x) & match_by %in% colnames(y),
                'match_by column should be a character or factor' = is.character(x[[match_by]]) | is.factor(x[[match_by]]),
                'match_by column should be a character or factor' = is.character(y[[match_by]]) | is.factor(y[[match_by]]))
      
      # Order rows
      x <- x[x[[match_by]] %in% y[[match_by]],]
      y <- y[y[[match_by]] %in% x[[match_by]],]
      y <- y[match(x[[match_by]], y[[match_by]]),]
      
      if (cor_type == 'rowwise') id_col <- x[[match_by]]
      
    }
    
  } else if (length(match_by) == 2) {
    
    stopifnot('Invalid match_by column(s)' = match_by[1] %in% colnames(x) & match_by[2] %in% colnames(y),
              'match_by column should be a character or factor' = is.character(x[[match_by[1]]]) | is.factor(x[[match_by[1]]]),
              'match_by column should be a character or factor' = is.character(y[[match_by[2]]]) | is.factor(y[[match_by[2]]]))
    
    # Order rows
    x <- x[x[[match_by[1]]] %in% y[[match_by[2]]],]
    y <- y[y[[match_by[2]]] %in% x[[match_by[1]]],]
    y <- y[match(x[[match_by[1]]], y[[match_by[2]]]),]
    
    if (cor_type == 'rowwise') id_col <- x[[match_by[1]]]
    
  }
  
  # Order columns
  if (!is.null(cols_to_keep)) {
    stopifnot('cols_to_keep columns not found in data' = all(cols_to_keep %in% colnames(x) | cols_to_keep %in% colnames(y)),
              'cols_to_keep should not contain the same columns as match_by' = !any(cols_to_keep %in% match_by),
              'cols_to_keep should only be used with rowwise correlations' = cor_type == 'rowwise')
    # Keep columns specified in cols_to_keep
    cols_to_keep <- bind_cols(dplyr::select(x, any_of(cols_to_keep)),
                              dplyr::select(y, any_of(cols_to_keep)),
                              .name_repair = 'universal')
  }
  x <- dplyr::select(x, where(is.numeric))
  x <- x[colnames(x) %in% colnames(y)]
  y <- y[colnames(x)]
  
  
  if (cor_type == 'total') {

    res <- my.cor.test(unlist(x), unlist(y), cor_method = cor_method)
    
    } 
  else if (cor_type == 'rowwise') {
    
    res <- lapply(1:nrow(x), function(i) my.cor.test(as.numeric(x[i,]), as.numeric(y[i,]), cor_method = cor_method))
    res <- data.frame(id_col, do.call(rbind.data.frame, res))
    colnames(res)[1] <- match_by[1]
    res <- add_column(res, 
                      Adj.p.value = p.adjust(res$P.value, method = 'fdr'),
                      .after = 'P.value')
    if (!is.null(cols_to_keep)) res <- add_column(res, cols_to_keep, .after = 1)
    
    }
  else if (cor_type == 'colwise') {
    
    res <- lapply(1:ncol(x), function(i) my.cor.test(as.numeric(x[,i]), as.numeric(y[,i]), cor_method = cor_method))
    res <- data.frame('Sample.ID' = colnames(x), do.call(rbind.data.frame, res))
    res <- add_column(res,
                      Adj.p.value = p.adjust(res$P.value, method = 'fdr'),
                      .after = 'P.value')
    }
  
  return(res)
}


cor_intervals <- function(x) {
  cut(x, breaks = c(-1, 0.3, 0.5, 0.7, 1), 
      labels = c('No correlation', 'Weak correlation', 
                 'Moderate correlation', 'Strong correlation'),
      include.lowest = TRUE, right = FALSE)
}
