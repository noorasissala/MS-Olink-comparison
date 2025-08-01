
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
  # All numeric columns present in both x and y are used to calculate the correlation.
  # match_by: How should the rows in the two data frames be matched and ordered? 
  #   Specify either the name of a column found in both data frames, 
  #   a vector of two column names where the first one is found in x and 
  #   the second one in y, or 'rownames' to use the row names.
  # cor_type: 'total' for a single correlation value for the whole data frame,
  #   'rowwise' for a correlation value for each row, 
  #   or 'colwise' for a correlation value for each column.
  # cor_method: The method to use for the correlation test, e.g. 'pearson', 
  #   'spearman', or 'kendall'.
  # cols_to_keep: A vector of additional column names to keep in the output when 
  #   cor_type is 'rowwise'.
  
  library(tidyverse) 
  
  cor_type <- match.arg(cor_type)
  
  x <- as.data.frame(x)
  y <- as.data.frame(y)
  
  if (!is.null(cols_to_keep)) {
    stopifnot('cols_to_keep should only be used with rowwise correlations' = cor_type == 'rowwise',
              'cols_to_keep should not contain the same columns as match_by' = !any(cols_to_keep %in% match_by),
              'cols_to_keep columns not found in data' = all(cols_to_keep %in% colnames(x) | cols_to_keep %in% colnames(y)))
  }
  
  # Match rows based on the specified match_by argument
  if (length(match_by) == 1 && match_by == 'rownames') {
    
    message('Matching rows by row names')
    
    id_col <- '.rowname'
    joined <- x |> 
      rownames_to_column(var = id_col) |>
      inner_join(y |> rownames_to_column(var = id_col), 
                 by = id_col,
                 suffix = c('.x', '.y'))
    
  } else if (length(match_by) == 1) {
    
    stopifnot('Invalid match_by column' = match_by %in% colnames(x) & match_by %in% colnames(y),
              'match_by column should be a character or factor' = is.character(x[[match_by]]) | is.factor(x[[match_by]]),
              'match_by column should be a character or factor' = is.character(y[[match_by]]) | is.factor(y[[match_by]]))
    
    message('Matching rows by shared column: ', match_by)
    
    id_col <- match_by
    joined <- inner_join(x, y, by = id_col, suffix = c('.x', '.y'))
    
    
  } else if (length(match_by) == 2) {
    
    stopifnot('Invalid match_by column(s)' = match_by[1] %in% colnames(x) & match_by[2] %in% colnames(y),
              'match_by column should be a character or factor' = is.character(x[[match_by[1]]]) | is.factor(x[[match_by[1]]]),
              'match_by column should be a character or factor' = is.character(y[[match_by[2]]]) | is.factor(y[[match_by[2]]]))
    
    message('Matching ', match_by[1], ' in x to ', match_by[2], ' in y')
    
    id_col <- match_by[1]
    joined <- inner_join(x, y, by = c(setNames(match_by[2], match_by[1])), suffix = c('.x', '.y'))
    
  } else {
    stop('match_by must be "rownames", a single column name, or a pair of column names.')
  }
  
  # Extract numeric columns for correlation calculation
  x_num <- joined |> 
    dplyr::select(ends_with('.x')) |>
    dplyr::select(where(is.numeric)) |> 
    dplyr::rename_with(~str_remove(., '\\.x$')) # Remove suffix from column names
  
  y_num <- joined |>
    dplyr::select(ends_with('.y')) |>
    dplyr::select(where(is.numeric)) |>
    dplyr::rename_with(~str_remove(., '\\.y$')) # Remove suffix from column names
  
  # Select common numeric columns for correlation calculation
  common_cols <- intersect(colnames(x_num), colnames(y_num))
  if (length(common_cols) == 0) {
    stop('No common numeric columns found in the two data frames.')
  }
  
  x_num <- x_num[, common_cols, drop = FALSE]
  y_num <- y_num[, common_cols, drop = FALSE]
  
  
  # Calculate correlations based on the specified cor_type
  if (cor_type == 'total') {

    return(my.cor.test(unlist(x_num), unlist(y_num), cor_method = cor_method))
    
  } else if (cor_type == 'rowwise') {
    
    res <- lapply(seq_len(nrow(x_num)), function(i) {
      my.cor.test(as.numeric(x_num[i, ]), as.numeric(y_num[i, ]), cor_method = cor_method)}) |>
      bind_rows() |>
      mutate(!!id_col := joined[[id_col]], .before = 1) |>
      mutate(Adj.p.value = p.adjust(P.value, method = 'fdr'),
             .after = 'P.value')
    
    if (!is.null(cols_to_keep)) {
      
      # check if cols_to_keep was in both data frames
      new_cols <- map(cols_to_keep, function(i) {
        if (i %in% colnames(x) & i %in% colnames(y)) {
          c(paste0(i, '.x'), paste0(i, '.y'))
        } else {i}
        }) |> 
        unlist()
        
      
      res <- bind_cols(res, 
        joined |> dplyr::select(all_of(new_cols))) |> 
        relocate(all_of(new_cols), .after = id_col)
      }
    
  } else if (cor_type == 'colwise') {
    
    res <- lapply(seq_along(common_cols), function(i) {
      my.cor.test(x_num[[i]], y_num[[i]], cor_method = cor_method)}) |> 
      bind_rows() |> 
      mutate(Sample.ID = common_cols, .before = 1) |>
      mutate(Adj.p.value = p.adjust(P.value, method = 'fdr'),
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


test_associations <- function(data, independent_vars, dependent_var) {
  
  library(dplyr)
  library(purrr)
  library(broom)
  library(stringr)
  library(tidyr)

  numeric_vars <- independent_vars[sapply(independent_vars, function(var) is.numeric(data[[var]]))]
  categorical_vars <- independent_vars[!independent_vars %in% numeric_vars]
  
  # Create one-vs-rest dummy variables for categorical predictors
  one_vs_rest_vars <- map_dfc(categorical_vars, function(var) {
    levels <- na.omit(unique(data[[var]]))
    map_dfc(levels, function(lvl) {
      col_name <- paste0(var, "_", lvl)
      tibble(!!col_name := as.integer(data[[var]] == lvl))
    })
  })
  
  colnames(one_vs_rest_vars) <- make.names(colnames(one_vs_rest_vars))
  
  # Z-score transform numeric variables (independent and dependent)
  data <- data %>%
    mutate(across(all_of(c(numeric_vars, dependent_var)), ~ scale(.) %>% as.numeric())) |> 
    bind_cols(one_vs_rest_vars) |> 
    dplyr::select(-all_of(categorical_vars))  # Remove original categorical variables

  vars <- c(numeric_vars, colnames(one_vs_rest_vars))
  
  # Train lm models for each independent variable
  lm_results <- tibble(Variable = vars) %>%
    mutate(
      Model = map(Variable, ~ lm(as.formula(paste(dependent_var, '~', .x)), data = data)),
      Tidy = map2(Model, Variable, ~ {
        ci <- confint(.x) # Get confidence intervals
        tidy(.x) %>%
          filter(term != '(Intercept)') %>%
          mutate(
            term = if (is.numeric(data[[.y]])) .y else str_remove(term, paste0('^', .y)),
            Conf.low = ci[term, 1],
            Conf.high = ci[term, 2])
        }),
      R.squared = map_dbl(Model, ~ glance(.x)$r.squared),
      Adj.r.squared = map_dbl(Model, ~ glance(.x)$adj.r.squared)) |> 
    unnest(Tidy) |> 
    dplyr::rename(Coef = estimate) |> 
    dplyr::rename_with(str_to_sentence) |> 
    dplyr::select(Variable, Coef, Conf.low, Conf.high, everything(), -Model, -Term) |> 
    mutate(Adj.p.value = p.adjust(P.value, method = 'fdr'), .after = 'P.value') |> 
    arrange(Adj.p.value)
  
}
