
rowwise.t.test <- function(data, test_grp, ref_grp) {
  # Performs an independent t-test for each row of a data frame
  # between columns specified in test_grp and ref_grp
  # test_grp is compared against ref_grp (i.e. test_grp - ref_grp)
  
  library(broom)
  library(dplyr)
  
  res <- apply(data, 1, 
               function(x, test_grp, ref_grp) {
                 tryCatch({
                   tidy(t.test(as.numeric(x[test_grp]), as.numeric(x[ref_grp]), 
                               na.action = 'na.omit'))},
                   error = function(e) {
                     as_tibble(as.list(
                       setNames(rep(NA, 10), c('estimate', 'estimate1', 'estimate2', 'statistic', 
                                               'p.value', 'parameter', 'conf.low', 'conf.high', 
                                               'method', 'alternative'))))}
                   )}, 
               test_grp = test_grp, ref_grp = ref_grp)
  
  n_df <- apply(data, 1, 
                 function(x, test_grp, ref_grp) {
                   n_test <- sum(!is.na(x[test_grp]))
                   n_ref <- sum(!is.na(x[ref_grp]))
                   c(n.testGrp = n_test, n.refGrp = n_ref, n.total = n_test + n_ref)
                 }, 
                 test_grp = test_grp, ref_grp = ref_grp) |> 
    t() |> 
    as.data.frame()
  
  res <- bind_rows(res) |> 
    mutate(adj.p = p.adjust(p.value, method = 'fdr'),
           variable = rownames(data)) |>
    dplyr::select(variable, estimate:p.value, adj.p, parameter:alternative) |> 
    dplyr::rename(log2FC = estimate,
                  mean.testGrp = estimate1,
                  mean.refGrp = estimate2) |> 
    bind_cols(n_df) |> 
    relocate(n.testGrp, n.refGrp, n.total, .after = mean.refGrp)
  
  return(res)
  
}

