
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
  
  res <- as.data.frame(do.call(rbind.data.frame, res))
  res$adj.p <- p.adjust(res$p.value, method = 'fdr')
  res$variable <- rownames(data)
  res <- dplyr::select(res, variable, estimate:p.value, adj.p, parameter:alternative)
  colnames(res)[2:4] <- c('log2FC', 'mean.testGrp', 'mean.refGrp')
  
  return(res)
  
}