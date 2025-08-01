
CV.log2 <- function(x) {
  
  sqrt(exp((log(2)*sd(x, na.rm = T))^2)-1)*100
  
}

calc_tech_CV <- function(ms_data, olink_data) {
  
  olink_CV <- olink_data %>% 
    group_by(UniProt, OlinkID) %>% 
    filter(str_detect(Sample.ID, 'Control')) %>% 
    summarize(Technical.CV = CV.log2(NPX)) |> 
    mutate(Analyzed = !is.na(Technical.CV) & !is.nan(Technical.CV),
           Technical.CV = ifelse(Technical.CV > 100, 100, Technical.CV), # cap CV at 100%
           Technical.CV = ifelse(is.nan(Technical.CV), NA, Technical.CV))
  
  replicates <- c('LCP27' , 'LCP106' , 'LCP65' , 'LCP76' , 'LCP109' , 'LCP90')
  
  ms_CV <- ms_data %>% 
    dplyr::select(UniProt, contains(replicates)) %>% 
    pivot_longer(cols = -UniProt,
                 names_to = 'Sample.ID',
                 values_to = 'Value') %>% 
    mutate(Sample.ID = str_remove(Sample.ID, 'a|b')) %>% 
    group_by(UniProt, Sample.ID) %>% 
    summarize(Technical.CV = CV.log2(Value)) %>% 
    group_by(UniProt) %>% 
    summarize(Technical.CV = mean(Technical.CV, na.rm = T)) |> 
    mutate(Analyzed = !is.na(Technical.CV) & !is.nan(Technical.CV),
           Technical.CV = ifelse(Technical.CV > 100, 100, Technical.CV), # cap CV at 100%
           Technical.CV = ifelse(is.nan(Technical.CV), NA, Technical.CV))
  
  CV <- bind_rows('MS' = ms_CV, 'Olink' = olink_CV, .id = 'Platform')
  
  return(CV)
  
}
