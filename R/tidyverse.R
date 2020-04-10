# separate (tidyr) is breaking up the column 'coverage' into two columns called 'covleft', 'covright' using the separator ','
# mutate (dplyr) is adding a new column that is calculated from values in other columns

igvjunctions = igvjunctions %>%
  separate(coverage, c("covleft","covright"), sep=",", extra="drop") %>%     
  mutate(firstintron = as.integer(start) + as.integer(covleft) + 1) %>%
  mutate(lastintron = as.integer(end) - as.integer(covright)) %>%
  mutate(welchintron = paste0(firstintron + 57, "-", lastintron + 57)) %>%
  select(depth, welchintron)




# How do I spread multiple columns? Need to create a new variable with `unite`
# Gather -> Unite -> Spread pipeline
https://stackoverflow.com/questions/30592094/r-spreading-multiple-columns-with-tidyr

