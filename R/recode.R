# Recoding recode example

library(sjmisc)
mtcars %>% select(gear, carb) %>%
  rec(rec = "min:3=1; 4:max=2")

