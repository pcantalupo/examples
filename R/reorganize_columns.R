# use tidyverse
library(tidyverse)

# read_csv function from tidyverse
data = read_csv("data.csv")

# move last column to the first column
data %>%
  select(last_column_name, everything())

