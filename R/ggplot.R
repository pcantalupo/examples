# Use 'gridExtra::grid.arrange' to place multiple ggplots in a plot
library(ggplot2)
library(gridExtra)
plots = lapply(names(mock), function(n) {
  g = TSNEPlot(mock[[n]], do.label = T, do.return = T)
  g + ggtitle(paste0("Mock.",n))
})
grid.arrange(grobs = plots, ncol = 2)


ggplot()
