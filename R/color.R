# Color fun in R
# see https://www.nceas.ucsb.edu/~frazier/RSpatialGuides/colorPaletteCheatsheet.pdf)
library(RColorBrewer)
display.brewer.all()

library(colorspace)
pal = choose_palette()
col2rgb("salmon")

# Emulate ggplot2 default color palette:
# https://stackoverflow.com/questions/8197559/emulate-ggplot2-default-color-palette
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
n=5
cols = gg_color_hue(5)
plot(1:n, pch = 16, cex = 2, col = cols)
