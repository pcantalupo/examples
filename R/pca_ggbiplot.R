# see https://www.datacamp.com/community/tutorials/pca-analysis-r


mtcars
cars = mtcars[,c(1:7,10,11)]
pc = prcomp(cars, scale = T, center = T)
summary(pc)


mtcars.country <- c(rep("Japan", 3), rep("US",4), rep("Europe", 7),rep("US",3), "Europe", rep("Japan", 3), rep("US",4), rep("Europe", 3), "US", rep("Europe", 3))
ggbiplot(pc, ellipse=TRUE,  labels=rownames(cars), groups=mtcars.country)


