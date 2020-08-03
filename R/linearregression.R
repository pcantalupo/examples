# Linear regression example
# Useful links https://rstudio-pubs-static.s3.amazonaws.com/65641_88a692252c6c4f2ab279d115e59e6767.html
#              https://www.youtube.com/watch?v=ZkjP5RJLQF4&list=PLIeGtxpvyG-LoKUpV0fSY8BGKIMIdmfCi

# Can think of weights as being gene expression values and height being number of UMI in the cells
weight = c(137, 130, 154, 155, 160, 190)
height = c(60, 65, 75, 80, 85, 90)

df = data.frame(w = weight, h = height)
model = lm(weight~height, data=df)

# the plot shows that the weight of 137 is higher than predicted by its height
plot(height,weight)
abline(model)

# Plot the Model
layout(matrix(c(1,2,3,4),2,2))
plot(model, main = "Height by Weight")

# The weight 137 person has a high residual (8.59)
resid(model)
plot(resid(model))
plot(model, which = 1)

# After scaling the residuals, this shows that weight 137 (low weight or low expression) and weight 190 (high weight or expression) have high weight (expression) compared to everybody else. This is because we regressed out the effect of height (UMI) on the weights (expression)
srw = scale(resid(model))
plot(srw, weight)


# Show summary of the model
msum = summary(model)
names(msum)
msum

# Get R^2 and Adj R^2
msum$r.squared
msum$adj.r.squared

# Model p-value: If you want to obtain the p-value of the overall regression model use this function
# from: https://stackoverflow.com/questions/5587676/pull-out-p-values-and-r-squared-from-a-linear-regression
lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
lmp(model)

# In the case of a simple regression with one predictor, the model p-value and the p-value for the coefficient will be the same.

# Coefficient p-values:
# If you have more than one predictor, then the above will return the model p-value, and the p-value for coefficients can be extracted using:
msum$coefficients[,4]  

# Show R^2 and Pvalue on plot
r2 = format(msum$adj.r.squared, digits = 3)
myp = format(lmp(model), format ="e", digits = 2)
plot(height, weight)
abline(model)
text(70, 170, labels = paste0("R2 = ", r2, " pval = ", myp))
dev.off()



