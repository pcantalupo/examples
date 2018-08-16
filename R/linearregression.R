# Linear regression example
# Can think of weights as being gene expression values and height being number of UMI in the cells
weight = c(137, 130, 154, 155, 160, 190)
height = c(60, 65, 75, 80, 85, 90)

df = data.frame(w = weight, h = height)
model = lm(weight~height, data=df)

# the plot shows that the weight of 137 is higher than predicted by its height
plot(height,weight)
abline(model)

# The weight 137 person has a high residual (8.59)
resid(model)
plot(resid(model))
plot(model, which = 1)

# After scaling the residuals, this shows that weight 137 (low weight or low expression) and weight 190 (high weight or expression) have high weight (expression) compared to everybody else. This is because we regressed out the effect of height (UMI) on the weights (expression)
srw = scale(resid(model))
plot(srw, weight)
