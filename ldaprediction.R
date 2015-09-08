rm(list=ls())

# 10 genes (rows) by 10 samples (columns)
d = data.frame(matrix(rnorm(100),10,10))
rownames(d) = c(LETTERS[1:10])

# 4 expt, 3 control, 3 other
groups = c(rep("E",4),rep("C",3),rep("O",3))

# for vector of values (x) and group labels for each value (g)
anovapval <- function (x, g) {
  a = anova(lm(x ~ g))
  return(a$"Pr(>F)"[1])
}

# for each gene (in rows) get anova pvalue
pvals = apply(d, 1, anovapval, groups)
pvals
signif_genes = names(head(sort(pvals),n=3))
fmla = as.formula(paste("groups ~ ", paste(signif_genes, collapse = "+")))
fmla

# transpose so genes are in columns and samples in rows
td = as.data.frame(t(d))

library(MASS)
model = lda(fmla, td)
# predict class of the first Sample with predict
predict(model,td[1,])$class

# test accuracy of prediction with CV=TRUE
model.cv = lda(fmla, td, CV=T)
class(model.cv)
model.cv$class

# Manually create a confusion matrix and get simple stats for prediction
correct = table(model.cv$class, groups)
correct
# percent ccorrect for each category
diag(prop.table(correct, 1))
# total percent correct
sum(diag(prop.table(correct)))

# CARET - create confusion matrix and many stats for prediction!
library(caret)
#creates a confusion matrix with sensitivity, specificity, accuracy,
#kappa and much more; please see documentation (?confusionMatrix) for
#details.
confusionMatrix(data=model.cv$class, reference=groups) 


# from BRB-array tools
