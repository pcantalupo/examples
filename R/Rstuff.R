#Log2 dataframe
  #set 0 values to 0.001
    r[r==0] = 0.001    # get rid of 0 values
  #after logging, change “-Inf” values to NAs and get rows with no NAs
    r[r == "-Inf"] = NA
    r = r[complete.cases(r),]

#Write table
write.table(df, “df.tsv”, sep="\t", row.names=F, quote=F, na="")

#grep
#Get values of words that start with ‘f’
x = c(“foo”,”fab”,”bar”)
grep(“^f”, x, value=T)
    foo  fab

#Get indices of words that start with ‘f’
grep(“^f”, x)
   1  2

#Box plots
library(reshape2)
m.melt = melt(rpkm[,c("cell_type", rbe2f)])
	creates df with: cell_type \t variable [gene] \t value [rpkm]
boxplot(log2(value) ~ cell_type * variable, data =m.melt, notch=T, col=(c("gold","darkgreen")), main = "lihc Rb-E2F Summary", cex.axis=0.8,las=2)

# normal vs tumor for set of genes
g = ggplot(m.melt, aes(x=cell_type, y=log2(value)))
g + geom_boxplot() + facet_grid(~ variable) + geom_jitter(alpha=I(1), aes(color=cell_type), size=3) + theme(legend.position = "none")

# if you want to color the boxplots, use aes(fill=VARIABLE) in geom_boxplot(). If you want to color the borders of a boxplot, use aes(color=VARIABLE)

#Pheatmap
annotvirus = subset(virus_noNA, select=c("hbv"))
annotclin = subset(clin_noNA, select=c("cell_type"))
a = cbind(annotvirus,annotclin)
r = rpkm_noNA[,rbe2f]
r[r==0] = 0.001    # get rid of 0 values
identical(rownames(r), rownames(a))
pheatmap(log2(t(r)), annotation=a, show_colnames=F)

#Correlations
source("~/local/lib/R/Correlator.R")
spam = merge(x=virusT_noNA, y=clin_noNA, by=0)   # analyze tumors only
c = correlator(spam,7,7,13,ncol(spam))
c = c[!is.nan(c$p.value),]
c = c[c$p.value < 0.05,]
c = c[order(c$p.value),]
saveCorrTables(c,spam)
write.table(c, "corr.tsv", quote=F, row.names=F, sep="\t")

#Dcast
# from VICALL table, show the number of species X in each aid
#  this creates a df with R rows equal to number of unique AIDs
#                     and C columns equal to number of unique species
dcast(v, aid ~ species, fun.aggregate = length, value.var=”aid”)



#Dplyr
library(dplyr)
# get virus with most number of alignments for each patient_barcode
vic_cesc_hpv_max_aligns = vicall %>%
  filter(keep == 'yes', cancer == 'CESC') %>%
  select(pbc, organism, aligns) %>%
  group_by(pbc) %>%
  filter(min_rank(desc(aligns)) <= 1)

#Other summarizations
df = aggregate(E2F1 ~ cell_type, data=mm, FUN=mean) # return Data frame with R rows equal to number of unique values of ‘cell_type’
		cell_type	E2F1
  1	N		382.2
	2	T		916.7
v = c(1,2,3,4,5,6,7,8);   i = c(rep("n",4),rep("t",4))
tapply(v, FUN=mean, INDEX=i)  # returns named vector (names = “n” and “t”) containing the means.


#Relabelling a factor
library(plyr)
# relabel Classes to shorter abbreviations
n$Class = revalue(n$Class, c("N"="N", "abwt"="abwt", "apcmut"="am","bcatmut"="bm","abmut"="abm"))

#Reordering a factor
# reorder Classes to show Normals first
n$Class = factor(n$Class, levels(n$Class)[c(5,2,3,4,1)])
table(n$Class,n$gene)


# How to merge dataframes without adding a “Row.names” column?
# http://bit.ly/2dosHll
transform(merge(df1,df2,by=0,all=TRUE), row.names=Row.names, Row.names=NULL)


# Statistical Testing
# packages: gmodels (CrossTable function; see http://stackoverflow.com/questions/32651253/how-can-i-extract-error-percentage-from-a-crosstable-into-a-variable)


