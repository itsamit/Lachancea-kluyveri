library(dplyr)
library(plyr)
library(tidyselect)

setwd('D:/ScientificReports/GeneExpressionIntegration/')

expression_array=read.csv('GSE48135.txt',sep='\t')
idmap=read.csv('UniversalIDmap.csv')

colnames(expression_array)[1]='ID'

expr=join(idmap,expression_array,type='left')

boxplot(expr[,5:20])

summary(expr[,5:20])

df = subset(expr, select = -c(X,YeastID))

df_scaled=as.matrix(scale(df[,3:18]))

rownames(df_scaled)=df$ModelID

write.csv(df_scaled,'GeneExpressionNitrogenSources.csv')