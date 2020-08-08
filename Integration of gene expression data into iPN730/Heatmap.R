library(ggplot2)
library(reshape2)
library(plyr)

setwd('D:/ScientificReports/GeneExpressionIntegration/')

fva=read.csv('EnrichedReactionsFluxRange.csv',header = T)
paths=read.csv('PathwaysEnrichedImp.csv',header = F)

fva=data.frame(fva)
paths=data.frame(paths)

colnames(paths)=c('Reactions','Pathway')

merged=join(fva,paths,by='Reactions',type='right')

dataA=as.matrix(merged[,2:3])
dataU=as.matrix(merged[,4:5])

rownames(dataA)=merged[,1]
rownames(dataU)=merged[,1]
colnames(dataA)=c('minFluxA','maxFluxA')
colnames(dataU)=c('minFluxU','maxFluxU')

data=cbind(dataA,dataU)

mat=data

mat=round(mat,6)

ggplot(melt(mat), aes(Var2, Var1, fill= value)) + 
  geom_tile(aes(width=0.98, height=.93)) + geom_text(aes(label=value),color='white', size=3.0)+
  theme_classic()






