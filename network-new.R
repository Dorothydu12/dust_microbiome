library(igraph)
library(Hmisc)
source("//Users/shicong/Desktop/R script/random network/1.Pairwise_correlations.R")
otu=read.table("SL.txt",header=T,row.names = 1)
pattern<-co_occurrence_network(otu,0.8,0.001) 
write.graph(pattern$graph3,"SL-0.8-0.001.gml",format='gml')

# 去掉列和为 0 的列
df <- df[, colSums(df) != 0]

library(WGCNA)
library(igraph)
CorrDF <- function(cormat, pmat) {
       ut <- upper.tri(cormat) # 由于相关性矩阵是对称的，取上三角矩阵计算即可
       data.frame(
             from = rownames(cormat)[col(cormat)[ut]],
             to = rownames(cormat)[row(cormat)[ut]],
             cor = (cormat)[ut],
             p = pmat[ut]
         )
}
occor <- corAndPvalue(df, use='pairwise', method='spearman')
cor_df <- CorrDF(occor$cor , occor$p)
cor_df <- cor_df[which(abs(cor_df$cor) >= 0.6),]
cor_df <- cor_df[which(cor_df$p < 0.01),]
igraph <- graph_from_data_frame(cor_df, direct=F)
igraph1=upgrade_graph(igraph)
hub=hub_score(igraph1)$vector%>%
       sort(decreasing = TRUE)%>%
   as.data.frame()
colnames(hub)="hub_sca"
ggplot(hub)+geom_bar(aes(x = hub_sca,y = reorder(row.names(hub),hub_sca)),stat = "identity",fill = "#4DAF4A")
igraph=simplify(igraph)
length(V(igraph))
length(E(igraph))
write.graph(igraph_s,"50.gml", format="gml")
