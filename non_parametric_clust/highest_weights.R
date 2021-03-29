#Script used to generate the best k9 clusters determined by computer farm

library('flexclust') # handles weighted clustering

setwd("~/Documents/GitHub/BioinfGroupProject/non_parametric_clust")

#set.seed(255)
#read in and format counts (z_scores)
z_scores = read.csv("data_expression_median.txt", sep = "\t", stringsAsFactors = FALSE, header = TRUE)

map = z_scores[,c(1,2)]

inst_names = colnames(z_scores)[3:length(colnames(z_scores))]

z_scores = z_scores[,-2]
z_scores = as.data.frame(t(z_scores), stringsAsFactors = FALSE)
names(z_scores) = map[,1]
z_scores = z_scores[-1,]

z_scores = as.data.frame(apply(z_scores, 2, as.numeric))

name_fix = c()
for(name in inst_names){name_fix = c(name_fix,substr(name,11,nchar(name)-3))}


row.names(z_scores) = name_fix

#read in ranked file, extract names

ranking = read.csv('ranked_gene_list.tsv', sep = "\t")
rank_scores = ranking$SCORE

ranking = ranking$NAME


rz = z_scores[,colnames(z_scores) %in% ranking] #get genes from count included in ranking

rk = ranking[ranking %in% colnames(z_scores)] # get ranked genes included in counts
rk_score = rank_scores[ranking %in% rk]

rz_order = rz[,rk] # order columns by ranking
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
wg = range01(abs(rk_score))

plot(1:length(rk_score), wg, type = 'l', xlab = "Rank in Gene List", ylab = "Weighting", main = "Distribution of Weights Across Ranked Gene List")




library('Rfast')


wgs = data.frame()
for (i in 1:20){
  hi_wg = nth(wg, i, descending = T)
  gene_wg = ranking[nth(wg, i, descending = T, index.return = T)]
  df = data.frame(gene_wg, hi_wg)
  wgs = rbind(wgs, df)
}

colnames(wgs) = c("Gene", "Weighting")

write.table(wgs, file = "Top_20_weighted.tsv", quote = F, row.names = F, sep = "\t")

