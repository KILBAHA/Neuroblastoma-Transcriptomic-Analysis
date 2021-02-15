setwd("~/Documents/GitHub/BioinfGroupProject/GSEA_rank_cluster")

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
ranking = ranking$NAME


rz = z_scores[,colnames(z_scores) %in% ranking] #get genes from count included in ranking

rk = ranking[ranking %in% colnames(z_scores)] # get ranked genes included in counts

rz_order = rz[,rk] # order columns by ranking


library('flexclust')

c("j","k", "a","c","b","d","e") %in% c("a","b")

cc = cclust(rz_order, method = "hardcl", k=2, weights = 1/(nrow(rz_order):1))

c1 = names(which(clusters(cc) == 1))
c2 = names(which(clusters(cc) == 2))

write.table(c1, file="c1.txt", sep = "\n", row.names = F, col.names = F, quote = F)
write.table(c2, file = "c2.txt", sep = "\n", row.names = F, col.names = F, quote = F)



              