setwd("~/Documents/GitHub/BioinfGroupProject/GSEA_MARKER")

### Read in z_scores and format ###
z_scores = read.csv("data_mRNA_median_all_sample_Zscores.txt", sep = "\t", stringsAsFactors = FALSE, header = TRUE)

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


### Function that splits data up into high and low expression for given gene. Returns names of individuals###
get_split = function(gene, threshold = 0.5, data = z_scores){
  select = z_scores[,gene]
  above = row.names(data)[which(select > threshold)]
  below = row.names(data)[which(select <= threshold)]
  return(list(above,below)) #writes to list. To index: > threshold = list[1], <= threshold = list[2]
}

marker_genes = c("TP53", "PTPN6", "NTRK1", "MYCN","IFITM3") #list of genes we wish to split the data by

marker_genes %in% colnames(z_scores) #check all genes are present in z_scores - used to check for differences between naming conventions

#write splits of all markers to files in markers folder
for(marker in marker_genes){
  spl = get_split(marker)
  fname = paste("markers/", marker, sep = "") #we want to put this in a separate folder called markers hence the paste
  write.table(spl[1],file = paste(fname, "_above.txt", sep="") , sep = "\n", row.names = F, col.names = F, quote = T)
  write.table(spl[2],file = paste(fname, "_below.txt", sep="") , sep = "\n", row.names = F, col.names = F, quote = T)
}