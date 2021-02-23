setwd("~/Documents/GitHub/BioinfGroupProject/GSEA/extract_genes")

GE_output = read.csv("HALLMARK_MYC_TARGETS_V2.tsv", sep = "\t")

GE_output = GE_output[,c("SYMBOL", "CORE.ENRICHMENT")]

core = GE_output[GE_output$"CORE.ENRICHMENT"== "Yes",]

count = read.csv("data_expression_median.txt", sep = "\t")


gene_names = count[,1]

length(gene_names) #check no duplicate entries
length(unique(gene_names)) 



count = count[,3:ncol(count)]
row.names(count) = gene_names

inst_names = colnames(count)

name_fix = c()
for(name in inst_names){name_fix = c(name_fix,substr(name,11,nchar(name)-3))}

colnames(count) = name_fix

select = count[gene_names %in% core$SYMBOL,]

write.csv(select, file="GE_selected.csv")
