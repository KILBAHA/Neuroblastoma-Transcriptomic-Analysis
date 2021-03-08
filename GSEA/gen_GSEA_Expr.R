setwd("~/Documents/GitHub/BioinfGroupProject/GSEA")


### Read in and Format DF ###
count = read.csv("data_expression_median.txt", sep = "\t")

gene_names = count[,2]

count = count[,3:ncol(count)]
row.names(count) = gene_names

inst_names = colnames(count)

name_fix = c()
for(name in inst_names){name_fix = c(name_fix,substr(name,11,nchar(name)-3))}

colnames(count) = name_fix

### Find missing values ###
miss = count[rowSums(is.na(count)) > 0,] #no missing variables!

### Remove values with no phenotype (run PheLab first!) ### 
#todel = c("PACPJG", "PAJXLE") #change these variables to the output of nodat in PheLab file

#count = count[,-which(colnames(count) %in% todel)]


### Format into output ### 

#Number of genes
n_genes = nrow(count)

#Number of samples 
n_samp = ncol(count)

#NAME and Description
NAME = gene_names
rownames(count) =c()
description = c(rep(NA,length(NAME)))

count = cbind(NAME,description, count)

app = cat("#1.2","\n",n_genes,"\t",n_samp, "\n",file='expr2.gct')

write.table(count,file = 'expr2.gct', sep='\t', append= T, row.names = F)
