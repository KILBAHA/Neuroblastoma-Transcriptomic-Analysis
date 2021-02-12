pos = read.table("pos.txt")
for (p in pos$V1){getset(p)}

neg = read.table("neg.txt")
for(n in neg$V1){getset(n)}


getset = function(GSEA_output_filename, only_core=TRUE, count_filename="data_expression_median.txt"){
  GE_output = read.csv(GSEA_output_filename, sep = "\t")
  
  GE_output = GE_output[,c("SYMBOL", "CORE.ENRICHMENT")]
  
  if (only_core == TRUE){core = GE_output[GE_output$"CORE.ENRICHMENT"== "Yes",]}
  
  count = read.csv("data_expression_median.txt", sep = "\t")
  
  gene_names = count[,1]
  
  length(gene_names) == length(unique(gene_names))
  
  count = count[,3:ncol(count)]
  row.names(count) = gene_names
  
  inst_names = colnames(count)
  
  name_fix = c()
  for(name in inst_names){name_fix = c(name_fix,substr(name,11,nchar(name)-3))}
  
  colnames(count) = name_fix
  
  select = count[gene_names %in% core$SYMBOL,]
  
  
  
  patient = as.data.frame(read.csv("data_clinical_patient.csv", sep = "\t", stringsAsFactors = FALSE, header = TRUE, skip = 4))
  p_nam = patient[,1]
  
  p_nam_fix = c()
  for(name in p_nam){p_nam_fix = c(p_nam_fix,substr(name,11,nchar(name)))}
  
  patient = patient[,2:ncol(patient)]
  
  row.names(patient) = p_nam_fix
  
  five_yr = 365*5
  
  patient$EFS_fct = ifelse(patient[,"EFS_TIME"] > five_yr, "EFS_G5", "EFS_L5")
  
  efs = data.frame(patient[,"EFS_fct"])
  row.names(efs) = p_nam_fix
  select = t(select)
  
  select = merge.data.frame(select,efs, by="row.names")
  
  names(select)[length(names(select))] = "EFS_fct"
  
  select = select[!(is.na(select$EFS_fct)),]
  
  write_filename = paste("redfiles/",GSEA_output_filename, "_dimred.csv", sep="")
  
  write.csv(select, write_filename)
}


GSEA_output_filename = "ASGHARZADEH_NEUROBLASTOMA_POOR_SURVIVAL_DN.tsv"
count_filename = "data_expression_median.txt"
output_count_reduced_filename = "nb_poor_survival.csv"
only_core = TRUE #set to false if you want to include all genes in set and not just core enrichment


GE_output = read.csv(GSEA_output_filename, sep = "\t")

GE_output = GE_output[,c("SYMBOL", "CORE.ENRICHMENT")]

if (only_core == TRUE){core = GE_output[GE_output$"CORE.ENRICHMENT"== "Yes",]}

count = read.csv("data_expression_median.txt", sep = "\t")

gene_names = count[,1]

length(gene_names) == length(unique(gene_names))

count = count[,3:ncol(count)]
row.names(count) = gene_names

inst_names = colnames(count)

name_fix = c()
for(name in inst_names){name_fix = c(name_fix,substr(name,11,nchar(name)-3))}

colnames(count) = name_fix

select = count[gene_names %in% core$SYMBOL,]



patient = as.data.frame(read.csv("data_clinical_patient.csv", sep = "\t", stringsAsFactors = FALSE, header = TRUE, skip = 4))
p_nam = patient[,1]

p_nam_fix = c()
for(name in p_nam){p_nam_fix = c(p_nam_fix,substr(name,11,nchar(name)))}

patient = patient[,2:ncol(patient)]

row.names(patient) = p_nam_fix

five_yr = 365*5

patient$EFS_fct = ifelse(patient[,"EFS_TIME"] > five_yr, "EFS_G5", "EFS_L5")

efs = data.frame(patient[,"EFS_fct"])
row.names(efs) = p_nam_fix
select = t(select)

select = merge.data.frame(select,efs, by="row.names")

names(select)[length(names(select))] = "EFS_fct"

select = select[!(is.na(select$EFS_fct)),]
  

