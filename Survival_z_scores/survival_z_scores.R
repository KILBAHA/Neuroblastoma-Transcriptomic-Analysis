### Read in z_scores and format ###
z_scores = read.csv("data_RNA_Seq_mRNA_median_all_sample_Zscores.csv", sep = "\t", stringsAsFactors = FALSE, header = TRUE)

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

### Deal with missing vals ###

na_count = colSums(is.na.data.frame(z_scores)) # NAs in each gene

thresh = 0 # change to vary number of allowed NAs
to_del = c()
for(name in names(z_scores)){
  if(na_count[name]>thresh){
    to_del = c(to_del,name)
  }
}

z_scores = z_scores[,!(names(z_scores) %in% to_del)]

### Omit Attributes with majority small z scores ###
instances = nrow(z_scores)

maj = instances %/% 2

sig_thresh = 1

smallz = c()
for(name in names(z_scores)){ #for each column
  current_col = z_scores[,name]
  current_count = 0
  for(val in current_col){ #find num abs(z scores) < sig_thresh
    if(abs(val) <= sig_thresh){current_count = current_count + 1}
  }
  if(current_count > maj){smallz = c(smallz,name)} #if count > maj -> remove
}

z_scores = z_scores[,!(names(z_scores) %in% smallz)]




### Add in Sample Histology ###

patient = as.data.frame(read.csv("data_clinical_patient.csv", sep = "\t", stringsAsFactors = FALSE, header = TRUE, skip = 4))
p_nam = patient[,1]

p_nam_fix = c()
for(name in p_nam){p_nam_fix = c(p_nam_fix,substr(name,11,nchar(name)))}

patient = patient[,2:ncol(patient)]

row.names(patient) = p_nam_fix

stat = data.frame(patient[,"OS_STATUS"])
colnames(stat) = "OS_STATUS"
rownames(stat) = p_nam_fix

rownames(z_scores) %in% rownames(stat)

mg = merge.data.frame(stat,z_scores, by="row.names")

cert = mg[!(mg$OS_STATUS == ""),]

write.csv(cert,file = "z_score_status.csv")
