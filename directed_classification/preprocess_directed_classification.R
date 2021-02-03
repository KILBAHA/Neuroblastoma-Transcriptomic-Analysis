
### Read in z_scores and format ###
z_scores = read.csv("data_RNA_Seq_mRNA_median_all_sample_Zscores.txt", sep = "\t", stringsAsFactors = FALSE, header = TRUE)

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

## Read in Clinical Data ###
patient = as.data.frame(read.csv("data_clinical_patient.txt", sep = "\t", stringsAsFactors = FALSE, header = TRUE, skip = 4))
p_nam = patient[,1]

p_nam_fix = c()
for(name in p_nam){p_nam_fix = c(p_nam_fix,substr(name,11,nchar(name)))}

patient = patient[,2:ncol(patient)]

row.names(patient) = p_nam_fix

status = data.frame(patient[,"OS_STATUS"])
colnames(status) = "STATUS"
rownames(status) = p_nam_fix

mg = merge.data.frame(status,z_scores, by="row.names")


check_attrib = function(df, attrib){ #Script to determine if an attribute is in the dataframe
  ret = attrib %in% colnames(df)
  return(ret)
}


attributes = c("NTRK1", "MYCN", "MDM2", "ALK", "CHD5", "CADM1",
               "CD44", "CD-133", "KIT", "NTRK2", "DLK1","STATUS")


attributes[!check_attrib(mg, attributes)]

#write.csv(colnames(mg), "gene_list.csv")

select_attributes = attributes[check_attrib(mg,attributes)]

dr = mg[mg$STATUS!="",select_attributes]
#unique(dr$STATUS)

write.csv(dr, "dir_class.csv")


