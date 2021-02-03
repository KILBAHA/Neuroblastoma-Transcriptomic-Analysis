setwd("~/Documents/GitHub/BioinfGroupProject/GSEA")


### Read in and Format count data ###
count = read.csv("data_expression_median.txt", sep = "\t")


count = count[,3:ncol(count)]


inst_names = colnames(count)

name_fix = c()
for(name in inst_names){name_fix = c(name_fix,substr(name,11,nchar(name)-3))}

colnames(count) = name_fix

count = t(count)[,1:2] #don't need the data, just need to transpose the df 


### Merge with Clinical data ###

patient = as.data.frame(read.csv("data_clinical_patient.csv", sep = "\t", stringsAsFactors = FALSE, header = TRUE, skip = 4))
p_nam = patient[,1]

p_nam_fix = c()
for(name in p_nam){p_nam_fix = c(p_nam_fix,substr(name,11,nchar(name)))}

patient= as.data.frame(patient[,"OS_STATUS"])

row.names(patient) = p_nam_fix
colnames(patient) = c("SURVIVAL")

row.names(patient)%in% row.names(count)
row.names(count)%in%row.names(patient)

mg = merge.data.frame(count,patient, by="row.names")
mg = mg[,c(1,4)]
nodat = mg[mg$SURVIVAL == "",][,1]
nodat #need to remove these columns from Expr file

mg = mg[!(mg$SURVIVAL ==""),]

#number of samples
num_samp = nrow(mg)
#number of groups
num_groups = length(unique(mg$SURVIVAL))


#first group name
nam_fg = unique(mg$SURVIVAL)[1] #can change order of indicies 
#second group name
nam_sg = unique(mg$SURVIVAL)[2]

cat(num_samp, " ",num_groups, " ", "1", "\n", "#"," ",nam_fg, " ", nam_sg, "\n", file="PheLab.cls")

cat(mg$SURVIVAL, file = "PheLab.cls", sep = "\t", append = T)



