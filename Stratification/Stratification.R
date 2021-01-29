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

### Read in Patient data and format ###

patient = as.data.frame(read.csv("data_clinical_patient.csv", sep = "\t", stringsAsFactors = FALSE, header = TRUE, skip = 4))
p_nam = patient[,1]

p_nam_fix = c()
for(name in p_nam){p_nam_fix = c(p_nam_fix,substr(name,11,nchar(name)))}

patient = patient[,2:ncol(patient)]

row.names(patient) = p_nam_fix

### Merge Frames together ###
patient = merge.data.frame(patient,z_scores, by="row.names")

### Create Subsets ###

five_yr = 365*5

surv_g_5 = patient[patient$OS_DAYS > five_yr,]

surv_g_5 = surv_g_5[(is.na(surv_g_5$OS_STATUS)),]

surv_u_5 = patient[!(patient$OS_DAYS) > five_yr,]

EFS_g_5 = patient[patient$EFS_TIME > five_yr,]
EFS_l_5 = patient[!(patient$EFS_TIME > five_yr),]

subs = c(surv_g_5,surv_u_5,EFS_g_5,EFS_l_5)






### Write to files ###

write.csv(surv_g_5, file="overall_survival_g5.csv")
write.csv(surv_u_5, file="overall_survival_l5.csv")

write.csv(EFS_g_5, file="overall_survival_")