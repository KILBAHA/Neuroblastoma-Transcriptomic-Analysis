setwd("~/Documents/GitHub/BioinfGroupProject/dim_red_definitive")


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

srv = mg[mg$STATUS == "0:LIVING",]
dec = mg[mg$STATUS == "1:DECEASED",]

get_mean = function(dat){
  atribs = colnames(dat)[colnames(dat)!= "Row.names" & colnames(dat)!= "STATUS"]
  ret = c()
  for(atrib in atribs){
    mea = mean(dat[,atrib])
    ret = cbind(ret, mea)
  }
  
  ret = as.data.frame(ret)
  colnames(ret) = atribs
  return(ret)
}

srv_mean = get_mean(srv)
dec_mean = get_mean(dec)

get_deg = function(d1, d2, thresh=0.5){ #change thresh to correct for multiple significance testing
  ab = abs(d1-d2)
  sig = apply(ab,2, function(imp){imp > thresh})
  deg = ab[,sig]
  return(deg)
}

deg = get_deg(srv_mean, dec_mean)
deg




#write.csv(cert,file = "z_score_hist.csv")