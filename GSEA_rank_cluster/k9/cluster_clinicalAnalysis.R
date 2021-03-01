clinDat = read.csv("data_clinical_patient.txt", sep = "\t", skip = 4)
clinNames= clinDat[,1] 
clinDat = clinDat[,2:ncol(clinDat)]

name_fix = c()
for(name in clinNames){name_fix = c(name_fix,substr(name,11,nchar(name)))}

rownames(clinDat) = name_fix

num_pie = function(dat, tit){
  lbs = unique(dat)
  vls = c()
  for (lb in lbs){
    vls = c(vls, sum(dat== lb))
  }
  pie(vls, labels = lbs, main = tit)
}

plot_cls = function(clst_num, cdat = clinDat){
  clst = as.vector(read.table(paste("c",clst_num,".txt", sep = "")))
  clstDat = clinDat[rownames(clinDat) %in% clst$V1,]
  
  par(mfrow=c(3,3))
  
  
  
  hist(clstDat$AGE, main = paste("Age cluster", clst_num), breaks = max(clstDat$AGE)-1)
  abline(v=median(clstDat$AGE), col = "blue", lwd = 2)
  num_pie(clstDat$SEX, tit = paste("Sex cluster", clst_num))
  hist(clstDat$EFS_TIME, xlab = "EFS Time (days)", main = paste("EFS Time cluster", clst_num))
  abline(v=median(clstDat$EFS_TIME, na.rm = TRUE), col = "blue", lwd = 2)
  hist(clstDat$AGE_IN_DAYS, xlab = "Age (days)", main = paste("Age (days) cluster", clst_num))
  abline(v = median(clstDat$AGE_IN_DAYS), col = "blue", lwd = 2)  
  num_pie(clstDat$RACE, tit = paste("Race cluster", clst_num))
  num_pie(clstDat$RISK_GROUP, tit = paste("Risk Group", clst_num))
  #mtext(paste("Cluster ", clst_num, ":", sep =""), side = 3, line = -21, outer = TRUE)
  
  #par(mfrow = c(1,3))
  num_pie(clstDat$INSS_STAGE, tit = paste("INSS STAGE Cluster", clst_num))
  num_pie(clstDat$TUMOR_SAMPLE_HISTOLOGY, tit = paste("Tumour Sample Histology Cluster", clst_num))
  num_pie(clstDat$DIAGNOSIS, tit = paste("Diagnosis Cluster", clst_num))
  
} 

for(i in 1:9){
  plot_cls(i)
}





