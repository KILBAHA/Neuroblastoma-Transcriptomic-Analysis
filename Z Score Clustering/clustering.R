### Import modules ###

### Set wd ###
setwd("~/Documents/Scripts")


### Read in z_scores and format ###
z_scores = read.csv("data_mRNA_median_Zscores.txt", sep = "\t", stringsAsFactors = FALSE, header = TRUE)

map = z_scores[,c(1,2)]

inst_names = colnames(z_scores)[3:length(colnames(z_scores))]

z_scores = z_scores[,-2]
z_scores = as.data.frame(t(z_scores), stringsAsFactors = FALSE)
names(z_scores) = map[,1]
z_scores = z_scores[-1,]

z_scores = as.data.frame(apply(z_scores, 2, as.numeric))

row.names(z_scores) = inst_names

### Add Prognosis Markers ###

markers = as.data.frame(read.csv('marker.csv'))


up_thresh = 1
good_prog = markers[markers[,"Prognosis"] ==1,]
good_progU = good_prog[good_prog$Expression ==1,]
for (mk in good_progU[,"Marker"]){
  ifelse(z_scores[,mk] >= up_thresh, print("true"), print("false"))# fix this!
}




good_progD = good_prog[good_prog$Expression == 0,]


bad_prog = markers[markers$Prognosis == 0,]
bad_progU = bad_prog[bad_prog$Expression ==1,]
bad_progD = bad_prog[bad_prog$Expression == 0,]


label_prog = function(data, prog_markers){
  #find good prognosis markers
  
  good_prog = 
  
  
  
  #for each instance in df:
    #note prognosis
  #find bad prognosis markers
  #for each instance in df:
    #note prognosis
  #bind prognosis to dataframe
}





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


### Perform PCA ###

z_scores_stand = as.data.frame(scale(z_scores))

#sapply(z_scores_stand,sd) #now, standard deviations are 1
#sapply(z_scores_stand,mean) #now, mean should be 0 (or very very close to 0)

pca = prcomp(z_scores_stand,scale=T) 

summary(pca)
pca = as.data.frame(pca$x[,1:180]) #contains 95% of variance of dataset


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

### dim reduced PCA ###

z_scores_stand = as.data.frame(scale(z_scores))

#sapply(z_scores_stand,sd) #now, standard deviations are 1
#sapply(z_scores_stand,mean) #now, mean should be 0 (or very very close to 0)

pca2 = prcomp(z_scores_stand,scale=T) 

summary(pca2)
pca2 = as.data.frame(pca2$x[,1:124]) #contains 95% of variance of dataset



### Clustering ###

library('cluster')

find_hierac_clust = function(data,linkage){
  dist_mat = dist(data, method = 'euclidian')
  hclust = hclust(dist_mat, method = linkage)
  plot(hclust)
  cut_avg = cutree(hclust, k = 2)
  return(cut_avg)
}


z_clustAv = find_hierac_clust(z_scores,'average')
z_clustComp = find_hierac_clust(z_scores, 'complete')
z_clustSing = find_hierac_clust(z_scores, 'single')
z_clustCent = find_hierac_clust(z_scores, 'centroid')

hist(z_clustAv)
hist(z_clustComp) #--> good deliniation
hist(z_clustSing)
hist(z_clustCent)


pc_clustAv = find_hierac_clust(pca,'average')
pc_clustComp = find_hierac_clust(pca, 'complete')
pc_clustSing = find_hierac_clust(pca, 'single')
pc_clustCent = find_hierac_clust(pca, 'centroid')


hist(pc_clustAv)
hist(pc_clustComp) #--> better deliniation!
hist(pc_clustSing)
hist(pc_clustCent)


pc2_clustAv = find_hierac_clust(pca2,'average')
pc2_clustComp = find_hierac_clust(pca2, 'complete')
pc2_clustSing = find_hierac_clust(pca2, 'single')
pc2_clustCent = find_hierac_clust(pca2, 'centroid')


hist(pc2_clustAv)
hist(pc2_clustComp) #--> best
hist(pc2_clustSing)
hist(pc2_clustCent)




### Attach Clustering lbls and export as .csv ###

z_scores$cluster = pc2_clustComp

write.csv(z_scores, file="dr_z.csv")
