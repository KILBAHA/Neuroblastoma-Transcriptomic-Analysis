setwd("~/Documents/GitHub/BioinfGroupProject/GSEA_rank_cluster/evaluate_clusters")

#Load Clusters
c1 = read.csv("c1.txt", sep = "\t", header= F)
c2 = read.csv("c2.txt", sep = "\t", header= F)
Eg5 = read.csv("EF_g5.txt", sep = "\t", header=F)
Eu5 = read.csv("EF_u5.txt", sep = "\t", header = F)

c("a", "b") %in% c("1", "a", "2", "b", "c")
c("1", "a", "2", "b", "c") %in% c("a", "b")

jacard = function(a,b){ #returns percentage similarity of two sets
  return(length(intersect(a,b))/length(union(a,b))*100)
}

overlap = function(a,b){
  return(length(intersect(a,b))/min(length(a),length(b))*100)
}


jframe = data.frame("HIGH_EFS" = c(jacard(c1$V1,Eg5$V1),jacard(c2$V1,Eg5$V1)),
                    "LOW_EFS" = c(jacard(c1$V1,Eu5$V1),jacard(c2$V1,Eu5$V1)))

rownames(jframe) = c("c1","c2")


oframe = data.frame("HIGH_EFS" = c(overlap(c1$V1,Eg5$V1),overlap(c2$V1,Eg5$V1)),
                    "LOW_EFS" = c(overlap(c1$V1,Eu5$V1),overlap(c2$V1,Eu5$V1)))
rownames(oframe) = c("c1","c2")


jacard(c1$V1,Eg5$V1)
jacard(c1$V1,Eu5$V1)

jacard(c2$V1,Eg5$V1)
jacard(c2$V1,Eu5$V1) #cluster 2 more than 50% identical to EFS < 5



overlap(c1$V1,Eg5$V1)
overlap(c1$V1,Eu5$V1)

overlap(c2$V1,Eg5$V1)
overlap(c2$V1,Eu5$V1) #cluster 2 more than 50% identical to EFS < 5








