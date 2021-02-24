#File used to understand foreach command for parallel processing of clusters

library('foreach')
library('parallel')
library('flexclust') # handles weighted clustering

setwd("~/Documents/GitHub/BioinfGroupProject/GSEA_rank_cluster")

#set.seed(255)
#read in and format counts (z_scores)
z_scores = read.csv("data_expression_median.txt", sep = "\t", stringsAsFactors = FALSE, header = TRUE)

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

#read in ranked file, extract names

ranking = read.csv('ranked_gene_list.tsv', sep = "\t")
ranking = ranking$NAME


rz = z_scores[,colnames(z_scores) %in% ranking] #get genes from count included in ranking

rk = ranking[ranking %in% colnames(z_scores)] # get ranked genes included in counts

rz_order = rz[,rk] # order columns by ranking



cl = makeCluster(detectCores() -1)
cl

get_WCSS = function(clustering, data){
  cents = as.data.frame(parameters(clustering))
  
  wcss = c()
  wcss_clst = c()
  for (k in 1:nrow(cents)){ # for each cluster
    clust_cent = cents[k,] 
    clust_dat = data[names(which(clusters(clustering) == k)),]    # need to make sure only select data from relevant cluster
    
    for (atrib in names(clust_cent)){ # for each attribute
      wcss_peratrib = sum(abs(as.numeric(clust_cent[atrib])-clust_dat[atrib])^2)
      wcss_clst = c(wcss_clst,wcss_peratrib)
    }
    wcss = c(wcss,sum(wcss_clst))
    wcss_clst = c()
  }
  
  return(sum(wcss))
} #returns WCSS for given cluster
convergence = function(num_seeds, k_clust, data, weights = wg){ #function to determine seed that optimises given k
  seeds = sample(1:1000, num_seeds, replace = FALSE) # Select random seeds to itterate through
  
  
  clust_frame = data.frame(seed = c(), wcss = c())
  
  for(seed in seeds){
    set.seed(seed)
    
    cls = cclust(data, method = "hardcl", k=k_clust, weights = weights)
    cls_WCSS = get_WCSS(cls,data)
    
    clust_frame = rbind(clust_frame, data.frame(seed = seed, wcss = cls_WCSS))
  }
  
  return(clust_frame)
}

#optimise sd - which sd value gives the lowest WCC?
sdev = seq(100,3000, by=100)

bestWCC = Inf
bestDev = 0 #1400

WCCs = foreach(s = sdev) %dopar% {
  wg <- dnorm(seq(1, ncol(rz_order)), mean = ncol(rz_order)/2, sd = sdev)
  range01 <- function(x){(x-min(x))/(max(x)-min(x))}
  wg = range01(wg)
  
  cv = convergence(100,2,rz_order, weights=wg) #change this to allow convergence
  WCC = min(cv$wcss)
  WCC
}

bestDev = sdev[which.min(WCCs)]

wg <- dnorm(seq(1, ncol(rz_order)), mean = ncol(rz_order)/2, sd = bestDev)

range01 <- function(x){(x-min(x))/(max(x)-min(x))}

wg = range01(wg)


elb = foreach(k = 2:10) %dopar% {
  cv = convergence(100,k,rz_order,weights=wg)
  WCC = min(cv$wcss)
  WCC
}

elb = unlist(elb)

plot(2:10, elb, type='l', main = "Elbow Plot", ylab = "WCSS", xlab = "K" )

elb_frame = data.frame(k_val = 2:10, WCSS = elb)

write.table(elb_frame, "elbow_para.txt", sep = "\t")





