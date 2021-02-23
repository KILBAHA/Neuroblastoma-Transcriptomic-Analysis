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



#optimise sd - which sd value gives the lowest WCC?
#sdev = seq(100,3000, by=100)


#bestWCC = Inf
#bestDev = 0 #1400
#for (s in sdev){
#  wg <- dnorm(seq(1, ncol(rz_order)), mean = ncol(rz_order)/2, sd = 1000)
#  range01 <- function(x){(x-min(x))/(max(x)-min(x))}
#  wg = range01(wg)
#  
#  cv = convergence(10,2,rz_order)
#  WCC = min(cv$wcss)
#  
#  if(WCC < bestWCC){
#    bestWCC = WCC
#    bestDev = s
#  }
#}

wg <- dnorm(seq(1, ncol(rz_order)), mean = ncol(rz_order)/2, sd = 1400)

range01 <- function(x){(x-min(x))/(max(x)-min(x))}

wg = range01(wg)
#wg = 1/(nrow(rz_order):1)










library('flexclust') # handles weighted clustering
library('beepr') # imports alert noise when function finished
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

#cls = cclust(rz_order, method = "hardcl", k=2, weights = wg)

#get_WCSS(cls,rz_order)




elbow = function(data,max_k,num_conv,weights = wg){ # Run elbow method to determine optimal k value
  k_rng = 2:max_k
  
  wcss_frame = data.frame(k_core = c(), wcss =c(), lw_seed = c())
  for (k in k_rng){ #for each value of K, find optimal seed for clustering. run optimal clustering and output WCSS
    cf = convergence(num_conv, k, data)
    lowest_seed = cf[which.min(cv$wcss),"seed"]
    set.seed(lowest_seed)
    cls = cclust(data, method = "hardcl", k=k, weights = weights)
    cls_WCSS = get_WCSS(cls,data)
    
    wcss_frame = rbind(wcss_frame, data.frame(k_core = k, wcss = cls_WCSS, lw_seed = lowest_seed))
  }
  
  plot(wcss_frame$k_core, wcss_frame$wcss, type='l', main = "Elbow Plot", ylab = "WCSS", xlab = "K" ) #plot WCSS for all k vals
  return(wcss_frame)
  
}
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

simpleClst = function(seed, k_clust,data, weights =wg){ #run clustering with specified seed and k
  set.seed(seed)
  cls = cclust(data, method = "hardcl", k=k_clust, weights = weights)
  return(cls)
}


#eb = elbow(rz_order, 10, 50)
#write.table(eb, file = "elbow.txt", sep = "\t")
#beep()

cv = convergence(50,2,rz_order)
write.table(cv, file="convergence.txt", sep ="\t")
hist(cv$wcss)

best_seed = cv[which.min(cv$wcss),"seed"]

best_clst = simpleClst(best_seed, 2, rz_order) #seed = 4


#cls = simpleClst(29, 2, rz_order) #best performing k=2



c1 = names(which(clusters(best_clst) == 1))
c2 = names(which(clusters(best_clst) == 2))

write.table(c1, file="c1.txt", sep = "\n", row.names = F, col.names = F, quote = T)
write.table(c2, file = "c2.txt", sep = "\n", row.names = F, col.names = F, quote = T)



              