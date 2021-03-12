#File used to analyse GSEA enrichment result outputs. Used to find unique signatures and overlaps in 

setwd("~/Documents/GitHub/BioinfGroupProject/rank_cluster_2/k5/analysis")

get_sig = function(set, cluster, pval=0.05, qval=0.25){ # get significant GSEA values from a specific cluster and gene set
  file = paste(set, "/gsea_report_for_C",cluster,".tsv", sep = "")
  res = read.csv(file, sep = "\t")
  return(res[(res[,"NOM.p.val"] <= pval) & (res[,"FDR.q.val"] <= qval),])
}

get_sig('Hallmark', 2)


find_overlap = function(c1, c2){ #find overlap between enrichment sets in two clusters
  ol = intersect(c1[,"NAME"], c2[,"NAME"])
  c1_ol = c1[c1[,"NAME"] %in% ol,]
  c2_ol = c2[c2[,"NAME"] %in% ol,]
  ol_frame = data.frame(c1_p = c1_ol[,"NOM.p.val"],
                        c2_p = c2_ol[,"NOM.p.val"],
                        c1_q = c2_ol[,"FDR.q.val"],
                        c2_q = c2_ol[,"FDR.q.val"])
  row.names(ol_frame) = ol
  return(ol_frame)
}

find_overlap_multi = function(clsts, set, pval= 0.05, qval= 0.25){ # Finds intersection of multiple clusters
  res = list()
  for(i in 1:length(clsts)){ # may cause error
    #print(clsts[i])
    res[i] = ifelse(length(get_sig(set, clsts[i], pval, qval)[,"NAME"]) > 0, get_sig(set, clsts[i], pval, qval), NA) #dim(find_unique(i,clsts[-i], set))[1]==0
  }
  
  return(Reduce(intersect, res))
  
  
  sec = c()
  
  
  
  
  for (i in 1:length(clsts)){
    if (i +1 <= length(clsts)){
      sec = c(sec,intersect(res[[i]], res[[i+1]]))
    }
  }
  return(sec)
}


#find_overlap(get_sig('Hallmark', 2), get_sig('Hallmark', 8))

find_overlap_multi(clsts,"Hallmark")




#find_overlap_multi(c(2,1), "Hallmark")
#find_overlap_multi(c(2,8,9), "Hallmark")

#find_overlap_multi(c(6,8), "Hallmark")
#find_overlap_multi(c(8,9), "Hallmark")

#find_overlap_multi(c(5,6), "Hallmark")
#find_overlap_multi(c(1,2,3), "Hallmark")
#find_overlap_multi(c(8,9), "Hallmark")


#find_overlap_multi(c(2,8), "COVID_PPI")

#find_overlap_multi(c(1,5), "Curated")

find_separate = function(c1, c2){ #find unique enrichment sets in two clusters
  ol = intersect(c1[,"NAME"], c2[,"NAME"])
  c1_sep = c1[!(c1[,"NAME"] %in% ol),]
  c2_sep = c2[!(c2[,"NAME"] %in% ol),]
  
  return(list(c1_sep, c2_sep))
}


find_unique = function(c, remain, set){ #find all unique gene sets in cluster
  r_tot = c() #all gene sets in remaining clusters
  for (r in remain){
    rem = get_sig(set,r)[,"NAME"]
    r_tot = c(r_tot, rem)
  }
  
  r_tot = unique(r_tot)
  c_sig = get_sig(set,c) # all gene sets in given cluster
  #print(c_sig[,"NAME"])
  
  c_uni = c_sig[!(c_sig[,"NAME"] %in% r_tot),] #return gene sets present only in cluster
  return(c_uni)
}

find_all_unique = function(set, clsts){ #finds all unique gene sets for each cluster
  ret = list()
  for (i in clsts){
    ret[i] = ifelse(dim(find_unique(i,clsts[-i], set))[1]==0, NA, find_unique(i, clsts[-i], set))
  }
  return(ret)
}

#Idea for potential function:
#for cluster, open each unique set, extract leading edge (store in list). Find overlap between each leading edge[like correlation matrix]
# use matrix to build network. Find most central nodes

read_ol = function(all_over, clst_a, clst_b){ 
  select = all_over[all_over[,1] == clst_a & all_over[,2] == clst_b,]
  return(strsplit(select[,3], split = ", "))
}

all_overlap = function(set, clsts, pval=0.05, qval=0.25){ #Returns df containing all clusters with overlaps between each other
  ol = data.frame() #initialize df
  
  for(i in clsts){ #get all pairs of GSEA results
    sig_i = get_sig(set, i, pval, qval)
    for(j in clsts){
      if(i != j){ #ensure don't find match between same GSEA result
        sig_j = get_sig(set,j, pval, qval)
        fo = find_overlap(sig_i, sig_j) #find the overlap
        if (dim(fo)[1]!=0){ # if there is a detectable overlap
          set_nms = paste(rownames(fo), collapse='","') #convert all sets to a single string
          ol = rbind(ol,c(i,j,nrow(fo),set_nms)) #create a df with cluster names, size of overlap and items included in overlap
        }
      }
    }
  }
  if (ncol(ol) >= 2){ # if an overlap was found
    ol = ol[!duplicated(apply(ol[,1:2],1,function(x) paste(sort(x),collapse=''))),] #removes duplicate reversed entries
    rwnm = paste(ol[,1],"&", ol[,2]) # create rowname out of first and second entries
    ol = ol[,3:ncol(ol)]
    row.names(ol) = rwnm
    colnames(ol) = c("Overlap Size", "Overlap Sets")
    
    
  }
  return(ol)
}

ao_hall =all_overlap('Hallmark', c(1,3))
#find_overlap(get_sig('Hallmark', 6),get_sig('Hallmark', 8))


#find_overlap_multi(c(9,4), 'Hallmark')

#find_overlap_multi(c(2,8), 'GO_BP', pval= 0.005, qval = 0.1)

#find_overlap_multi(c(2,8,6,9), 'GO_BP', pval = 0.054)



clsts = c(1:5)




#sets = c("Hallmark", "Curated", "COVID", "COVID_PPI", "GO_BP")

#for (set in sets){ #Iterate through each set, find all unique values for each cluster, write to a file in "Unique" folder
#  fnm = paste("Unique/Unique_",set,".txt", sep = "")
#  sink(fnm)
#  uni = find_all_unique(set, clsts)
#  print(uni)
#  sink()
  
  # Create a summary file in same directory that lists the number of unique gene sets
#  len_u = c()
#  for(u in uni){len_u = c(len_u, length(u))}
  
#  write.table(data.frame("Cluster" =  clsts, "Length" = len_u), paste("Unique/UniqueSummary_", set, ".tsv", sep = ""), sep = "\t", row.names = F)
#}

#for (set in sets){ #Iterate through each set, find all overlap for each cluster, write to a file
#  ao = all_overlap(set, clsts)
#  fnm = paste("Overlap/",set, "_overlap.csv", sep ="")
#  write.csv(ao,fnm)
#}
