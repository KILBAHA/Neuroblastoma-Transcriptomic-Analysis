setwd("~/Documents/GitHub/BioinfGroupProject/rank_cluster_2/k5/analysis/interaction")

get_ol = function(compare, sets){
  cmp = read.csv(paste(compare, ".tsv", sep = ""), sep = "\t")
  
  lst = list()
  for(set in sets){
    current = read.csv(paste(set,".tsv", sep = ""), sep="\t")
    lst[set] = ifelse(length(current$SYMBOL[current$SYMBOL %in% cmp$SYMBOL]) > 0 , list(current$SYMBOL[current$SYMBOL %in% cmp$SYMBOL]), NA)
  }
  return(lst)
}

uni2 = c("HALLMARK_HEDGEHOG_SIGNALING")

get_ol("COVID19-M PROTEIN HOST PPI FROM KROGAN", uni2)



E2 = read.csv("E2.tsv", sep = "\t")

lst = list()
for(file in files){
  current = read.csv(file, sep = "\t")
  lst[file] = ifelse(length(current$SYMBOL[current$SYMBOL %in% E2$SYMBOL]) > 0 , list(current$SYMBOL[current$SYMBOL %in% E2$SYMBOL]), NA)
}
lst





