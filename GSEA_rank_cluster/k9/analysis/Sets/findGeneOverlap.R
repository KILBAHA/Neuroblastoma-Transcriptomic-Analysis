setwd("~/Documents/GitHub/BioinfGroupProject/GSEA_rank_cluster/k9/analysis/Sets")


files = read.table('files.txt')

files = files$V1



get_ol = function(compare, sets){
  cmp = read.csv(paste(compare, ".tsv", sep = ""), sep = "\t")
  
  lst = list()
  for(set in sets){
    current = read.csv(paste(set,".tsv", sep = ""), sep="\t")
    lst[set] = ifelse(length(current$SYMBOL[current$SYMBOL %in% cmp$SYMBOL]) > 0 , list(current$SYMBOL[current$SYMBOL %in% cmp$SYMBOL]), NA)
  }
  return(lst)
}

get_ol("E2", c("HALLMARK_APOPTOSIS", "HALLMARK_COMPLEMENT", "HALLMARK_IL2_STAT5_SIGNALING", "HALLMARK_IL6_JAK_STAT3_SIGNALING", "HALLMARK_INFLAMATORY_RESPONSE", "HALLMARK_INTERFERON_ALPHA_RESPONSE", "HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_P53_PATHWAY"))

uni2 = c("HALLMARK_PI3K_AKT_MTOR_SIGNALING", "HALLMARK_ALLOGRAFT_REJECTION", "HALLMARK_MITOTIC_SPINDLE", "HALLMARK_UV_RESPONSE_UP")
uni3 = c("HALLMARK_MYC_TARGETS_V1", "HALLMARK_MYC_TARGETS_V2", "HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT")
uni6 = c("HALLMARK_PROTEIN_SECRETION", "HALLMARK_PEROXISOME", "HALLMARK_KRAS_SIGNALING_DN")


get_ol("SARS CORONAVIRUS NUCLEOCAPSID PROTEIN FROM VIRUS-HOST PPI P-HIPSTER 2020", uni2)
get_ol("SARS CORONAVIRUS NSP4-PP1A FROM VIRUS-HOST PPI P-HIPSTER 2020", uni2)
#get_ol("SARS-COV-2/HUMAN INTERACTOME GENE SET FROM GUZZI", uni2)
get_ol("COVID19-NSP9 PROTEIN HOST PPI FROM KROGAN", uni2)

get_ol("COVID19-NSP8 PROTEIN HOST PPI FROM KROGAN", uni3)

get_ol("SARS CORONAVIRUS LEADER PROTEIN FROM VIRUS-HOST PPI P-HIPSTER 2020", uni6)
get_ol("SARS CORONAVIRUS FORMERLY KNOWN AS GROWTH-FACTOR-LIKE PROTEIN FROM VIRUS-HOST PPI P-HIPSTER 2020", uni6)



E2 = read.csv("E2.tsv", sep = "\t")

lst = list()
for(file in files){
  current = read.csv(file, sep = "\t")
  lst[file] = ifelse(length(current$SYMBOL[current$SYMBOL %in% E2$SYMBOL]) > 0 , list(current$SYMBOL[current$SYMBOL %in% E2$SYMBOL]), NA)
}
lst





