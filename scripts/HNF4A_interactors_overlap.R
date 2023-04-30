setwd("/media/user/Edison1/GeneRank_Ruth/Rahul_analysis/jupyter-ttest/Step-3 top131(123mappedCPDB)Genes")

write.table(gene.data[order(gene.data$p_value),][,1][1:131]
            , "/media/user/Edison1/GeneRank_Ruth/Rahul_analysis/jupyter-ttest/Step-3 top131(123mappedCPDB)Genes/Top131genes.txt"
            , quote = F
            , row.names= F)

IG_high <- read.table("/media/user/Edison1/GeneRank_Ruth/Rahul_analysis/jupyter-ttest/Step-3 top131(123mappedCPDB)Genes/Top131genes.txt"
                      ,as.is=T, header = T)$x

######### InnateDB
#innateDB_HNF4A = read.table("/media/user/Edison/GeneRank_Ruth/Rahul_analysis/Information Gain Only/HNF4A_target_Genes/InnateDB_HNF4A.txt"
#                            , sep = "\t"
#                            , as.is = T
#                            , fill = NA)

innateDB_HNF4A = read.table("/media/user/Edison1/GeneRank_Ruth/Rahul_analysis/Information Gain Only/HNF4A_target_Genes/InnateDB_HNF4A_fromWeb.txt"
                            , sep = "\t"
                            , as.is = T
                            , fill = NA)


innateDB_HNF4A_genes = innateDB_HNF4A[,1]
targets_HNF4A <- sub(innateDB_HNF4A_genes,
                     pat="  interacts with HNF4A",rep="")

targets_HNF4A <- sub(targets_HNF4A,
                     pat="HNF4A  interacts with ",rep="")

targets_HNF4A <- sub(targets_HNF4A,
                     pat="HNF4A  physically interacts with ",rep="")
targets_HNF4A <- sub(targets_HNF4A,
                     pat="  physically interacts with HNF4A",rep="")



targets_HNF4A <- sub(targets_HNF4A,
                     pat="Colocalization of ",rep="")

targets_HNF4A <- sub(targets_HNF4A,
                     pat="  and HNF4A",rep="")

targets_HNF4A <- sub(targets_HNF4A,
                     pat="Colocalization of HNF4A  and ",rep="")

targets_HNF4A <- sub(targets_HNF4A,
                     pat="Genetic interaction between HNF4A and ",rep="")

targets_HNF4A <- sub(targets_HNF4A,
                     pat="HNF4A  and ",rep="")

targets_HNF4A <- sub(targets_HNF4A,
                     pat="HNF4A physically associates with ",rep="")


targets_HNF4A <- sub(targets_HNF4A,
                     pat="Genetic interaction between ",rep="")# Genetic interaction between HNF4A  and 

## Remove character vector "fullname"
targets_HNF4A  = targets_HNF4A[-2333]

# Find the unique genes
targets_HNF4A = unique(targets_HNF4A)

# remove HNF4A from the target gene list
targets_HNF4A = targets_HNF4A[-1197]

write.table(file = "HNF4A_targets.txt", targets_HNF4A, quote = F, row.names = F)
targets_HNF4A = read.table("HNF4A_targets.txt", header = T)



HNF4A_target_genes = intersect(as.character(targets_HNF4A$x), as.character(gene.data$Gene[1:131]))


length(HNF4A_target_genes)
write.table(file = "26targetgenesofHNF4A_intheNetwork.txt"
		, HNF4A_target_genes
		, quote = F
		, row.names = F)

## Get the top 100 HNF4A target genes based on Info Gain
#gene.data = read.csv(file = "/media/user/Edison1/GeneRank_Ruth/Rahul_analysis/jupyter-ttest/Step -2 Probe2Gene/ProbeSummarized2gene_usingTtest.csv")[, -1]

Top100HNF4Atargets = gene.data[which(as.character(gene.data$Gene) %in% as.character(targets_HNF4A$x) == TRUE),][27:126, "Gene"]
Top200HNF4Atargets = gene.data[which(as.character(gene.data$Gene) %in% as.character(targets_HNF4A$x) == TRUE),][27:226, "Gene"]
setwd("/media/user/Edison1/GeneRank_Ruth/Rahul_analysis/jupyter-ttest/Step 5 - HNF4A interactors")
write.table(Top100HNF4Atargets, "Top100targetsofHNF4A.txt")
write.table(Top200HNF4Atargets, "Top200targetsofHNF4A.txt")


SeedlistBackground = read.table("/media/user/Edison1/GeneRank_Ruth/Rahul_analysis/jupyter-ttest/Step-3 top131(123mappedCPDB)Genes/backgroundgenes_13126.txt"
                            , header = T)$x
# background overlap
length(intersect(as.character(SeedlistBackground), as.character(targets_HNF4A$x)))

write.table(intersect(SeedlistBackground, as.character(targets_HNF4A$x)), "ttest-BackgroundList_HNF4A-targets.txt"
            , row.names = F, col.names = F, quote = F)

### TOtal seedlist background = 13126  #  length(gene.data[order(-gene.data$p_value),][,1])
### Seed gene list = 131   # length(IG_high)

### Total db background = 2281  # dim(targets_HNF4A)[1]
### OVerlap btwn TOtal seedlist background and Total db background = 1878

### Overlap btwn seedlistgene list and Total db background = 26 # length(intersect(levels(droplevels(targets_HNF4A$x)), IG_high))

## chi
## 10778	1554
## 500		65
## chi - 0.49647, p < 0.4811
M1 <- chisq.test(as.table(rbind(c(13126, 1878), c(131, 26))))



## Save HNF4A target genes' expression ##
save(gene_expression
     , file = "/media/user/Edison1/GeneRank_Ruth/Rahul_analysis/jupyter-ttest/Step -2 Probe2Gene/gene_expression.RData")
# load file gene_expression
load(file = "/media/user/Edison1/GeneRank_Ruth/Rahul_analysis/jupyter/Step -2 Probe2Gene/gene_expression.RData")

# Read an old file containing clinical data 
probe_expr_clinical = read.csv("/media/user/Edison1/GeneRank_Ruth/Rahul_analysis/jupyter-ttest/Step-3 top131(123mappedCPDB)Genes/ProbeExpressionClinical.csv")
probe_expr_clinical = probe_expr_clinical[, c("Patient", "RFS_event")]
names(probe_expr_clinical) = c("Patient", "RFS_event")

gene_expression$Patient = row.names(gene_expression)

de <- merge(gene_expression, probe_expr_clinical, by.x = "Patient", by.y = "Patient")
de = de[, -1] # remove ID column
write.csv(file = "/media/user/Edison1/GeneRank_Ruth/Noa Deep learning/DeepHNF4A/DeepHNF4A/forPaper_ttest/Gene_expression.csv"
          , de
          , row.names = F)

# get the HNF4A targets ordered by p-value and filter out the top 26 genes which are present in top 1%
topgenes_filtered = gene.data[which(as.character(gene.data$Gene) %in% as.character(targets_HNF4A$x) == TRUE),][-c(1:26), "Gene"]


HNF4A_target_genes = intersect(as.character(topgenes_filtered), as.character(targets_HNF4A$x))
HNF4A_targets_expression = de[, c(HNF4A_target_genes, "RFS_event.x")]     

HNF4A_targets_expressionTop100 = de[, c(as.character(Top100HNF4Atargets), "RFS_event.x")]     

write.csv(file = "/media/user/Edison1/GeneRank_Ruth/Noa Deep learning/DeepHNF4A/DeepHNF4A/forPaper_ttest/HNF4A_targets_expression.csv"
          , HNF4A_targets_expression
          , row.names = F)

write.csv(file = "/media/user/Edison1/GeneRank_Ruth/Noa Deep learning/DeepHNF4A/DeepHNF4A/forPaper_ttest/HNF4A_targets_expression100.csv"
          , HNF4A_targets_expressionTop100
          , row.names = F)


#### Calculate chi-square test ####
## From Agresti(2007) p.39
calculate_chi <- function(M){
  chisq = list()
  for(i in 1:length(M)){
    k = M[[i]]
    dimnames(k) <- list(Group = c("Group1", "Group2"),
                        GenesCount = c("Background", "Signature"))
    (Xsq <- chisq.test(k))  # Prints test summary
    print(Xsq)
    chisq[[i]] = Xsq$p.value
  }
  return(chisq)
}
M1 <- as.table(rbind(c(13101, 1878), c(131, 25)))
M2 <- as.table(rbind(c(13101, 1878), c(100, 21)))
M3 <- as.table(rbind(c(13101, 1878), c(65, 16)))
M4 <- as.table(rbind(c(13101, 1878), c(44, 11)))
M4 <- as.table(rbind(c(13101, 1878), c(33, 7)))
M5 <- as.table(rbind(c(13101, 1878), c(26, 7)))
M6 <- as.table(rbind(c(13101, 1878), c(22, 4)))

M = list(M1, M2, M3, M4, M5, M6)
mine = calculate_chi(list(M1, M2))


plot(c(1, 0.5, 0.25, 0.125, 0.06, 0.03), c(0.2, 0.1, 0.07, 0.47, 0.21, 0.88))


Xsq$observed   # observed counts (same as M)
Xsq$expected   # expected counts under the null
Xsq$residuals  # Pearson residuals
Xsq$stdres     # standardized residuals

# X-squared = 3.3561, df = 1, p-value = 0.06696


boxplot_randomizedTrial <- function(data){
  library(reshape)
  meltData <- melt(data)
  png(filename = "/media/user/Edison/GeneRank_Ruth/Noa Deep learning/DeepHNF4A/DeepHNF4A/for Paper/RandomSample.png"
      , width = 1196
      , height = 510)
  boxplot(data=meltData, value~variable
          , ylim = c(0.6,0.85) 
          , col = c("royalblue2", "cyan", "blue", "magenta")
          , main = paste("AUC and accuracy of randomly sampled ML",  "\n" , "models compared to HNF4A based models")
          , xlab = "ML metrics"
          , ylab = "Fitness of the model"
          , cex.axis = 1.3
          , font = 2
          , font.lab = 2
          , cex.lab = 1.3
  )
  points(1, 0.78, type = , col = "red", pch = 20,cex = 1.5 )
  points(2, 0.72, type = , col = "green", pch = 20, cex = 1.5)
  points(3, 0.76, type = , col = "red", pch = 20, cex = 1.5)
  points(4, 0.71, type = , col = "green", pch = 20, cex = 1.5)
  dev.off()
}

data1 <- read.table("/media/user/Edison1/GeneRank_Ruth/Noa Deep learning/DeepHNF4A/DeepHNF4A/forPaper_ttest/Targets_1852_randomSample100/statistics.txt") #[, c(1, 2)]

data <- read.table("/media/user/Edison1/GeneRank_Ruth/Noa Deep learning/DeepHNF4A/DeepHNF4A/forPaper_ttest/RandomSample1852/statistics.txt")#[, c(3, 4)]
data2 = cbind(data, data1[,1])
data3 = cbind(data2, data1[,2])

names(data3) = c(  "AUC-random"
                   , "Accuracy-random"
                   , "AUC-HNF4A targets"
                   , "Accuracy-HNF4A targets")

boxplot_randomizedTrial(data3)