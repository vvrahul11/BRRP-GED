setwd("/media/user/Edison1/GeneRank_Ruth/Rahul_analysis/jupyter-ttest")

#
#' @ Author: Develped by Rahul Valiya Veettil and Eitan Rubin
#' @ email: rahul.epfl@gmail.com, erubin@bgu.ac.il
#' @ University: Ben-Gurion University of the negev, Israel
#' 
#' The script takes 1809 breast cancer patient samples obtained from 
#' the scientific paper  “Gyo¨rffy et. al, An online survival analysis
#' tool to rapidly assess the effect of 22,277 genes on breast cancer
#' prognosis using microarray data of 1,809 patients, Breast Cancer 
#' Research and Treatment, 2010.  (Link : https://www.ncbi.nlm.nih.gov/pubmed/20020197)
#' 
#' ### The exact location of the data is at 
#' #   http://kmplot.com/analysis/studies/@MAS5_1000_1809_rounded_final.zip
#' 
#' 
#' 
#' Multiple correction: Idea 
#' https://stats.stackexchange.com/questions/128894/p-value-correction-for-multiple-t-tests
#'  lets take the generic example of doing 100 comparisons using a significance threshold of 0.05.
#' Now, a p-value of 0.05 means there is a 5% chance of getting that result when the null hypothesis
#' is true. Therefore, if you do these 100 comparisons, you would expect to find 5 genes significant
#' just by random chance.
#' The choice in correction can vary too. Bonferroni is a common correction but if you have 1000s
#' of genes, it is going to be exceedingly unlikely you will find anything significant because it
#' will be so conservative. In that case, you may use the FDR (False Discovery Rate) correction. 
#' 
#' set.seed(8)
#' df <- data.frame(expression=runif(1000), 
#'                  gene=rep(paste("gene", seq(250)), 4), 
#'                  treatment = rep(c("A","A","B","B"), each=250))
#' 
#' 
#' out <- do.call("rbind", 
#'     lapply(split(df, df$gene), function(x) t.test(expression~treatment, x)$p.value))
#' 
#' length(which(out < 0.05))  This was a random data and still 9 significant genes...Therefore you should
#' apply FDR
#

# Load libraries 
require(ggplot2)     # for data visualization
library(ggpubr)      # for data visualization
library(reshape2)    # Transform data between wide and long format
library(survival)    # Survival analysis
library(survplot)    # Survival analysis visualization
library(maigesPack)  # microarray
library(plyr)        # Data preprocessing
library(mixtools)    # for mixture models
library(cowplot)     # Add-on for ggplot
library(limma)       # Package for differential gene expression 
require(Biobase)     # convert data to expressionset object

# Set working directory
setwd("/media/user/Edison1/GeneRank_Ruth/Rahul_analysis/jupyter-ttest")

removeNAvalues <- function(genes){
  # The probes obtained from the Gyoffry data set was mapped to their corresponding genes 
  # using gprofiler tool. There are several probes for which gprofiler couldnot find a match
  # which were annotated with NA values.
  # 
  # This funtion is used to remove those NA values
  
  # List of all N/A probes (no result in gprofiler) and removing them from the data:
  NA_list_size = dim(genes[which(genes[,2] == 'N/A'), ])[1] #  17037 -> 2023
  gene_list_size = dim(genes[which(genes[,2] != 'N/A'), ])[1] #17037 -> 15014
  
  list_na=genes[1:NA_list_size,]
  genes_fixed=genes[1:gene_list_size,]
  a=1
  b=1
  for(i in 1:dim(genes)[1]){
    if(genes[i,2] =='N/A'){
      list_na[a,]=genes[i,]
      a=a+1
    }else{
      genes_fixed[b,]=genes[i,]
      b=b+1
    }
  }
  
  ### After removing NA values we have 15014 probes out of 17037 probes 
  genes=genes_fixed
  rm (a,b,i,genes_fixed)
  return(genes)
}

probeSummation <- function(genes, t.test.result){
  # 
  # For probe selection we used an approach where we selected the probe that 
  # is the best separator between recurred and non-recurred patients, from all possible 
  # probes that maps to a gene
  # 
  
  gene.level.data=NULL;
  genes.to.process= unique(genes$Gene);
  for (this.gene in genes.to.process) {
    #this.gene = "ATG3"
    this.probes= genes$Probe[genes$Gene==this.gene];
    #this.probes
    this.data.slice= t.test.result[t.test.result[,"Probe"] %in% this.probes,];
    #this.data.slice
    this.data.slice.means =   this.data.slice[which(this.data.slice$p_value == min(this.data.slice$p_value)), ]
    #this.data.slice.means  
    gene.level.data=rbind(gene.level.data,this.data.slice.means)
  }
  
  dim(gene.level.data)
  gene.level.data[1:5,]
  ### Make sure that the number of probes ie the 
  ### length of genes.to.process and gene.level.data$Probe are the same
  length(gene.level.data$Probe)
  length(genes.to.process)
  # We were able to map 19878 unique probes to 13101 unique genes
  
  # Combine the unique gene names to their corresponsing probe name
  gene.data=cbind(genes.to.process, gene.level.data)
  rownames(gene.data) <- 1:dim(gene.data)[1]
  # Give column names
  colnames(gene.data) <- c('Gene', 'Probe','t_test', 'p_value')
  #rm (gene.level.data, genes, ig_r_5, this.data.slice, this.data.slice.means, this.gene, this.probes, genes.to.process);
  write.csv(gene.data, file = "/media/user/Edison1/GeneRank_Ruth/Data-Noa-GeneBased/ProbeSummarized2gene_usingTtest.csv")
  #write.csv(list_na, file = "NA probes_of the expr data.csv")
  
  # ## I was able to retrieve 14114 genes and probes from the gprofiler
  return(gene.data)
}

readExpressiondataBreastCancer <- function(probe_expr){
  #
  # Read the breast cancer data
  # Seperate the clinical features from expression features
  # Transpose the data so that each row is a patient and each column is a gene
  # 
  probe.df.rownames = probe_expr$Patient
  probe_expr_clinical = probe_expr[, c(1:6, 22222:22228)]
  write.csv(file = "/media/user/Edison1/GeneRank_Ruth/Rahul_analysis/jupyter-ttest/Step-3 top131(123mappedCPDB)Genes/ProbeExpressionClinical.csv"
            , probe_expr_clinical)
  probe_expr = probe_expr[, -c(1:6, 22222, 22224:22228)]
  
  # Save the file without clinincal features
  # There was some problem with converting the character matrix to 
  # numeric matric due to memory problem. To avoid that this is one 
  # of the easiest way
  write.csv(probe_expr, "/media/user/Edison1/GeneRank_Ruth/Rahul_analysis/jupyter-ttest/Step -2 Probe2Gene/probe_expr_numeric.csv")
  probe_expr = read.csv("/media/user/Edison1/GeneRank_Ruth/Rahul_analysis/jupyter-ttest/Step -2 Probe2Gene/probe_expr_numeric.csv", header = T)
  
  RFS = probe_expr$RFS
  
  probe_expr = as.matrix(probe_expr)
  probe_expr = probe_expr[, -1]
  probe_expr_transpose = t(probe_expr)
  
  colnames(probe_expr_transpose) = probe.df.rownames
  
  probe_expr_transpose = as.data.frame(probe_expr_transpose)
  
  probe_expr_transpose$Probe = toupper(row.names(probe_expr_transpose))
  row.names(probe_expr_transpose) = NULL
  
  
  for(i in 1:dim(probe_expr_transpose)[1]){
    new_probe = paste0(strsplit(probe_expr_transpose$Probe[i], '')[[1]][-1], collapse = '')
    probe_expr_transpose$Probe[i] = new_probe
    
  }
  
  return(probe_expr_transpose)
}

HNF4A_interactors_identify <- function(){
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
  return(targets_HNF4A)
}

find_HNF4A_interactors <- function(gene.data){
  setwd("/media/user/Edison1/GeneRank_Ruth/Rahul_analysis/jupyter-ttest/Step-3 top131(123mappedCPDB)Genes")
  targets_HNF4A = HNF4A_interactors_identify()
  Top_genes = gene.data[order(-gene.data$p_value),][,1][1:500]
  HNF4A_target_genes = intersect(levels(droplevels(targets_HNF4A$x)), IG_high)
  length(HNF4A_target_genes)
  write.table(file = "65targetgenesofHNF4A_intheNetwork.txt"
              , HNF4A_target_genes
              , quote = F
              , row.names = F)
  
  ## Get the top 100 HNF4A target genes based on t-test
  gene.data = gene.data
  Top100HNF4Atargets = gene.data[gene.data$Gene %in% HNF4A_target_genes,][1:100, "Gene"]
  
  SeedlistBackground = read.table("backgroundgenes_13101.txt"
                                  , header = T)$x
  length(intersect(SeedlistBackground, as.character(targets_HNF4A$x)))
  
  write.table(intersect(SeedlistBackground, as.character(targets_HNF4A$x)), "ttest-BackgroundList_HNF4A-targets.txt"
              , row.names = F, col.names = F, quote = F)
  
  M1 <- chisq.test(as.table(rbind(c(10778, 1554), c(500, 65))))
  
  
  
  ## Save HNF4A target genes' expression ##
  save(gene_expression
       , file = "/media/user/Edison1/GeneRank_Ruth/Rahul_analysis/jupyter-ttest/Step -2 Probe2Gene/gene_expression.RData")
  # load file gene_expression
  load(file = "/media/user/Edison1/GeneRank_Ruth/Rahul_analysis/jupyter-ttest/Step -2 Probe2Gene/gene_expression.RData")
  
  # Read an old file containing clinical data 
  probe_expr_clinical = read.csv("/media/user/Edison1/GeneRank_Ruth/Rahul_analysis/jupyter-ttest/Step-3 top131(123mappedCPDB)Genes/ProbeExpressionClinical.csv")
  probe_expr_clinical = probe_expr_clinical[, c("Patient", "RFS_event")]
  names(probe_expr_clinical) = c("ID", "RFS_event")
  
  gene_expression$Patient = row.names(gene_expression)
  names(probe_expr_clinical)[1] = "Patient"
  
  de <- merge(gene_expression, probe_expr_clinical, by.x = "Patient", by.y = "Patient")
  de = de[, -1] # remove ID column
  write.csv(file = "/media/user/Edison1/GeneRank_Ruth/Noa Deep learning/DeepHNF4A/DeepHNF4A/forPaper_ttest/Gene_expression.csv"
            , de
            , row.names = F)
  
  HNF4A_target_genes = intersect(SeedlistBackground, as.character(targets_HNF4A$x))
  HNF4A_targets_expression = de[, c(HNF4A_target_genes, "RFS_event")]     
  
  HNF4A_targets_expressionTop100 = de[, c(levels(droplevels(Top100HNF4Atargets)), "RFS_event")]     
  
  write.csv(file = "/media/user/Edison1/GeneRank_Ruth/Noa Deep learning/DeepHNF4A/DeepHNF4A/forPaper_ttest/HNF4A_targets_expression.csv"
            , HNF4A_targets_expression
            , row.names = F)
  
  write.csv(file = "/media/user/Edison1/GeneRank_Ruth/Noa Deep learning/DeepHNF4A/DeepHNF4A/forPaper_ttest/HNF4A_targets_expression100.csv"
            , HNF4A_targets_expressionTop100
            , row.names = F)
  
}
  
savefiles <- function(gene_expression, probe_expression, gene.data){
  setwd("/media/user/Edison1/GeneRank_Ruth/Rahul_analysis/jupyter-ttest/Step -2 Probe2Gene")
  write.csv(gene_expression, "geneLevel_expression.csv")
  write.csv(probe_expression, "probeLevel_expression.csv")
  
  # Change directory and save the gene names and probe names
  setwd("/media/user/Edison1/GeneRank_Ruth/Data-Noa-GeneBased")
  write.table(file = "backgroundgenes_13126.txt"
              , names(gene_expression)
              , row.names = F
              , quote = F)
  
  write.table(file = "backgroundprobes_13126.txt"
              , names(probe_expression)
              , row.names = F
              , quote = F)
  
  
  ### Create a new csv file containg probename, genename, t-test and pvalue 
  setwd("/media/user/Edison1/GeneRank_Ruth/Rahul_analysis/jupyter-ttest/Step-2.1 Probe2Gene/")
  write.csv(gene.data[order(gene.data$p_value),]
            , "List.probe.gene.5y.csv")
  
  
  ## Save top 1% genes
  setwd("/media/user/Edison1/GeneRank_Ruth/Rahul_analysis/jupyter-ttest/Step-3 top131(123mappedCPDB)Genes/")
  write.table(gene.data[order(gene.data$p_value),][,1][1:131]
              , "Top131genes.txt"
              , quote = F
              , row.names= F)
  ## Save top 1% probes
  write.table(gene.data[order(gene.data$p_value),][,2][1:131]
              , "Top131probes.txt"
              , quote = F
              , row.names= F)
  
  ## Save background gene list
  write.table(gene.data[order(gene.data$p_value),][,1]
              , "backgroundgenes_13126.txt"
              , quote = F
              , row.names= F)
  
  ## Save background probe list
  write.table(gene.data[order(gene.data$p_value),][,2]
              , "backgroundprobes_13126.txt"
              , quote = F
              , row.names= F)
  
}

arrangeAndsaveData <- function(df.new){
  df.new.Gene = df.new[, -c(1, 3, 4)]
  df.new.Probe = df.new[, -c(2, 3, 4)]
  
  gene_expression = as.data.frame(t(df.new.Gene[, -1]))
  probe_expression = as.data.frame(t(df.new.Probe[, -1]))
  
  names(gene_expression) = df.new.Gene[, 1]
  names(probe_expression) = df.new.Probe[, 1]
  
  # save the files
  savefiles(gene_expression, probe_expression, gene.data)
  
  return(gene_expression)
  
}

# TODO: A funciton to calculate t-test 
# Perform a t-test on every probe expression in both 
t.test_function <- function( group1, group2) {
  p_value_list = NULL
  t_value_list = NULL
  for(i in 1:dim(group1)[1]){
    p_value_list[i] = t.test(group1[i, -1], group2[i, -1],paired = F)$p.value
    t_value_list[i] = t.test(group1[i, -1], group2[i, -1], paired = F)$statistic
  }
  
  statistics = list(p = p_value_list, t = t_value_list)
  return(statistics)
}

limma_analysis <- function(group1, group2, dim1, dim2){
  
  data = merge(group1
               , group2
               , by = "row.names"
  )
  row.names(data) = data$Row.names
  data$Row.names = NULL  
  
  # The dataframe needs to be converted to an ExpressionSet object to be used with limma
  rma95<-new("ExpressionSet", exprs=as.matrix(data))
  
  
  fac <- factor(rep(1:2,c(dim1, dim2)))
  
  
  fit <- lmFit(rma95, design=model.matrix(~ fac))
  colnames(coef(fit))
  
  fit <- eBayes(fit)
  tt <- topTable(fit, coef=2)
  write.csv(tt, "../limma-test.csv"
            , stringsAsFactors = FALSE)
  
  # view first few lines of the top genes
  # topTable(fit, coef=2, number=Inf, sort.by="none")
  
  return(tt)
  
}

probeSummation_limma <- function(genes, limma_result){
  # 
  # For probe selection we used an approach where we selected the probe that 
  # is the best separator between recurred and non-recurred patients, from all possible 
  # probes that maps to a gene
  # 
  
  gene.level.data=NULL;
  genes.to.process= unique(genes$Gene);
  for (this.gene in genes.to.process) {
    #this.gene = "ATG3"
    this.probes= genes$Probe[genes$Gene==this.gene];
    #this.probes
    this.data.slice= limma_result[limma_result[,"Probe"] %in% this.probes,];
    #this.data.slice
    this.data.slice.means =   this.data.slice[which(this.data.slice$p_value == min(this.data.slice$p_value)), ]
    #this.data.slice.means  
    gene.level.data=rbind(gene.level.data,this.data.slice.means)
  }
  
  dim(gene.level.data)
  gene.level.data[1:5,]
  ### Make sure that the number of probes ie the 
  ### length of genes.to.process and gene.level.data$Probe are the same
  length(gene.level.data$Probe)
  length(genes.to.process)
  # We were able to map 19878 unique probes to 13101 unique genes
  
  # Combine the unique gene names to their corresponsing probe name
  gene.data=cbind(genes.to.process, gene.level.data)
  rownames(gene.data) <- 1:dim(gene.data)[1]
  # Give column names
  colnames(gene.data) <- c('Gene', 'Probe','t_test', 'p_value')
  #rm (gene.level.data, genes, ig_r_5, this.data.slice, this.data.slice.means, this.gene, this.probes, genes.to.process);
  write.csv(gene.data, file = "/media/user/Edison1/GeneRank_Ruth/Data-Noa-GeneBased/ProbeSummarized2gene_usingTtest.csv")
  #write.csv(list_na, file = "NA probes_of the expr data.csv")
  
  # ## I was able to retrieve 14114 genes and probes from the gprofiler
  return(gene.data)
}

Figure_one <- function(gene_expression, gene.data){
  top10_gene= gene.data$Gene[1:10]
  topIGgenes = gene_expression[, c(as.character(top10_gene), "RFS_event")]
  
  tmp.df = melt(topIGgenes,id.vars = "RFS_event")  
  tmp.df[1:12,]
  #tmp.df = tmp.df[tmp.df$variable=="CACYBPP2",]
  tmp.df$RFS_event[tmp.df$RFS_event==1] = "red" 
  tmp.df$RFS_event[tmp.df$RFS_event==0] = "blue"
  #tmp.df$variable = as.character(tmp.df$variable)
  summary(tmp.df)
  
  # DATA #####
  df = abs(gene.data$t_test)
  # Fit the normal mixture(s)
  mixmdl <- normalmixEM(df, k = 3)
  
  #################Plot1 ########################
  # GLOBAL THEME AND GLOBAL AESTHETICS
  old <- theme_set(theme_bw() +
                     theme(text = element_text(size=12),
                           axis.title = element_text(size = 14, face="bold"),
                           title = element_text(face = "bold")
                     )
  )
  update_geom_defaults("point", list(size = 2.5))
  
  ###
  # FUNCTIONS ###
  ###
  # Motivated by http://tinyheero.github.io/2016/01/03/gmm-em.html
  #' Plot a Mixture Component
  #' 
  #' @param x Input data.
  #' @param mu Mean of component.
  #' @param sigma Standard of component.
  #' @param lam Mixture weight of component.
  plot_mix_comps <- function(x, mu, sigma, lam) {
    lam * dnorm(x, mu, sigma)
  }
  
  
  # Motivated by http://tinyheero.github.io/2016/01/03/gmm-em.html
  p1 <- data.frame(x = df) %>%
    ggplot(aes(x=x)) +
    geom_histogram(aes(x, ..density..), binwidth = 1, colour = "black", 
                   fill = "white") +
    stat_function(fun = plot_mix_comps, aes(colour="1"),
                  args = list(mixmdl$mu[1], mixmdl$sigma[1], 
                              lam = mixmdl$lambda[1]), lwd = 1.5) +
    stat_function(fun = plot_mix_comps, aes(colour="2"),
                  args = list(mixmdl$mu[2], mixmdl$sigma[2], 
                              lam = mixmdl$lambda[2]), lwd = 1.5) +
    stat_function(fun = plot_mix_comps, aes(colour="3"),
                  args = list(mixmdl$mu[3], mixmdl$sigma[3], 
                              lam = mixmdl$lambda[3]), lwd = 1.5) +
    scale_colour_manual("Component", values = c("1"="green", "2"="blue", "3"="red")) +
    ylab("Density") +
    xlab("abs(t test statistic)") +
    theme(legend.justification = c(0.98, 0.98)
          , legend.position = c(0.98, 0.98)
          , legend.text=element_text(size=10)
          , plot.title = element_text(size = 12, face = "bold")
          , legend.title=element_text(size=10) 
          , axis.title.x = element_text(size=14, face="bold")
          , axis.title.y = element_text(size=14, face="bold")
          , axis.text.x = element_text(size = 15
                                       , color = "black"
                                       , angle = 90, hjust = 1)
          , axis.text.y = element_text(size = 15
                                       , color = "black"))
  
  
  ################Plot 2 #######################
  p2 <- ggplot(tmp.df
               , aes(x = variable
                     , y = log10(value)
                     , fill= as.factor(tmp.df$RFS_event)
               )
  )  + xlab(label = 'Genes') + ylab(label = "ln(Expression)") +
    geom_boxplot(position=position_dodge(1)) +
    
    #stat_compare_means(
    #  aes(label = paste0(..p.format..))
    #)  +
    scale_fill_discrete(name = "RFS"
                        #, values = c("blue", "red")
                        #, guide = guide_legend(reverse = TRUE)
                        , labels=c("Recurred", "Non Recurred")
    )   +
    theme(legend.justification = c(0.55, 0.03)
          , legend.position = c(0.55, 0.03)
          , legend.text=element_text(size=7)
          , plot.title = element_text(size = 7, face = "bold")
          , legend.title=element_text(size=7) 
          , axis.title.x = element_text(size=14, face="bold")
          , axis.title.y = element_text(size=14, face="bold")
          , axis.text.x = element_text(size = 14
                                       , color = "black"
                                       , angle = 90, hjust = 1)
          , axis.text.y = element_text(size = 14
                                       , color = "black"))
  
  
  
  ####### save the plot ######
  
  p = cowplot::plot_grid(p1, p2, labels = "auto" )
  save_plot("Figure1.pdf", p, ncol = 2)
  
}

#### mixdistribution plot t-test
analyze_mixdist <- function (gene.data){
  wait = abs(gene.data$t_test)
  set.seed(123)
  mixmdl = normalmixEM(wait, k = 3) #, lambda = .5, sigma = 1)
  
  png("DensityPlot_mixedmodel.png", width = 1000, height = 600)
  plot(mixmdl, density = T, xlab2 = "absolute t statistic"
       #,main2 = main_m
       , whichplots = 2
       #, cex = 1.5
       #, cex.names = 1.5
       , cex.axis = 1.3
       , font = 2
       , font.lab = 2
       , cex.lab = 1.3
       , breaks = 30) 
  dev.off()
}

### Function to save files
save_expression_sorted_files <- function(expr, clin, filename1, filename2){
  setwd("/media/user/Edison1/GeneRank_Ruth/Rahul_analysis/jupyter-ttest/Step 6 - kmplotter/")
  write.table(expr
              , filename1
              , sep = "\t"
              , quote = F
              , row.names = F)
  write.table(clin
              , filename2
              , sep = "\t"
              , quote = F
              , row.names = F)
}

### Functions to write file
sortAffy <- function(expr, clin, type){
  names(expr)[1] = "AffyID"
  names(clin)[2] = "AffyID"
  clin = clin[order(clin$AffyID), ]
  expr = expr[order(expr$AffyID), ]
  filename1 = "@GEO expression data_sorted.txt"
  filename2 = "@GEO clinical data_sorted.txt"
  save_expression_sorted_files(expr, clin, filename1, filename2)
  
  return(list(expr, clin))
  
}

checkData = function(expr, clin){
  
  affyid_expr = as.character(expr[[1]]);
  affyid_clin = as.character(clin[[1]]);
  
  # check number of entries
  if(length(affyid_expr) != length(affyid_clin)){
    stop("The dimensions of data is not equal.");
  }
  
  for(i in 1:length(affyid_expr)){
    if(affyid_expr[[i]] != affyid_clin[[i]]){
      stop( paste("STOP: The ", i, "th ID of data is different.", sep="") )
    }
  }
  
}

auto_cutoff = function(row, gene_db, time_index, event_index){
  
  ordered_row = order(row);
  q1 = round(length(ordered_row)*0.25);
  q3 = round(length(ordered_row)*0.75);
  
  # m = sortedrow[round(i)];
  surv = Surv(gene_db[,time_index], gene_db[,event_index]);
  
  p_values = vector(mode="numeric", length = q3-q1+1)
  min_i = 0
  min_pvalue=1
  
  for(i in q1:q3){
    
    gene_expr = vector(mode="numeric", length=length(row))
    gene_expr[ordered_row[i:length(ordered_row)]] = 1
    
    cox = summary(coxph(surv ~ gene_expr))
    
    pvalue = cox$sctest['pvalue']
    
    p_values[i-q1+1] = pvalue
    
    if(pvalue < min_pvalue){
      min_pvalue = pvalue
      min_i = i
    }
    
  }
  
  gene_expr = vector(mode="numeric", length=length(row))
  gene_expr[ordered_row[min_i:length(ordered_row)]] = 1
  
  
  # overwrite m (median) and gene_expr
  m = row[ordered_row[min_i]]
  
  m
}

mySurvplot = function(surv, gene_expr, xlab="Time (years)", ylab="Probability", snames = c('low', 'high'), stitle = "Expression", hr.pos=NA){
  survplot(surv ~ gene_expr, xlab=xlab, ylab=ylab, snames = snames, stitle = stitle, hr.pos=hr.pos);
  
  cox = summary(coxph(surv ~ gene_expr))
  
  pvalue=cox$sctest['pvalue'];
  hr = round(cox$conf.int[1],2)
  hr_left = round(cox$conf.int[3],2)
  hr_right = round(cox$conf.int[4],2)
  
  conf_int =  paste(" (", hr_left, " - ", hr_right, ")", sep=""); 
  
  txt = paste("HR = ", hr, conf_int, "\nlogrank P = ", signif(pvalue, 2), sep="")
  text(grconvertX(0.98, "npc"), grconvertY(.97, "npc"),
       labels=txt,
       adj=c(1, 1))
  
  list(pvalue, hr, hr_left, hr_right)
}

createDirectory = function(base){
  i="";
  while(file.exists(paste(base, i, sep=""))){
    if(i==""){
      i=1;
    }else{
      i=i+1;
    }
  }
  toDir = paste(base, i, sep="")
  dir.create(toDir)
  
  toDir
}

getCutoff = function(quartile, median_row, manual_cutoff="false", verbose=FALSE){
  
  sortedrow=order(median_row);
  minValue = median_row[sortedrow[1]]
  maxValue = median_row[sortedrow[length(sortedrow)]]
  
  # if manual_cutoff is true, then the user can specify a discrete Cutoff value
  # 
  if(manual_cutoff=="true"){
    m = as.numeric(quartile);
    indices = which(median_row>m)
  }else{
    quartile = as.numeric(quartile);
    
    if(is.na(quartile) || !is.numeric(quartile)){
      stop("The quartile parameter isn't numeric.")
    }else if(quartile<5){
      quartile = 5
    } else if(quartile > 95){
      quartile = 95
    }
    
    i=length(sortedrow)*quartile/100;
    m=median_row[sortedrow[round(i)]];
    
    indices = which(m<median_row)
    #sortedrow[round(i):length(sortedrow)]
  }
  
  if(verbose){
    print(m)
    print(minValue)
    print(maxValue)
  }
  
  list(m, minValue, maxValue, indices)
}

getParameter = function(c_args, id){
  res = ""
  for(i in 1:length(c_args)){
    
    tmp = strsplit(c_args[i], "=")
    filter_id = tmp[[1]][1]
    value = tmp[[1]][2]
    
    if(filter_id == paste("-", id, sep="")){
      res = value
      break
    }
    
  }
  
  res
  
}

kmplot = function(expr, clin, event_index=3
                  , time_index=4, affyid="ESR1"
                  , auto_cutoff="true", quartile=50
                  , file_name){
  # checks the input: if the expression data and clinical data don't match, the script will fail.
  checkData(expr, clin);
  pvalueVector = list()
  gene_names  = names(expr)[-1] # omit the first affyID
  hrVector = list()
  survival_data = cbind(as.numeric(clin[[time_index]]), as.numeric(clin[[event_index]]));
  
  toDir = createDirectory("res");
  resTable=rbind();
  
  index_arr = 2:dim(expr)[2];
  if(affyid != "230772_at"){
    index = which(colnames(expr) == paste("X", affyid, sep=""));
    if(length(index) > 0){
      index_arr = c(index);
    }
  }
  for(j in 1:length(index_arr)){
    i=index_arr[j]
    print(paste("Processing...", colnames(expr)[i]));
    
    row = as.numeric(expr[[i]]);
    
    # --------------------- CUTOFF ----------------------
    if(auto_cutoff == "true"){
      m = auto_cutoff(row, survival_data, 1, 2)
    }else{
      # calculates lower quartile, median, or upper quartile
      tmp = getCutoff(quartile, row)
      
      m = tmp[[1]]
      minValue = tmp[[2]]
      maxValue = tmp[[3]]
      indices = tmp[[4]]
      
    }
    
    # gene_expr consists 1 if m smaller then the gene expression value,
    # 0 if m bigger then the gene expression value
    gene_expr=vector(mode="numeric", length=length(row))
    gene_expr[which(m < row)] = 1
    
    # --------------------- KMplot ----------------------
    tryCatch({
      # draws the KM plot into a png file
      
      png(paste(toDir, "/", colnames(expr)[i], ".png", sep=""));
      
      # Surv(time, event)
      surv<-Surv(survival_data[,1], survival_data[,2]);
      
      res = mySurvplot(surv, gene_expr)
      pvalue = res[[1]];
      hr = res[[2]]
      hr_left = res[[3]]
      hr_right = res[[4]]
      
      resTable = rbind(resTable, c(pvalue, hr, hr_left, hr_right));
      
      dev.off();
      
    }, interrupt = function(ex){
      cat("Interrupt during the KM draw");
      print(ex);
      dim(gene_expr);
    }, error = function(ex){
      cat("Error during the KM draw");
      print(ex);
      dim(gene_expr);
    }
    );
    
    pvalueVector[[j]] = pvalue
    hrVector[[j]] = hr
    print(j)
    
  }
  print(unlist(pvalueVector))
  print(unlist(hrVector))
  mystat_table = data.frame(gene = gene_names , pvalue= unlist(pvalueVector), hr = unlist(hrVector))
  write.csv(mystat_table, "mystat_table.csv")
  
}

clique_calculation <- function(expr){
  library(igraph)
  cor_mat = cor(expr[,-1])
  cor_g <- graph_from_adjacency_matrix(cor_mat
                                       , mode='undirected'
                                       , weighted = 'correlation')
  cor_edge_list <- as_data_frame(cor_g, 'edges')
  only_sig <- cor_edge_list[abs(cor_edge_list$correlation) > 0.45, ]
  new_g  <- graph_from_data_frame(only_sig, F)
  
  png(filename = "Clique.png", 
      width = 1200, height= 1000)
  plot(new_g, vertex.size=10)
  dev.off()
  
  clique_num(new_g)
  cliques(new_g, min=12)
  
  largest_cliques(new_g)
  
  closeness(new_g)
  
}

generate_clique_network = function(expr){
  
  ## Create correlation map ##
  #### create correlation map ####
  size = 27
  d.selected.names = expr[, -1]
  cormat = matrix(nrow=size,ncol=27)
  for (i in 1:size) {
    for (j in 1:size) {
      cormat[i,j] = cor.test(d.selected.names[,i], d.selected.names[,j])$estimate
    }
  }
  colnames(cormat) = colnames(d.selected.names)
  rownames(cormat) = colnames(d.selected.names)
  
  library(pheatmap)
  png(filename = "Correlation-heatmap.png", 
      width = 1200, height= 1000)
  pheatmap(data.matrix(cormat)
           , cellwidth=30, cellheight=30,
           , fontsize_row=25
           , fontsize_col=25
           , color = rainbow(size, s = 1, v = 1, 
                             start = 0, 
                             end = max(1, size - 1)/size, 
                             alpha = 1))
  dev.off()
  
  ## Identify the cliques ##
  d.selected = t(expr[, -1])
  
  all.cors <- NULL
  for (x.name in rownames(d.selected)) {
    for (y.name in rownames(d.selected)) {
      if (x.name >= y.name) {
        next;
      }
      #print(c(x.name,y.name))
      x <- d.selected[x.name,]
      y <- d.selected[y.name,]
      this.cor <- cor(x=x,y=y,use = "pair",method="spearman")
      this.cor.test=cor.test(x=x,y=y,use = "pair",method="spearman")
      if(this.cor.test$p.value>0.05) {
        this.cor=-999
      }
      #all.cors<-rbind(all.cors,c(x.name,y.name,jetset[x.name,"symbol"],jetset[y.name,"symbol"],this.cor))
      all.cors<-rbind(all.cors,c(x.name,y.name, this.cor))
    }
  }
  
  selected.probes <- all.cors[as.numeric(all.cors[,3])>=0.60, 1]
  selected.probes <- c(selected.probes,all.cors[as.numeric(all.cors[,3])>=0.60, 2])
  
  selected.probes <- unique(selected.probes)
  final.graph <- all.cors[all.cors[,1] %in% selected.probes &  all.cors[,2] %in% selected.probes,]
  
}


#### Main #####

# Input files
# TODO: Input data probe based expression
combined = read.csv("/media/user/Edison1/GeneRank_Ruth/Data-Noa-GeneBased/OriginalData and clinical/CombinedExpression.csv")

# Remove clnical data
combined_filetered = combined[, c(7:22221,22223)]

# TODO: Filter recurred and non-recurred event
combined_event1 = combined_filetered[combined_filetered$RFS_event==1, -22216]#[1:10, 1:10]
combined_event0 = combined_filetered[combined_filetered$RFS_event==0, -22216]#[1:12, 1:10]

# TODO: transpose the dataframe so that each row is a probe
# and we do t-test for each pair for both groups (recurred and non-recurred)
group1 = t(combined_event1)
group2 = t(combined_event0)

# TODO: Perform t-test
m = t.test_function(group1, group2)

t.test.result = data.frame(Probe = row.names(group1), t_test = m$t, p_value = m$p)
#t.test.result = t.test.result[which(t.test.result$p_value < 0.05),]
t.test.result = t.test.result[order(t.test.result$p_value),]

removeXcharacter <- function(df.result){
  # Remove 'X' character from the probe name
  df.result$Probe = as.character(df.result$Probe)
  for(i in 1:dim(df.result)[1]){
    new_probe = paste0(strsplit(as.character(df.result$Probe[i]), '')[[1]][-1], collapse = '')
    df.result$Probe[i] = new_probe
  }

  row.names(df.result) = NULL
  return(df.result)
}

t.test.result = removeXcharacter(t.test.result)

write.csv(t.test.result
          , file = "/media/user/Edison1/GeneRank_Ruth/Data-Noa-GeneBased/t.test.result.csv"
          , row.names = F)


# Use the above probe list and find their gene names
genes <- read.csv("/media/user/Edison1/GeneRank_Ruth/Data-Noa-GeneBased/gprofiler_results_ttest.csv", header=T, stringsAsFactors=FALSE)
colnames(genes) <- c('Probe','Gene')

# Find all duplicated probes index
dup_index = which(duplicated(genes$Probe) == TRUE)

# remove the duplicated names
if(length(dup_index) > 0){
  genes = genes[-dup_index,]
}

# Remove all NA values meaning probes that are not mapped to any gene
genes = removeNAvalues(genes)

# Make the probe names to case
t.test.result$Probe = tolower(t.test.result$Probe)
genes$Probe = tolower(genes$Probe)

# After removing NA values from 22215 probes we have 19878 probes remaining
# We will merge this df with ig_r_5 containg calculated IG value so that
# both file will have the same number of probes

# Select only those unique genes which has a probe value and
# probes which doent has NA values
t.test.result = t.test.result[which(t.test.result$Probe %in% genes$Probe == TRUE),]

#> dim(t.test.result)
#[1] 15014     3
#> dim(genes)
#[1] 15014     2

# Probe summation of the probes to get uniquely mapped gene values
gene.data = probeSummation(genes, t.test.result)

write.csv(gene.data, "gene.data.csv")


# ----- Limma validation #
limma_analysis(group1[, 1:689]
               , group2[, 1:830])

limma_result = read.csv("/media/user/Edison1/GeneRank_Ruth/Rahul_analysis/jupyter-ttest/limm_analysis/limma-test.csv", stringsAsFactors = F)
# Remove 'X' character from the probe name
limma_result$X = as.character(limma_result)
for(i in 1:dim(limma_result)[1]){
  new_probe = paste0(strsplit(as.character(limma_result$X), '')[[1]][-1], collapse = '')
  limma_result$X[i] = new_probe
}

row.names(limma_result) = NULL
limma_result = limma_result[, c(1, 4, 5, 6)]

colnames(limma_result) = c("Probe", "t_test", "p_value", "adj.P.Val")
# After removing NA values from 22215 probes we have 19878 probes remaining
# We will merge this df with ig_r_5 containg calculated IG value so that
# both file will have the same number of probes

# Select only those unique genes which has a probe value and
# probes which doent has NA values
t.test.resultk = limma_result[which(unique(limma_result$Probe) %in% genes$Probe == TRUE),]

t.test.resultk = t.test.result[which(t.test.result$Probe %in% genes$Probe == TRUE),]


#> dim(t.test.result)
#[1] 15014     3
#> dim(genes)
#[1] 15014     2

# Probe summation of the probes to get uniquely mapped gene values
gene.data = probeSummation(genes, t.test.result)




# probe summation
probeSummation_limma <- probeSummation_limma(genes, limma_result)
######### remove unused data ##############
rm(combined)
rm(combined_event0)
rm(combined_event1)
rm(combined_filetered)
rm(group1)
rm(group2)

# Read expression data of breast cancer patients
probe_expr_transpose = readExpressiondataBreastCancer(combined)

dim(probe_expr_transpose)
#[1] 22216  1520

dim(gene.data)
#[1] 10778     4

probe_expr_transpose$Probe = tolower(probe_expr_transpose$Probe)
probe_expr_transpose[1:5, c("Probe","GSM107072", "GSM107073")]
#Probe GSM107072 GSM107073
#1 1007_S_AT      4033      4905
#2   1053_AT       229       482
#3    117_AT       332       375
#4    121_AT      1669      2161
#5 1255_G_AT        99       103

### We have one matrix with probe name and expression 
#### another matrix with the probename, gene name, infogain and GR

gene.data = as.data.frame(gene.data)

head(gene.data)
dim(gene.data)
# Gene       Probe   t_test      p_value
# 1     ATG3 221492_s_at 15.76298 4.673388e-50
# 2 HSP90AA1 210211_s_at 15.30039 1.433973e-48
# 3 HNRNPCP2 200751_s_at 15.39926 6.339125e-48
# 4     CBX3 201091_s_at 15.38553 1.134624e-47
# 5    ACTR2 200729_s_at 15.29278 2.289181e-47
# 6    PSMA4   203396_at 15.24054 4.428939e-47



df.new = merge(gene.data, probe_expr_transpose, by.x = "Probe", by.y = "Probe")
df.new[1:5, 1:8]
dim(df.new)



# Save gene/probe expressions, save top 1% genes for downstream analysis
gene_expression = arrangeAndsaveData(df.new)
gene_expression$Patient = row.names(gene_expression)
#gene_expression = read.csv("/media/user/Edison1/GeneRank_Ruth/Rahul_analysis/jupyter-ttest/Step -2 Probe2Gene/geneLevel_expression.csv")

# clinical data 
GEOclin = read.csv("/media/user/Edison1/GeneRank_Ruth/Rahul_analysis/jupyter-ttest/Step-3 top131(123mappedCPDB)Genes/ProbeExpressionClinical.csv")
clin = GEOclin[, c("Patient", "RFS_event")]

gene_expression = merge(gene_expression
                         , clin
                         , by.x = "Patient"
                         , by.y = "Patient")

# Figure 1
Figure_one(gene_expression, gene.data )

# Figure 2 b
read_centrality = read.table("/media/user/Edison1/GeneRank_Ruth/Publication/Final drafts/Supplementary Files/centrality.txt"
                             , nrows=10)
plot(read_centrality$V2[1:8] 
     , type = "o"
     #, names.arg = as.character(read_centrality$V1[1:8])
     , cex.names = 4
     , font = 2
     , xlab = "Genes"
     , ylab = "Degree"
     )
#### Microarray data
GEOexpr = gene_expression
GEOexpr = subset(GEOexpr, select = -c(RFS_event))
rm(GEOexpr)
rm(GEOclin)

## FOR GEO data:
expression_clin_sorted = sortAffy(GEOexpr, GEOclin, type = "GEO")
expr  = expression_clin_sorted[[1]]
expr = subset(expr, select = c("AffyID", "ESR1", as.character(HNF4A_genes_30$Gene)))

clin  = expression_clin_sorted[[2]]

### Get only the survival and event time
#setwd("/media/user/Edison1/GeneRank_Ruth/Rahul_analysis/jupyter-ttest/Step 6 - kmplotter/")
#expr = read.table("@GEO expression data_sorted.txt"
#                  , header = T)

#expr = subset(expr, select = c("AffyID", "ESR1", as.character(HNF4A_genes_30$Gene)))
#clin = read.table("@GEO clinical data_sorted.txt"
#                  , header = T
#                  , sep = "\t")

###### ER positive and negative
clin = clin[, c("AffyID", "RFS_event", "RFS_time", "ER.status")]
names(clin) = c("AffyId", "Survival event", "Survival time", "ER.status")

clinERp = clin[which(clin$ER.status ==1),]
clinERn = clin[which(clin$ER.status ==0),]

# Get only those patient IDs that belongs to ER+ and ER- patients
# Subset the corresponding gene expression for the ER+ and ER- patients
bg2011missingFromBeg <- setdiff(x= as.character(expr$AffyID), y= as.character(clinERn$AffyId))
expr_ERn <- expr[!expr$AffyID %in% bg2011missingFromBeg, ]

bg2011missingFromBeg <- setdiff(x= as.character(expr$AffyID), y= as.character(clinERp$AffyId))
expr_ERp <- expr[!expr$AffyID %in% bg2011missingFromBeg, ]


names(clin[1]) = "AffyID"
names(expr[1]) = "AffyID"
names(expr_ERn[1]) = "AffyID"
names(expr_ERp[1]) = "AffyID"

# validation data
kmplot(expr
       , clin
       , event_index=2
       , time_index=3
       #,  affyid="ESR1"
       , auto_cutoff="true"
       , file_name = "All patients")

kmplot(expr_ERp
       , clinERp
       , event_index=2
       , time_index=3
       #,  affyid="ESR1"
       , auto_cutoff="true"
       , file_name = "ERp")

kmplot(expr_ERn
       , clinERn
       , event_index=2
       , time_index=3
       #,  affyid="ESR1"
       , auto_cutoff="true"
       , file_name = "ERn")



demo = getParameter(c_args, "demo");


####### correlation map of HNF4A targets #########
clique_calculation(expr)

# generate_clique_network #
generate_clique_network(expr)


### ssave the probe names of HNF4A targets ##
HNF4A_targets_probes= gene.data[gene.data$Gene %in% c("HNF4A", as.character(HNF4A_genes_30$Gene)), ]
write.csv(HNF4A_targets_probes
          , file = "/media/user/Edison1/GeneRank_Ruth/Rahul_analysis/jupyter-ttest/Step-2.1 Probe2Gene/HNF4A_targets_probes.csv"
          )








