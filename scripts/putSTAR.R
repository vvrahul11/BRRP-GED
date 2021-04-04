###############################################################
############ Put stat in the gene name ########################
put_Star <- function(gene_vector){
  
  IG_high = as.character(read.table("Top131genes.txt", header =F )$V1)
  for(i in 1:length(gene_vector)){
    print(gene_vector[i])
    if(gene_vector[i] %in% IG_high){
      next
    }
    else{
      gene_vector[i] = paste("*"
                             , gene_vector[i]
                             , sep = "")
      
    }
    
  }
  return(gene_vector)
  
}
### read file in ###
networkfile = read.csv("CPDB_inducedModules131.csv", sep = "\t")

## extract gene name column and call the function

geneA = as.character(networkfile[,2])
geneA.A = put_Star(geneA)
geneB = as.character(networkfile[,7])
geneB.B = put_Star(geneB)

#### Rename the gene names column with new vector names
networkfile[,2] = geneA.A
networkfile[,7] = geneB.B

write.table(networkfile
          , "HNF4AInducedNetwork_renamed.txt"
          , quote = F, row.names = F, sep = "\t")

write.csv(networkfile
          , "HNF4AInducedNetwork_renamed.csv")



