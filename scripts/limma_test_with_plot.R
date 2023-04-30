# http://genomicsclass.github.io/book/pages/using_limma.html
# http://www.sthda.com/english/wiki/unpaired-two-samples-t-test-in-r

limma_analysis <- function(group1, group2){
  #browser()
  library(SpikeInSubset)
  library(genefilter)
  
  data = merge(group1
               , group2
               , by = "row.names"
  )
  row.names(data) = data$Row.names
  data$Row.names = NULL  
  
  # The dataframe needs to be converted to an ExpressionSet object to be used with limma
  rma95<-new("ExpressionSet", exprs=as.matrix(data))
  fac <- factor(rep(1:2,c(ncol(group1), ncol(group2))))
  
  # ttest analysis
  rtt <- rowttests(exprs(rma95),fac)
  mask <- with(rtt, abs(dm) < .2 & p.value < .01)
  spike <- rownames(rma95) %in% colnames(pData(rma95))
  cols <- ifelse(mask,"red",ifelse(spike,"dodgerblue","black"))
  
  with(rtt, plot(-dm, -log10(p.value), cex=.8, pch=16,
                 xlim=c(-1,1), ylim=c(0,5),
                 xlab="difference in means",
                 col=cols))
  abline(h=2,v=c(-.2,.2), lty=2)
  
  
  rtt$s <- apply(exprs(rma95), 1, function(row) sqrt(.5 * (var(row[1:3]) + var(row[4:6]))))
  with(rtt, plot(s, -log10(p.value), cex=.8, pch=16,
                 log="x",xlab="estimate of standard deviation",
                 col=cols))
  
  
  
  # limma analysis
  fit <- lmFit(rma95, design=model.matrix(~ fac))
  colnames(coef(fit))
  
  fit <- eBayes(fit)
  tt <- topTable(fit, coef=2)
  write.csv(tt, "Output/limma-test.csv"
            , stringsAsFactors = FALSE)
  
  # view first few lines of the top genes
  topTable(fit, coef=2, number=Inf, sort.by="none")
  
  limmares <- data.frame(dm=coef(fit)[,"fac2"], p.value=fit$p.value[,"fac2"])
  with(limmares, plot(dm, -log10(p.value),cex=.8, pch=16,
                      col=cols,xlab="difference in means",
                      xlim=c(-5000, 8000), ylim=c(0,90)))
  abline(h=2,v=c(-.2,.2), lty=2)
  return(tt)
  
}


limma_analysis(group1[, 1:689]
               , group2[, 1:830])