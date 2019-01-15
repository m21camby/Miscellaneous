# Written by SJK (15 Jan 2019) 
# This file is about short function and ggplot for Seurat clustering analysis
# After running Seurat and findcluster, I want to know how many cells and median of UMI per cluster and show by barplot
# This file is written under Seurat(2.3.4), reshape2(1.4.3), ggplot2(3.1.0)
# In the future when Seurat version 3.0.0 comes, The code might need modification

#################################################

# Function of creating number of cell and median UMIs data frame
cluster_analysis <- function(Seurat, number_clusters, DF){
  iterations = number_clusters+1
  variables = 2
  output <- matrix(ncol=variables, nrow=iterations)
  for(i in c(0:number_clusters)){
    foo1 <- Seurat@ident[Seurat@ident == i] # extract C.B. by each cluster 
    foo2 <- rownames(as.data.frame(foo1)) # extract C.B. by each cluster
    foo3 <- DF[, foo2] # extract DF by each cluster C.B
    output[i+1,1] <- length(foo1) # 1st: how many cells in each cluster
    output[i+1,2] <- round(median(colSums(foo3))) # 2nd: median of UMIs per cluster 
  }
  output <- data.frame(output) # matrix to data frame 
  colnames(output) <- c("cell_number", "med_UMIs") # change col names
  rownames(output) <- paste0("cluster", "_", seq(from = 0, to = number_clusters)) # add row names by cluster
  return(output)
}

# Run function
foo <- cluster_analysis(ds079_200_PC10.SO, 11, ds079) # input as Seurat object, the number of clusters(e.g. 0 ~ 11), and orignal Data frame from Drop-seq pipeline

# data modfication for ggplot
foo$group <- rownames(foo)
foo$group <- factor(foo$group, rownames(foo))
foo2 <- melt(foo)

# bar plot
ggplot(data=foo2, aes(x=group, y=value, fill=variable)) +
  geom_bar(stat="identity", position=position_dodge()) + 
  geom_text(aes(label=value), position=position_dodge(width=0.9), vjust=-0.25) + 
  scale_fill_manual(values=c('#990033','#3399FF')) +
  theme(axis.title = element_blank(), axis.text = element_text(face = "bold"), legend.position = c(0.8, 0.8), legend.title = element_blank())


