# This file is for after running DESeq2 and various ways to do GO and GSEA
# Brief default running of the analysis

library(DESeq2)
library(clusterProfiler)
library(topGO)
library(gage)
library(org.Mm.eg.db)


# ************************************* #
# 1. Loading Digital gene expression matrix
# ************************************* #

DGEm <- readRDS("/data/home/R_Scripts/DGEm_final_stranded.rda")



# ************************************* #
# 2. DESeq2
# ************************************* #

sample_info <- data.frame(condition = c("Ctrl","Ctrl","Ctrl","Ctrl",
                                        "KO","KO","KO","KO"), row.names = names(DGEm))

dds <- DESeqDataSetFromMatrix(DGEm, colData = sample_info, design = ~ condition)

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds$condition <- relevel(dds$condition, ref = "Ctrl")
dds.ds <- DESeq(dds)
dds.ds <- estimateSizeFactors(dds.ds)

res <- results(dds.ds, name="condition_KO_vs_Ctrl", test="Wald")
res.df <- as.data.frame(res)
res.df$gene <- rownames(res.df)

# shrinkage apply for GSEA

resLFC <- lfcShrink(dds.ds, coef="condition_KO_vs_Ctrl", type="apeglm")
resLFC.df <- as.data.frame(resLFC)
resLFC.df$gene <- rownames(resLFC.df)



# ************************************* #
# 3. GO by topGO (topGO)
# ************************************* #

# Function to run topGO

runTopGO <- function(DESeq2_results.df, InterestGenes = "all", significant = 0.01){
  
  # gene universe can be all the genes expressed in experiments 
  geneUniverse <- rownames(DESeq2_results.df) 

  # extract gene of interest
  if(InterestGenes == "up"){
  genesOfInterest <- rownames(DESeq2_results.df[which(DESeq2_results.df$log2FoldChange > 0.5 & DESeq2_results.df$padj < significant), ])  
  }
  if(InterestGenes == "down"){
  genesOfInterest <- rownames(DESeq2_results.df[which(DESeq2_results.df$log2FoldChange < -0.5 & DESeq2_results.df$padj < significant), ])
  }
  if(InterestGenes == "all"){
  genesOfInterest <- c(rownames(DESeq2_results.df[which(DESeq2_results.df$log2FoldChange > 0.5 & DESeq2_results.df$padj < significant), ]), rownames(DESeq2_results.df[which(DESeq2_results.df$log2FoldChange < -0.5 & DESeq2_results.df$padj < significant), ]))
  }

  # export Genes of Interest
  genesOfInterest <<- genesOfInterest
  
  # create gene list
  geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
  names(geneList) <- geneUniverse  

  # create topGOdata object run GO analysis
  onts = c( "MF", "BP", "CC" )
  tab <- as.list(onts)
  names(tab) <- onts

  for(i in 1:3){
  sampleGOdata <- new("topGOdata",
                    description = "GO_analysis", 
                    ontology = onts[i],
                    allGenes = geneList, 
                    nodeSize = 10,
                    annot=annFUN.org, mapping="org.Mm.eg.db", ID = "symbol")

  # export GO ID from BP
  if(i == 2){
  allGO <<- genesInTerm(sampleGOdata)
  }
  
  # run tests
  resultTopGO.elim <- runTest(sampleGOdata, algorithm = "elim", statistic = "Fisher" )
  resultTopGO.classic <- runTest(sampleGOdata, algorithm = "classic", statistic = "Fisher" )

  # save to table
   ## look at results
  tab[[i]] <- cbind(data.frame(category = onts[i]), GenTable(sampleGOdata, Fisher.elim = resultTopGO.elim, 
                        Fisher.classic = resultTopGO.classic,
                        orderBy = "Fisher.elim" , topNodes = 20))

}

  # save to dataframe
  topGOResults <- plyr::rbind.fill(tab)
  topGOResults.df <- as.data.frame(topGOResults)

  # calculate gene ratio
  # the ratio be-tween significantly differentially regulated genes inthe respective gene set measured by sequencingand the total number of genes belonging in that GOcategory. 
  # ref: https://www.cell.com/cell-reports/pdfExtended/S2211-1247(18)31242-7 
 
   topGOResults.df$gene_ratio <- topGOResults.df$Significant / topGOResults.df$Annotated
  topGOResults.df$gene_ratio <- round(topGOResults.df$gene_ratio, 4)

  # modification appropriate for plot
  topGOResults.df$Fisher.elim <- as.numeric(topGOResults.df$Fisher.elim)
  topGOResults.df$Fisher.classic <- as.numeric(topGOResults.df$Fisher.classic)
  topGOResults.df$Term <- factor(topGOResults.df$Term, levels = rev(unique(topGOResults.df$Term)))

  return(topGOResults.df)
}



# ************************************* #
# 4. GO by topGOtable (clusterprofile)
# ************************************* #

# set target genes and background genes

# up-regulated (log2FC > 1)
GO_targets1 <- rownames(res.df[which(res.df$padj < .05 & res.df$log2FoldChange > 1), ])

# down-regulated (log2FC < -1)
GO_targets2 <- rownames(res.df[which(res.df$padj < .05 & res.df$log2FoldChange < -1), ])

# background genes
bg_ids <- rownames(dds.ds)[rowSums(counts(dds.ds)) > 0]


# run topGOtable

topgo_BP_KO <- topGOtable(GO_targets1, bg_ids,
                            ontology = "BP",
                            mapping = "org.Mm.eg.db",
                            geneID = "symbol")



# ************************************* #
# 5. GSEA by fgsea (fgsea)
# ************************************* #
# fgsea tutorial: https://stephenturner.github.io/deseq-to-fgsea/

# Creating rnk file
resLFC.df$Gene <- rownames(resLFC.df)
resLFC.df$fcsign <- sign(resLFC.df$log2FoldChange)
resLFC.df$logP=-log10(resLFC.df$pvalue)
resLFC.df$metric= resLFC.df$logP/resLFC.df$fcsign
y <- resLFC.df[,c("Gene", "metric")]

# Load the pathways into a named list
pathways.hallmark <- gmtPathways("/data/home/ReactomePathways.gmt")

# briefly convert mouse genes to human genes (should use biomaRt for conversion)
res_rank <- y[, c("Gene", "metric")]
res_rank$Gene <- toupper(res_rank$Gene)

ranks <- tibble::deframe(res_rank)

# run fgsea
fgseaRes <- fgsea(pathways.hallmark, ranks, minSize=15, maxSize=500, nperm=1000)

# check top pathways
topUp <- fgseaRes %>% 
    filter(ES > 0) %>% 
    top_n(10, wt=-padj)
topDown <- fgseaRes %>% 
    filter(ES < 0) %>% 
    top_n(10, wt=-padj)
topPathways <- bind_rows(topUp, topDown) %>% 
    arrange(-ES)



# ************************************* #
# 6. GSEA by gseKEGG (clusterprofile)
# ************************************* #
# Run this code to run KEGG than Reactome or Hallmark


# As KEGG uses Entrez ID for GSEA, extract genes and Entrez ID

kg.mouse <- kegg.gsets("mouse")
kegg.gs <- kg.mouse$kg.sets[kg.mouse$sigmet.idx]
 
gene.symbol.eg <- select(org.Mm.eg.db, keys=resLFC.df$gene, columns="ENTREZID", keytype="SYMBOL")
colnames(gene.symbol.eg) <- c("gene", "entrezID")

# set metric 
# Here sign of log2FC x -log10(pvalue) is applied (widely used)

resLFC.df$res_pre_rank <- sign(resLFC.df$log2FoldChange) * -log10(resLFC.df$pvalue)

# add Entrez ID to DESeq2 results
kegg.df <- inner_join(resLFC.df, gene.symbol.eg, by = "gene")
kegg.df <- kegg.df[!is.na(kegg.df$entrezID), ]
kegg_sub.df <- kegg.df[,c(8, 9)]

# set input geneList 

# feature 1: numeric vector
 geneList1 = kegg1_sub.df[,1]
# feature 2: named vector
names(geneList1) = as.character(kegg1_sub.df[,2])
# feature 3: decreasing orde
geneList1 = sort(geneList1, decreasing = TRUE)

# run gseKEGG
 
KEGG <- gseKEGG(geneList, organism = "mmu", keyType = "kegg", exponent = 1, nPerm = 10000, minGSSize = 10, maxGSSize = 500, pvalueCutoff = 0.25, pAdjustMethod = "BH", verbose = TRUE, use_internal_data = FALSE, seed = 1, by = "fgsea")









