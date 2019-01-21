# Written by SJK (17 Jan 2019) 
# This file is to create rRNA mapping QC plot
# After rRNA STAR mapping, I want to look each sample rRNA mapping rate
# In here, reshape2(1.4.3) dplyr(0.7.8) and ggplot2(3.1.0)
# For usage, the number of samples has to be known in advance

#########################################################################

# rRNA mapping QC data load
r_RNA_mapping <- read.table("path/to/file")

########################################
# Function for data manipulation rRNA QC
########################################

rRNA_mapping_QC_DF <- function(n_samples, DF){  # two arguments: the number of samples and Data Frame
  
  DF$V1 <- as.character(DF$V1) # changing DF from factor to character
  DF$V2 <- as.character(DF$V2)
  
  #######################################
  # reshaping matrix
  output <- matrix(ncol=n_samples, nrow=4) # empty matrix
  j = 1
  for(i in c(1:n_samples)){ # for loop for reshaping matrix from (n x 2) to (4 x n_samples)  
    output[1, i] <- DF[j,2]
    output[2, i] <- DF[j+1,2]
    output[3, i] <- DF[j+2,2]
    output[4, i] <- DF[j+3,2]
    
    j = j + 4
  }
  output <- data.frame(output) # matrix to data frame 
  
  ########################################  
  # This is function for colname of DF
  output_name <- matrix(ncol=1, nrow=n_samples) # empty matrix
  j = 1
  for(i in c(1:n_samples)){ # extracting only sample ID 
    output_name[i, 1] <- DF[j,2]
    j = j + 4
  }
  
  output_name_list <- data.frame(output_name)
  output_name_list$output_name <- as.character(output_name_list$output_name)
  rownames(output) <- c("Sample","Number_of_input_reads","Uniquely_mapped_reads(%)","Multi_mapped_reads(%)")
  colnames(output) <- output_name_list$output_name 
  
  #########################################
  # Should modify by experiment (This is for KCl data)
  # output <- output %>% dplyr::mutate_all(as.character) # change all columns to character
  # output[1, ] <- substring(output[1, ], 3, 10) # trimming sample name
  # output[1, ] <- sub("*_|*_S","", output[1, ]) # remove _ or _S in the 1 row
  # colnames(output) <- output[1, ] # remove first row
  # output <- output[-c(1), ] # remove first row which has sample name
  
  # Reshaping for ggplot2
  output$values <- rownames(output)
  output <- melt(output, id.vars = "values")
  output$variable <- gsub("_", "", output$variable)
  output$value <- gsub("_", "", output$value)
  output$value <- as.numeric(output$value)
  return(output)
}


##########################
# Creating plot function
##########################

Plot_rRNA_mapped <- function(DF, exp_list){ # Needs 2 arguments Data Frame & exp_list
  
  num1 <- seq(from = 1, to = nrow(DF), by = 1) # All number until number of exp * 3
  num2 <- seq(from = 1, to = nrow(DF), by = 3) # only 1st number each exp
  num3 <- num1[!num1 %in% num2] # In here only remains 2,3, 5,6, 8,9...
  DF <- DF[num3, ] # Extract every 2nd and 3rd rows each sample
  DF <- DF %>% group_by(variable) %>% summarise(value = sum(value)) # Add uniquely mapped + multi mapped
  DF$variable <- factor(DF$variable, exp_list) # To look plot in order
  
  return(ggplot(data=DF, aes(x=variable, y=value)) + geom_bar(stat="identity", fill = "#990033") + ylab("rRNA mappability (%)") + scale_y_continuous(expand = c(0,0), limits = c(0, 100))) 
}


# e.g. KCl Data creating plot

rRNA_Matrix <- rRNA_mapping_QC_DF(18, r_RNA_mapping) # Create matrix

exp_list <- c("HET1ctr", "HET1KCl", "WT1ctr", "WT2ctr", "WT4ctr", "WT5ctr", "WT1KCl", "WT2KCl", "WT4KCl", "WT5KCl", "KO2ctr", "KO3ctr", "KO4ctr", "KO5ctr", "KO2KCl", "KO3KCl", "KO4KCl", "KO5KCl")
plot1 <- Plot_rRNA_mapped(rRNA_M, exp_list) # Create plot

plot1 <- plot1 + theme(axis.text.y = element_text(size = 12, face = "bold", color = "black"), 
                   axis.text.x = element_text(angle = 90, hjust = .8, size = 10, vjust = 0.4, face = "bold"), 
                   axis.title.y = element_text(size = 12, face = "bold"), 
                   axis.title.x = element_blank(), 
                   legend.position = c(0.95, 0.7), legend.title = element_blank(), legend.text = element_text(face = "bold", size = 12),
                   panel.background = element_rect(fil = "white"), axis.line = element_line()) + 
  geom_vline(xintercept = 2.5, colour = "#CC0033", linetype="dashed", size= 0.2) +
  geom_vline(xintercept = 6.5, colour = "#CC0033", linetype="dashed", size= 0.2) +
  geom_vline(xintercept = 10.5, colour = "#CC0033", linetype="dashed", size= 0.2) + 
  geom_vline(xintercept = 14.5, colour = "#CC0033", linetype="dashed", size= 0.2) +
  annotate("text", x = 4.5, y = 95, label = "WT(ctr)", colour = "#003366") + 
  annotate("text", x = 8.5, y = 95, label = "WT(KCl)", colour = "#003366") +
  annotate("text", x = 12.5, y = 95, label = "KO(ctr)", colour = "#CC3300") + 
  annotate("text", x = 16.5, y = 95, label = "KO(KCl)", colour = "#CC3300")



