# Written by SJK (17 Jan 2019) 
# This file is to create mapping QC plot
# After STAR mapping, I want to look each sample input reads and uniquely mapping rate
# In here, reshape2(1.4.3) and ggplot2(3.1.0)
# For usage, the number of samples has to be known in advance

#################################################

mapping <- read.table("path/to/file")

########################################
# Function for data manipulation QC
########################################

mapping_QC_DF <- function(n_samples, DF){  # two arguments: the number of samples and Data Frame
  
  DF$V1 <- as.character(DF$V1) # changing DF from factor to character
  DF$V2 <- as.character(DF$V2)
  
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
  rownames(output) <- c("Sample","Number_of_input_reads","Number_of_mapped_reads","Uniquely_mapped_reads(%)")
  colnames(output) <- output_name_list$output_name 
  
  #########################################
  # Reshaping for ggplot2
  
  output$values <- rownames(output)
  output <- melt(output, id.vars = "values")
  output$variable <- gsub("_", "", output$variable)
  output$value <- gsub("_", "", output$value)
  
  return(output)
}

########################################
# Function creating input read plot
########################################

Plot_Input_Reads <- function(DF, exp_list){ # Needs 1 argument as Data Frame

  num <- seq(from = 2, to  = nrow(DF), by = 4) 
  DF <- DF[num, ] # Extract every 2nd row out of 4 rows
  colnames(DF) <- c("info","exp", "values") 
  DF$values <- as.integer(DF$values) # change from character to integer
  DF$exp <- factor(DF$exp, exp_list)
  return(ggplot(data=DF, aes(x=exp, y=values)) + geom_bar(stat="identity", fill = "#0066CC") + ylab("# of Input reads")) 
}

########################################
# Function creating mapped reads plot
########################################

Plot_Mapped_Reads <- function(DF, exp_list){
  num <- seq(from = 3, to  = nrow(DF), by = 4) 
  DF <- DF[num, ] # Extract every 3rd row out of 4 rows
  colnames(DF) <- c("info","exp", "values") 
  DF$values <- as.integer(DF$values) # change from character to integer
  DF$exp <- factor(DF$exp, exp_list)
  return(ggplot(data=DF, aes(x=exp, y=values)) + geom_bar(stat="identity", fill = "#336633") + ylab("# of mapped reads")) 
}

########################################
# Function creating mappability plot
########################################

Plot_Mappability <- function(DF, exp_list){
  num <- seq(from = 4, to  = nrow(DF), by = 4) 
  DF <- DF[num, ] # Extract every 4th row out of 4 rows
  colnames(DF) <- c("info","exp", "values") 
  DF$values <- as.integer(DF$values) # change from character to integer
  DF$exp <- factor(DF$exp, exp_list)
  return(ggplot(data=DF, aes(x=exp, y=values)) + geom_bar(stat="identity", fill = "#CC6600") + ylab("Mappability (%)")+ scale_y_continuous(expand = c(0,0), limits = c(0, 100))) 
}


# e.g. ploting example of KCl data
map_Matrix <- mapping_QC_DF(18, mapping)
exp_list <- c("HET1ctr", "HET1KCl", "WT1ctr", "WT2ctr", "WT4ctr", "WT5ctr", "WT1KCl", "WT2KCl", "WT4KCl", "WT5KCl", "KO2ctr", "KO3ctr", "KO4ctr", "KO5ctr", "KO2KCl", "KO3KCl", "KO4KCl", "KO5KCl")
plot1 <- Plot_Input_Reads(map_Matrix, exp_list)

plot1  + scale_y_continuous(expand = c(0,0), limits = c(0,4e+07)) + 
  theme(axis.text = element_text(size = 12, face = "bold", color = "black"), axis.text.x = element_text(angle = 90, hjust = .8, size = 10, vjust = 0.4), 
        axis.title.y = element_text(size = 12, face = "bold"), axis.title.x = element_blank(), 
        legend.position = c(0.95, 0.7), legend.title = element_blank(), legend.text = element_text(face = "bold", size = 12),
        panel.background = element_rect(fil = "white"), axis.line = element_line()) + 
  geom_vline(xintercept = 2.5, colour = "#CC0033", linetype="dashed", size= 0.2) +
  geom_vline(xintercept = 6.5, colour = "#CC0033", linetype="dashed", size= 0.2) +
  geom_vline(xintercept = 10.5, colour = "#CC0033", linetype="dashed", size= 0.2) + 
  geom_vline(xintercept = 14.5, colour = "#CC0033", linetype="dashed", size= 0.2) +
  annotate("text", x = 4.5, y = 3.5e+07, label = "WT(ctr)", colour = "#003366") + 
  annotate("text", x = 8.5, y = 3.5e+07, label = "WT(KCl)", colour = "#003366") +
  annotate("text", x = 12.5, y = 3.5e+07, label = "KO(ctr)", colour = "#CC3300") + 
  annotate("text", x = 16.5, y = 3.5e+07, label = "KO(KCl)", colour = "#CC3300")


