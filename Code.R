install.packages("dplyr")
install.packages("ggplot2")
install.packages("reshpae2")
install.packages("gtools") # Necessary for sorting position step
library(dplyr)
library(ggplot2)
library(reshape2)
library(purrr)
library(readr)
library(gtools)

#set working directory
setwd("/Users/supantha/Documents/R Assignment/")

#Data Inspection

#read in fang_genotypes text file and set it as genotype_fang

genotype_fang <- read.table("fang_et_al_genotypes.txt", header = T)
head(genotype_fang)
str(genotype_fang)
ncol(genotype_fang)
nrow(genotype_fang)

## Observe the SNP data

SNP <- read.delim("snp_position.txt",header=T)
head(SNP)
str(SNP)
ncol(SNP)
nrow(SNP)


## Data Processing

SNP_processed <- SNP[ ,c(1,3,4)] ##Since SNP_ID, Chromosome, and Position are in these columns

# Maize

#filtering for MAIZE Data
maize <- filter(genotype_fang, Group == "ZMMIL" | Group == "ZMMLR" | Group == "ZMMMR") 


#transpose maize data  so that it reads horizontally
maize_transposed <- t(maize)

#use merge function to combine both files together using the SNP_ID column in processed
##SNP file and the row names in  maize_transposed after doing the last step.
maize_merged <- merge(SNP_processed, maize_transposed, by.x = "SNP_ID", by.y = "row.names")

#Position in increasing order
maize_merged_ascend <-maize_merged[mixedorder(maize_merged$Position),]

#process for creating output
chromosomes <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10")

# Function to filter, arrange, and write CSV
write_ordered_csv <- function(chrom) {
  filtered_data <- filter(maize_merged_ascend, Chromosome == chrom) 
  
  
  # Write filtered data to a CSV file
  write_csv(filtered_data, path = paste0("Maize/maize_chrom", chrom, "_increase.csv"))
}

  # Apply the function to each chromosome level
  walk(chromosomes, write_ordered_csv)

####Maize Descending order
  #Position in decreasing order
  maize_merged_descend <- maize_merged[(mixedorder(maize_merged$Position,decreasing = T)),]
  
  #Convert matrix to data frame and replace (?) with (-)

  # gsub is used here to to replace "?" with "-" across all cells
  maize_merged_descend <- as.data.frame(apply(maize_merged_descend, 2, function(x) gsub("\\?", "-", x)))
  
  
  # Function to filter, arrange, and write CSV
  write_ordered_csv <- function(chrom) {
    filtered_data <- filter(maize_merged_descend, Chromosome == chrom) 
    
    
    # Write filtered data to a CSV file
    write_csv(filtered_data, path = paste0("Maize/maize_chrom", chrom, "_decrease.csv"))
  }
  
  # Apply the function to each chromosome level
  walk(chromosomes, write_ordered_csv)
  
  
###############Teosinte

teosinte <- filter(genotype_fang, Group == "ZMPBA" | Group == "ZMPIL" | Group == "ZMPJA")  

#transpose 
teosinte_transposed <- t(teosinte)

#merge
teosinte_merged <- merge(SNP_processed, teosinte_transposed, by.x = "SNP_ID", by.y = "row.names")

#####Position in increasing order
teosinte_merged_ascend <-teosinte_merged[mixedorder(teosinte_merged$Position),]


# Function to filter, arrange, and write CSV
write_ordered_csv <- function(chrom) {
  filtered_data <- filter(teosinte_merged_ascend, Chromosome == chrom) 
  
  # Write filtered data to a CSV file
  write_csv(filtered_data, path = paste0("Teosinte/teosinte_chrom", chrom, "_increase.csv"))
}

# Apply the function to each chromosome level
walk(chromosomes, write_ordered_csv)

####Teosinte Descending order
#Position in decreasing order
teosinte_merged_descend <- teosinte_merged[(mixedorder(teosinte_merged$Position,decreasing = T)),]

#Convert matrix to data frame and replace (?) with (-)

# gsub is used here to to replace "?" with "-" across all cells
teosinte_merged_descend <- as.data.frame(apply(teosinte_merged_descend, 2, function(x) gsub("\\?", "-", x)))


# Function to filter, arrange, and write CSV
write_ordered_csv <- function(chrom) {
  filtered_data <- filter(teosinte_merged_descend, Chromosome == chrom) 
  
  
  # Write filtered data to a CSV file
  write_csv(filtered_data, path = paste0("teosinte/teosinte_chrom", chrom, "_decrease.csv"))
}

# Apply the function to each chromosome level
walk(chromosomes, write_ordered_csv)



######Visualization
#Melting the original genotypes.fang file

both_long <- filter(genotype_fang, Group == "ZMMIL" | Group == "ZMMLR" | Group == "ZMMMR" | Group == "ZMPBA" | Group == "ZMPIL" | Group == "ZMPJA")
both <- melt(both_long, measure.vars = colnames(genotype_fang)[4:986])
colnames(both)[4:5] <- c("SNP_ID", "Homozygous")
colnames(both)


# Change all homozygous SNPs to TRUE
both <- mutate(both, Homozygous = ifelse(Homozygous %in% c("A/A", "C/C", "G/G", "T/T"), TRUE, Homozygous))

# Change all heterozygous SNPs to FALSE
both <- mutate(both, Homozygous = ifelse(Homozygous %in% c("A/C", "A/G", "A/T", "C/G", "C/T", "G/T"), FALSE, Homozygous))
# Change all missing values to NA
both <- mutate(both, Homozygous = ifelse(Homozygous %in% c("?/?"), NA, Homozygous))

#Sorting the dataframe

both <- arrange(both, Sample_ID, Group)

#my favorite palette

npg_palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", 
                 "#D55E00", "#CC79A7", "#999999", "#000000", "#FFFFFF",
                 "#FF5733", "#FFC300", "#DAF7A6", "#FF5733", "#C70039")


# Plot maize data with some changes and better settings
ggplot(data = maize_merged) +
  geom_bar(mapping = aes(x = Chromosome, fill = Chromosome)) +  # Fill bars with chromosome colors
  scale_x_discrete(limits = c(1:10, "unknown", "multiple")) +
  scale_fill_manual(values = npg_palette) +  # Set color palette to npg
  ggtitle(label = "SNPs per chromosome") +
  xlab(label = "Chromosome") +
  ylab(label = "Number of SNPs") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),  # Center the plot title
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12)
  )

ggsave(filename = "Total_SNPs_across_Maize.png", width = 12, height = 8, dpi = 300)

# Plot teosinte data with some changes and better settings
ggplot(data = teosinte_merged) +
  geom_bar(mapping = aes(x = Chromosome, fill = Chromosome)) +  # Fill bars with chromosome colors
  scale_x_discrete(limits = c(1:10, "unknown", "multiple")) +
  scale_fill_manual(values = npg_palette) +  # Set color palette to npg
  ggtitle(label = "SNPs per chromosome") +
  xlab(label = "Chromosome") +
  ylab(label = "Number of SNPs") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),  # Center the plot title
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12)
  )

ggsave(filename = "Total_SNPs_across_teosinte.png", width = 12, height = 8, dpi = 300)
## In both cases, the chromosome 1 has highest number of SNPs, followed by chromosome 2. It makes
   #sense since chromosome 1 is the largests in size.

#Plotting heterozygotes-homozygotes by sample and group

ggplot(data = both) +
  geom_bar(mapping = aes(x = Sample_ID, fill = Homozygous), stat = "count") +
  ggtitle(label = "SNPs by Ordered Sample_ID") +
  ylab(label = "Number of SNPs") +
  scale_fill_manual(values = npg_palette) +  # Set color palette to npg
  ggtitle(label = "SNPs across sample") +
  xlab(label = "Sample") +
  ylab(label = "Number of SNPs") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),  # Center the plot title
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12)
  )
ggsave(filename = "SNPs by ordered samples.png", width = 12, height = 8, dpi = 300)

ggplot(data = both) +
  geom_bar(mapping = aes(x = Group, fill = Homozygous), stat = "count") +
  ggtitle(label = "SNPs by  groups") +
  ylab(label = "Number of SNPs") +
  scale_fill_manual(values = npg_palette) +  # Set color palette to npg
  ggtitle(label = "SNPs across groups") +
  xlab(label = "Group") +
  ylab(label = "Number of SNPs") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),  # Center the plot title
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12)
  )
ggsave(filename = "SNPs by ordered Groups.png", width = 12, height = 8, dpi = 300)

# Most SNPs are homozygous among samples, withe around 20 percent being heterozygotes.
#A good number is NA too. Among groups, ZMMLR has the highest number of SNPs, with most them 
#being homozygous, same for all other groups too.

#My own analysis
#Here I plotted only ZMMLR group from the group sections. The purpose was to observe what happens in this group specifically. 

ggplot(data = both) +
  geom_bar(mapping = aes(x = 'ZMMLR', fill = Homozygous), stat = "count") +
  ggtitle(label = "SNPs of ZMMLR") +
  ylab(label = "Number of SNPs") +
  scale_fill_manual(values = npg_palette) +  # Set color palette to npg
  ggtitle(label = "SNPs across groups") +
  xlab(label = "Group") +
  ylab(label = "Number of SNPs") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),  # Center the plot title
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12)
  )
ggsave(filename = "Own Analysis_SNPs ZMMLR.png", width = 12, height = 8, dpi = 300)

#We see around 2000000 homozygous, and 500000 as heterozygous and NA.
