
install.packages('Seurat')
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("infercnv")
BiocManager::install("biomaRt")

#load libraries
library(Seurat)
library(dplyr)

#clear existing environment
rm(list = ls())

#set the directory
setwd("C:/Users/email/Isha/columbia/azizi_lab/project_with_joy/actual code/wgs and scrna-seq mapping")

#read the infercnv object (contains list of genes)
infercnv_obj <- readRDS("ribas_310_pre_run.final.infercnv_obj")
given_genes <- infercnv_obj@gene_order

#load the table with genomic regions from the seg file
genomic_regions <- read.table(file = "Ribas-1Pre_S102_L001_tumor.seg")
names(genomic_regions) <- as.matrix(genomic_regions[1, ])
genomic_regions <- genomic_regions[-1, ]

sno <- 1:28
genomic_regions <- cbind(genomic_regions, sno)

#add column to store the identified genes for each genomic region
identified_gene <- rep('', times = 28)
genomic_regions <- cbind(genomic_regions, identified_gene)

#add column to store the identified genomic region for each gene
identified_region <- rep(0, times = 7661)
given_genes <- cbind(given_genes, identified_region)


for (k in 1:7661){
  
  #store the start coordinate, end coordinate, and the chromosome number
  chr_str <- given_genes$chr[k]
  start_gene <- strtoi(given_genes$start[k])
  end_gene <- strtoi(given_genes$stop[k])

  #filter the genomic_regions based on chromosome number
  chr_num <- substring(chr_str, 4)
  which_array <- which(genomic_regions$chr == chr_num)
  end <- which_array[1] + length(which_array) - 1
  genomic_lst_chr <- genomic_regions[which_array[1]:end,]
  
  num_regions <- length(genomic_lst_chr$start)
  lengths <- rep(0, times = num_regions)
  
  for (l in 1:(num_regions)){
    start_genomic <- strtoi(genomic_lst_chr$start[l])
    end_genomic <- strtoi(genomic_lst_chr$end[l])
    
    if (start_gene >= start_genomic & end_gene <= end_genomic){
      given_genes$identified_region[k] = genomic_lst_chr$sno[l]
    } 
    else if (start_gene >= start_genomic & end_gene > end_genomic){
      len <- end_genomic - start_gene
      
    }
    else if (start_gene < start_genomic & end_gene <= end_genomic){
      len <- end_gene - start_genomic
    } 
    
    lengths[l] <- len
    
  }
  
  max_len <- max(lengths)
  max_pos <- which.max(lengths)
  
  if (max_len > 0){
    given_genes$identified_region[k] = genomic_lst_chr$sno[max_pos]
  }
  
} 


for (i in 1:27){
  
  which_array <- which(given_genes$identified_region == strtoi(i))
  end <- which_array[1] + length(which_array) - 1
  gene_lst_chr <- given_genes[which_array[1]:end,]
  
  #obtain the names of the genes corresponding to each genomic region
  genes <- row.names(gene_lst_chr)
  
  #add this list to the genomic_regions df
  genomic_regions$identified_gene[i] <- paste(unlist(genes), sep='', collapse=', ')
  
}

#write genomic_regions df to a csv
write.csv(genomic_regions,"Ribas-1Pre_S102_L001_tumor_genes-identified-multiple-diffnew.csv", row.names = FALSE)
