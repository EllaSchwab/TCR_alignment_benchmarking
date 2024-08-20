#Author: Ella Schwab
#This script will generate 5 "healthy donor" repertories and generate 5 fasta files each containing
#5,000 sequences which reads will be generated from 

#install the immuneSIM package
install.packages("immuneSIM")

#load the package
library(immuneSIM)

#install and load readr
install.packages("readr")
library(readr)

install.packages("dplyr")
library(dplyr)

# Read the text files
reference <- read_lines("IMGT_ref.txt")
C_region_reference <- read_lines("TRBC1_ref.txt")


# Initialize variables
trbc_sequence <- ""
recording <- FALSE

# Loop through each line of the C region reference file
for (line in C_region_reference) {
  if (startsWith(line, ">M12887|TRBC1")) {
    # Start recording after the header line
    recording <- TRUE
  } else if (recording && startsWith(line, ">")) {
    # Stop recording if another header line is encountered
    break
  } else if (recording) {
    # Append the line to the sequence
    trbc_sequence <- paste0(trbc_sequence, line)
  }
}

#Grab V alleles from reference file 
# Initialize vectors to store the extracted data
genes <- c()
alleles <- c()
sequences <- c()

# Process the file line by line
for (i in seq_along(reference)) {
  line <- reference[i]
  
  if (startsWith(line, ">")) {
    # Split the header line by '|'
    parts <- strsplit(sub("^>", "", line), "\\|")[[1]]
    
    # Extract the gene and allele information
      gene_allele <- parts[2]
      gene_parts <- strsplit(gene_allele, "\\*")[[1]]
      gene <- gene_parts[1]
      allele <- gene_parts[2]
      
      # Extract the sequence from subsequent lines until the next ">"
      sequence <- ""
      for (j in (i + 1):length(reference)) {
        if (startsWith(reference[j], ">")) break
        sequence <- paste0(sequence, reference[j])
      }
      
      # Append to vectors
      genes <- c(genes, gene)
      alleles <- c(alleles, allele)
      sequences <- c(sequences, sequence)
    }
  
}

# Create the dataframe
NewGeneList <- data.frame(
  gene = genes,
  allele = alleles,
  sequence = sequences,
  stringsAsFactors = FALSE
)

# Merge the two data frames
NewGeneList <- NewGeneList %>%
  left_join(list_germline_genes_allele_01$hs$tr$b$V %>% select(gene, allele, species, frequency), by = c("gene", "allele"))

# Process the data to drop the genes that aren't in Emerson 2017
NewGeneList <- NewGeneList %>%
  group_by(gene) %>% 
  # Filter out groups where all frequencies are NA
  filter(any(!is.na(frequency))) %>%
  # Replace all frequencies in a group with the first non-NA frequency
  mutate(frequency = first(frequency[!is.na(frequency)]),
         species = first(species[!is.na(species)])) %>%
  ungroup()  # Ungroup to remove the grouping structure


#just make a new list by adding onto the reference dataset in ImmuneSIM for TRB samples 
final_gene_list <-list_germline_genes_allele_01

#change the default V gene list to be our new list of V genes and alleles
final_gene_list$hs$tr$b$V <- NewGeneList

#removing the empty D gene from the list of possible D genes:
final_gene_list$hs$tr$b$D <-final_gene_list$hs$tr$b$D[-1, ]

# For V, D, and J alleles - combine 'gene' and 'allele' into 'gene_name'
# List of segments to process
segments <- c("V", "D", "J")

# Loop through each segment and combine 'gene' and 'allele'
for (segment in segments) {
  final_gene_list$hs$tr$b[[segment]]$gene <- paste(
    final_gene_list$hs$tr$b[[segment]]$gene, 
    final_gene_list$hs$tr$b[[segment]]$allele, 
    sep = "_"
  )
}

#simulate a repertoire for HD sample 
#equal cc will ensure that the clonal abundance is equal across all clones

# Set seed for reproducibility
set.seed(1234)

# Define the number of sequences and number of samples
num_seqs <- 5000
num_samples <- 5


# Define the output directory
output_dir <- "."

# Initialize a list to store the data frames
hd_samples <- list()


# Loop to simulate each HD sample and save as CSV
for (i in 1:num_samples) {
  # Simulate the repertoire for the current HD sample
  hd_samples[[i]] <- immuneSIM(
    number_of_seqs = num_seqs,
    vdj_list = final_gene_list,
    species = "hs",
    receptor = "tr",
    chain = "b",
    max_cdr3_length = 20,
    min_cdr3_length = 6,
    equal_cc = TRUE,
    verbose = TRUE,
    airr_compliant = TRUE
  )
  
  # Append TRBC sequence to the sequences in the simulated data frame
  hd_samples[[i]]$sequence <- paste0(hd_samples[[i]]$sequence, trbc_sequence)
  
  # Save the current HD sample as a CSV file
  output_path <- file.path(output_dir, paste0("HD", i, ".csv"))
  write.csv(hd_samples[[i]], file = output_path, row.names = TRUE)
}


#create a .fasta file for each sample/repertoire so that reads can be generated  

charVectorToFasta <- function(sequences, outputName) {
  fastaContent <- vector(mode = "character", length = length(sequences) * 2)
  
  for (i in 1:length(sequences)) {
    fastaContent[2*i - 1] <- as.character(paste(">seq", i, sep = ""))
    fastaContent[2*i] <- as.character(sequences[i])
  }
  
  outFile <- file(outputName)
  writeLines(fastaContent, outFile)
  close(outFile)
}

# List of data frames and corresponding file names
data_frames <- lapply(hd_samples, function(df) df$sequence)
output_files <- c("HD1.fasta", "HD2.fasta", "HD3.fasta", "HD4.fasta", "HD5.fasta")

# Directory to save FASTA files
output_dir <- "./fasta"

# Generate FASTA files for each data frame
for (i in seq_along(data_frames)) {
  # Define the full file path
  output_path <- file.path(output_dir, output_files[i])
  
  # Save the FASTA file
  charVectorToFasta(data_frames[[i]], output_path)
  
  # Print the first few lines of the generated FASTA file for verification
  fastaContent <- readLines(output_path)
  cat("First few lines of", output_files[i], ":\n")
  print(head(fastaContent))
  cat("\n")
}


#blah blah 

