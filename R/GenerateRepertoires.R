#Author: Ella Schwab
#This script will generate 5 "healthy donor" repertories and generate 5 fasta files each containing
#250,000 sequences which reads will be generated from 

#install the immuneSIM package
install.packages("immuneSIM")

#load the package
library(immuneSIM)

#install and load readr
install.packages("readr")
library(readr)

install.packages("dplyr")
library(dplyr)

# Read the text file
reference <- read_lines("IMGT_ref.txt")

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

# Define the number of sequences and number of samples
num_seqs <- 5000
num_samples <- 5

# Define the output directory
output_dir <- "/project/mangul_341/eschwab/TCR_alignment/"

# Initialize a list to store the data frames
hd_samples <- list()

# Loop to simulate each HD sample and save as CSV
for (i in 1:num_samples) {
  # Simulate the repertoire for the current HD sample
  hd_samples[[i]] <- immuneSIM(
    number_of_seqs = num_seqs,
    vdj_list = new_gene_list,
    species = "hs",
    receptor = "tr",
    chain = "b",
    equal_cc = TRUE,
    verbose = TRUE,
    airr_compliant = TRUE
  )
  
  # Save the current HD sample as a CSV file
  output_path <- file.path(output_dir, paste0("HD", i, ".csv"))
  write.csv(hd_samples[[i]], file = output_path, row.names = TRUE)
}


# Install viridis if not already installed
install.packages("viridis")

# Load viridis
library(viridis)
# Load the libraries
library(ggplot2)



# Calculate frequencies for V alleles
frequencyV_HD1 <- HD1 %>%
  count(v_call) %>%             # Count occurrences of each unique value
  rename(frequency = n)         # Rename the count column for clarity

# Sort the data frame by 'v_call'
frequencyV_HD1 <- frequencyV_HD1 %>%
  arrange(v_call)  # Change to arrange(desc(v_call)) for descending order

frequencyV_HD1 <- frequencyV_HD1 %>%
  mutate(percentage = (frequency / sum(frequency)) * 100)  # Convert to percentage


# Calculate frequencies
frequencyD_HD1 <- HD1 %>%
  count(d_call) %>%             # Count occurrences of each unique value
  rename(frequency = n)         # Rename the count column for clarity

# Sort the data frame by 'd_call'
frequencyD_HD1 <- frequencyD_HD1 %>%
  arrange(d_call)  # Change to arrange(desc(v_call)) for descending order

frequencyD_HD1 <- frequencyD_HD1 %>%
  mutate(percentage = (frequency / sum(frequency)) * 100)  # Convert to percentage

# Calculate frequencies
frequencyJ_HD1 <- HD1 %>%
  count(j_call) %>%             # Count occurrences of each unique value
  rename(frequency = n)         # Rename the count column for clarity

# Sort the data frame by 'd_call'
frequencyJ_HD1 <- frequencyJ_HD1 %>%
  arrange(j_call)  # Change to arrange(desc(v_call)) for descending order

frequencyJ_HD1 <- frequencyJ_HD1 %>%
  mutate(percentage = (frequency / sum(frequency)) * 100)  # Convert to percentage

# Create the bar plot for V alleles
plotV <- ggplot(frequencyV_HD1, aes(x = v_call, y = percentage, fill = v_call)) +
  geom_bar(stat = "identity") +    # Use 'identity' to plot frequencies
  theme_minimal() +      # Use a minimal theme
  scale_fill_viridis_d(option = "mako") +  
  labs(x = "V allele", y = "Frequency (%)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for readability

# Create the bar plot for J alleles
plotJ <-ggplot(frequencyJ_HD1, aes(x = j_call, y = percentage, fill = j_call)) +
  geom_bar(stat = "identity") +    # Use 'identity' to plot frequencies
  theme_minimal() +      # Use a minimal theme
  scale_fill_viridis_d(option = "mako") +  
  labs(x = "J allele", y = "Frequency (%)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for readability


#plot the V call and J call side by side
# Install patchwork if not already installed
install.packages("patchwork")

# Load patchwork
library(patchwork)

# Combine the plots side by side
stacked_plot <- plotV/plotJ

# Save the combined plot to a file
ggsave(filename = "stacked_plot.png", plot = stacked_plot, width = 12, height = 8, dpi = 300)


#Plot CDR3 length distribution for HD1

# Calculate CDR3 lengths and their frequencies
HD1_CDR3freq <- HD1 %>%
  mutate(CDR3_length = nchar(junction_aa)) %>%  # Calculate CDR3 lengths
  group_by(CDR3_length) %>%                     # Group by length
  summarize(count = n()) %>%                    # Count occurrences
  mutate(frequency = count / sum(count) * 100)  # Calculate frequency as a percentage

# Calculate TCR sequence lengths and their frequencies
HD1_TCRfreq <- HD1 %>%
  mutate(TCR_length = nchar(sequence)) %>%  # Calculate CDR3 lengths
  group_by(TCR_length) %>%                     # Group by length
  summarize(count = n()) %>%                    # Count occurrences
  mutate(frequency = count / sum(count) * 100)  # Calculate frequency as a percentage



# Create the bar plot with frequency on the y-axis
CDR3_lengths<- ggplot(HD1_CDR3freq, aes(x = CDR3_length, y = frequency)) +
  geom_bar(stat = "identity", fill = "skyblue", color = "black") +
  theme_minimal() +
  labs(x = "CDR3 Length (Amino Acids)", y = "Frequency (%)", 
       title = "Distribution of CDR3 Lengths")

# Create the bar plot with frequency on the y-axis and sequence length on the x
ggplot(HD1_TCRfreq, aes(x = TCR_length, y = frequency)) +
  geom_bar(stat = "identity", fill = "red", color = "black") +
  theme_minimal() +
  labs(x = "TCR Length (Nucleotides)", y = "Frequency (%)", 
       title = "Distribution of TCR Lengths")


# Save the combined plot to a file
ggsave(filename = "CDR3_lengths.png", plot = CDR3_lengths, width = 12, height = 8, dpi = 300)
