#Author: Ella Schwab
#This script will plot repertoire statistics for each of the 
#5 samples as a sanity check before generating reads

# Install viridis if not already installed
install.packages("viridis")

# Load necessary libraries
library(dplyr)
library(ggplot2)
library(patchwork)  # For combining plots
library(viridis)    # For color scales

# Set the directory where your .csv files are located
directory <- "."

# Get the names of the .csv files in the directory
csv_files <- list.files(path = directory, pattern = "HD\\d+\\.csv", full.names = TRUE)

# Read each .csv file and assign it to a data frame with the appropriate name
for (file in csv_files) {
  # Extract the base name without extension
  file_name <- tools::file_path_sans_ext(basename(file))
  
  # Read the .csv file into a data frame
  df <- read.csv(file)
  
  # Create a variable name based on the file name
  # Assign the data frame to a variable in the global environment
  assign(file_name, df, envir = .GlobalEnv)
}



# List of data frame names
data_frames <- list(HD1, HD2, HD3, HD4, HD5)

# Create an empty list to store results
frequency_list <- list()

# Calculate frequencies for each data frame
for (i in seq_along(data_frames)) {
  df <- data_frames[[i]]
  
  frequency_df <- df %>%
    count(v_call) %>%             # Count occurrences of each unique value
    rename(frequency = n)         # Rename the count column for clarity
  
  frequency_df <- frequency_df %>%
    arrange(v_call) %>%           # Sort the data frame by 'v_call'
    mutate(percentage = (frequency / sum(frequency)) * 100)  # Convert to percentage
  
  # Store the result in the list with a name like "frequencyV_HD1", "frequencyV_HD2", etc.
  frequency_list[[paste0("frequencyV_HD", i)]] <- frequency_df
}


# Access individual data frames
frequencyV_HD1 <- frequency_list[["frequencyV_HD1"]]
frequencyV_HD2 <- frequency_list[["frequencyV_HD2"]]
frequencyV_HD3 <- frequency_list[["frequencyV_HD3"]]
frequencyV_HD4 <- frequency_list[["frequencyV_HD4"]]
frequencyV_HD5 <- frequency_list[["frequencyV_HD5"]]

create_bar_plot <- function(df, sample_name) {
  ggplot(df, aes(x = v_call, y = percentage, fill = v_call)) +
    geom_bar(stat = "identity") +    # Use 'identity' to plot frequencies
    theme_minimal() +      # Use a minimal theme
    scale_fill_viridis_d(option = "mako") +  
    labs(x = "V allele", y = "Frequency (%)", title = paste("V Alleles -", sample_name)) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),  # Adjust angle and size
      axis.title.x = element_text(size = 12),  # Adjust x-axis title size
      axis.title.y = element_text(size = 12),  # Adjust y-axis title size
      legend.position = "right",  # Position legend to the right
      legend.title = element_text(size = 10),  # Reduce legend title size
      legend.text = element_text(size = 8),    # Reduce legend text size
      legend.box.margin = margin(0, 0, 0, 0)  # Remove margin around legend box
    )
}

# Generate plots for each data frame
plot_list <- list (
plotV_HD1 = create_bar_plot(frequencyV_HD1, "HD1"),
plotV_HD2 = create_bar_plot(frequencyV_HD2, "HD2"),
plotV_HD3 = create_bar_plot(frequencyV_HD3, "HD3"),
plotV_HD4 = create_bar_plot(frequencyV_HD4, "HD4"),
plotV_HD5 = create_bar_plot(frequencyV_HD5, "HD5")

)


# Specify the directory and file name
output_dir <- "./figures"

# Loop through the plot list to save each plot
for (name in names(plot_list)) {
  # Define the file path
  file_name <- paste0(name, ".jpeg")
  output_path <- file.path(output_dir, file_name)
  
  # Save the plot
  ggsave(filename = output_path, plot = plot_list[[name]], width = 22, height = 8, dpi = 300)
}




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
