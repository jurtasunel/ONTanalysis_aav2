### This script produces a coverage plot pdf and a tsv with low depth positions.
### It requires two files as argument inputs:
### 1: consensus_genomes.fasta with all the medaka consensus.
### 2: quality_report.csv with headers "prefix_barcode", "average_depth", "perc_over_30X".

### Libraries:
library(seqinr)
library(stringr)
library(ggplot2)

### Read inputs:
# Allow argument usage.
args = commandArgs(trailingOnly = TRUE)

# Get input file with depth from command line and print it.
input_file = args[1] # consensus_genomes.fasta.
input_file2 = args[2] # quality report csv.
print("Reading input files...")
# Read in the quality csv and add extra column for the coverage.
quality_report <- read.csv(input_file2, header = TRUE, stringsAsFactors = FALSE)
quality_report <- cbind(quality_report, coverage = NA)
# Read in the fasta with all the sequences and rename names to remove "Artic/Medaka" tag.
fasta_file <- read.fasta(input_file, as.string = TRUE, forceDNAtolower = TRUE, set.attributes = FALSE)

### Quality filter:
print("Applying quality filter...")
# Make vectors to store the coverage of each sequence, and coverage and ID of low coverage sequences.
nt_perc <- c()
low_quality_seqs_ID <- c()
low_quality_seqs_Nperc <- c()
# Loop through the fasta file.
for (i in 1:length(fasta_file)){
  
  # Save the current ID.
  current_ID <- names(fasta_file)[i]
  print(current_ID)
  # Get the percentage of N of each sequence and the percente of nt.
  N_perc <- (str_count(fasta_file[[i]], pattern = "n") * 100) / nchar(fasta_file[[i]])
  coverage <- (100 - N_perc)
  nt_perc <- c(nt_perc, coverage)
  print(nt_perc)
  # Add the current coverage to the coverage column on the quality report row matching the current ID. 
  quality_report[quality_report$prefix_barcode == current_ID, "coverage"] <- coverage
  print(paste0("Coverage of ", current_ID, " = ", coverage))
  
  # If N percentage is higher than 20, get the value and the Id to separate vectors.
  if (N_perc > 10){
    low_quality_seqs_ID <- c(low_quality_seqs_ID, names(fasta_file)[i])
    low_quality_seqs_Nperc <- c(low_quality_seqs_Nperc, N_perc)
  }
}

# Overwrite the existing quality report with the new one. 
print("Appending quality report with samples coverage...")
write.csv(quality_report, "quality_report.csv", row.names = FALSE)
quality_report

# Make a dataframe to plot the coverage of all IDs
df_plot_quality <- data.frame(names(fasta_file), nt_perc)
colnames(df_plot_quality) <- c("ID", "Coverage")

print("Plotting coverage...")
# Plot barplot
quality_plot <- ggplot(df_plot_quality, aes(x = ID, y = Coverage)) +
  geom_bar(stat = "identity", fill = "lightblue") +
  geom_hline(yintercept = 90, color = "red") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.title.y = element_text(),
        axis.text.y = element_text(size = 12),
        axis.ticks.y = element_line(),
        axis.title.x = element_text(),
        axis.text.x = element_text(angle = 90, size = 4, )) +
  ggtitle("Quality filter: Coverage > 90%") +
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold")) +
  labs(caption = paste0("Sequences excluded from analysis: ", length(low_quality_seqs_Nperc)))

ggsave("Coverage_barplot", quality_plot, device = "png", dpi = 700, bg = "white")

# Remove sequences of low quality from fasta file and metadata file.
#fasta_file <- fasta_file[! names(fasta_file) %in% low_quality_seqs_ID]







