### This script produces a depth plot pdf, a csv with the average depth and a tsv with low depth positions.
### It requires two files as argument inputs:
### 1: _depth.txt file output from samtools depth.
### 2: quality_report.csv with headers "prefix_barcode", "average_depth", "perc_over_30X".

### Libraries:
library(ggplot2)

### Read inputs:
# Allow argument usage.
args = commandArgs(trailingOnly = TRUE)
# Print required input file if typed help.
if (args[1] == "-h" || args[1] == "help"){
  print("Syntax: Rscript.R depth_file.txt quality_report.csv")
  q()
  N
}

# Get input file with depth from command line and print it.
input_file = args[1] # _depth.txt file from samtools depth.
input_file2 = args[2] # quality report csv.

# Separate the file name from the path and remove file tag.
ID = tail(unlist(strsplit(input_file, "/")), n = 1)
ID = gsub("_depth.txt", "", ID)
print(paste0("Analyzing ", ID, " depth..."))
# Read the depth file as table and change column names.
depth_table <- read.table(input_file,  sep = '\t', header = FALSE)
colnames(depth_table) <- c("Chrom", "Position", "Depth")
# Read in the quality report csv.
quality_report.csv <- read.csv(input_file2, header = TRUE, stringsAsFactors = FALSE)

# Calculate the average depth and the percentage of bases with depth over 30X and print it.
average_depth <- mean(depth_table$Depth)
perc_over30 <- ((nrow(depth_table[depth_table$Depth > 30,])) * 100) / nrow(depth_table)
print(paste0("Average depth = ", average_depth, "; percentage of bases over 30X = ", perc_over30))
# Add the ID and values to the next row of the quality report, then overwrite the existing report with the new one.
quality_report.csv[nrow(quality_report.csv) + 1,] <- c(ID, average_depth, perc_over30)
write.csv(quality_report.csv, "quality_report.csv", row.names = FALSE)

print(paste0("Finding ", ID, " low depth positions..."))
# Fill a vector with a depth score for each position, and save bad score position.
Depth_score <- c()
BC_positions <- c()
for (i in 1:nrow(depth_table)){
  
  # 30X is routinary cutoff for ONT sequencing depth.
  if (depth_table$Depth[i] < 30){
    Depth_score <- c(Depth_score, "C (<30)")
    BC_positions <- c(BC_positions, depth_table$Position[i])
    
  } else if (depth_table$Depth[i] < 60){
    Depth_score <- c(Depth_score, "B (30-60)")
    BC_positions <- c(BC_positions, depth_table$Position[i])
    
  } else {Depth_score <- c(Depth_score, "A (>60)")}
}

# Add a column to the depth table with the depth scores vector.
depth_table <- cbind(depth_table, Depth_score)

# Get ranges of values of B an C positions.
BC_range <- split(BC_positions, cumsum(c(1, diff(BC_positions) != 1)))
# Get first and last value of each range.
Low_depth_positions <- c()
for (i in 1:length(BC_range)){
  current_range <- as.character(unlist(BC_range[i]))
  if (length(current_range) == 1){
    Low_depth_positions <- c(Low_depth_positions, current_range)
  } else{
    Low_depth_positions <- c(Low_depth_positions, paste0(head(current_range, 1), "-", tail(current_range, 1)))
  }
}
# Write the ranges of values for low depth positions to a csv file.
Low_depth_positions <- data.frame(Low_depth_positions)
write.csv(Low_depth_positions, file = paste0(ID,"_low_depth_positions.csv"), row.names = FALSE)

print(paste0("Plotting ", ID, " depth..."))
# Plot depth table.
p <- ggplot(depth_table, aes(x = Position, y = Depth, fill = Depth_score)) +
  geom_bar(stat = "identity") + # Fill identity for ggplot to accept the y axis data on geom bar.
  theme_minimal() + # Remove background grid.
  scale_fill_manual(values = c("#339966", "#3333FF", "#CC0033")) +
  theme(legend.position = "top") +
  ggtitle(ID) +
  theme(plot.title = element_text(hjust = 0.5, size = 14))

# Create pdf to save the plot.
pdf(paste0(ID,"_depth.pdf"))
print(p) # Save plot on first page
#print(p2)
#print(p3).. for saving on next pages of pdf.
  
