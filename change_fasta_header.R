# This script reads in a fasta file and an input string and renames the header to the new string. Will only work with one sequence fastas.
library(seqinr)

### Functions:
# Write fasta file from seqinr fasta object. Inputs are fasta object and string for file name.
write_fasta_seqinr <- function(fasta_seqs, file_name){
  # Get all sequences on a vector.
  seqs <- lapply(1:length(fasta_seqs), function(x) fasta_seqs[[x]])
  write.fasta(sequences = seqs, names = names(fasta_seqs), file.out = (file_name))
  # Print location.
  print(paste0("Fasta write to: ", getwd()))
}

### Read inputs:
# Allow argument usage.
args = commandArgs(trailingOnly = TRUE)
# Print required input file if typed help.
if (args[1] == "-h" || args[1] == "help"){
  print("Syntax: Rscript.R file.fata string")
  q()
  N
}
input_file = args[1]
input_file2 = args[2]

# Read the fasta file and change its name to the input string.
fasta_file <- read.fasta(input_file, as.string = TRUE, forceDNAtolower = TRUE, set.attributes = FALSE)
names(fasta_file)[1] <- input_file2

print(paste0(input_file, "_cns.fasta header renamed to ", input_file2))
# Write out new fasta with same name to overwrite the existing one.
write_fasta_seqinr(fasta_file, paste0(input_file2, "_cns.fasta"))
