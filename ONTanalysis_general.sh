#!bin/bash
##################################################################################################################################################################################################
  ### This script maps ONT reads from tile-amplicons to a reference with minimap2 and produces depth plots for each sample and a quality report csv.
  ### On data path should be ONLY directories with the Sample names, and within each sample directory should be ONLY barcodes directories containing fastq.gz files
  ### The pipeline goes as follows:
  # - get all sample names
  # - Move to each sample dir and get all barcodes for that sample.
  # - Remove adapters with porechop and merge all adpter-removed fastqs into a single sample.fastq
  # - align with minimap2, make sam>bam>srt_bam>consensus.fasta
  # - call Rscript to rename consensus fasta header to be the sample name (by default is the reference genome header name).
  # - get depth and move all file to a sample_results directory.
  # When all samples have been done, create a quality_report.csv and go through each results folder to fill it and plot depth with Rscript. Also concatenate all consensus into a single fasta.
  # Call Rscript to add coverage to the quality report and plot coverage of all samples.
##################################################################################################################################################################################################

# Get paths of data, reference and scripts directory.
data_path="/path/to/data"
reference="/path/to/reference/reference.fasta"
scripts="/path/to/scripts_dir/"

# Get all the sample names on the data path.
files=(${data_path}/*) # files will store the full path of the directories.
# Create emtpy array to store only filenames without path extension.
sample_names=() 
for i in "${!files[@]}"; do
    filename="$(basename "${files[i]}")" # Get only name without path.
    sample_names+=("$filename") # Append the filenames array.
done
# Print the sample names.
echo -e "Sample names: ${sample_names[@]}\n"


# Loop through the barcode names.
for i_sample_name in ${sample_names[@]}; do

# Get the path of the current sample.
i_data="${data_path}/${i_sample_name}"
# Change directory to current sample.
cd ${i_data}
# Get all the file names on the data path.
files=(${i_data}/*) # files will store the full path of the directories.
# Create emtpy array to store only filenames without path extension.
file_names=() 
for i in "${!files[@]}"; do
    filename="$(basename "${files[i]}")" # Get only name without path.
    file_names+=("$filename") # Append the filenames array.
done
# Print the barcode names
echo -e "Sample ${i_sample_name} input barcodes: ${file_names[@]}\n"

# Loop through the barcode names.
for i in ${file_names[@]}; do
# Remove adapters with porechop.
echo -e "Removing adapters from ${i_sample_name} ${i} with porechop...\n"
porechop -i ${i} -o ${i}_${i_sample_name}_adaprm.fastq
done

# Concatenate all amplicon fastqs of one sample into one fastq.
echo -e "Merging all amplicon adapter removed fastqs of ${i_sample_name} into ${i_sample_name}_fastq\n"
cat *${i_sample_name}_adaprm.fastq > ${i_sample_name}.fastq
rm *${i_sample_name}_adaprm.fastq

# Align with minimap2.
echo -e "Aligning ${i_sample_name}.fastq to reference with minimap2...\n"
minimap2 -a ${reference} ${i_sample_name}.fastq > ${i_sample_name}.sam
# Process sam file.
echo -e "Converting sam to bam file...\n"
samtools view -b ${i_sample_name}.sam > ${i_sample_name}.bam
rm ${i_sample_name}.sam
echo -e "Sorting bam file...\n"
samtools sort ${i_sample_name}.bam -o ${i_sample_name}_srtd.bam
rm ${i_sample_name}.bam
# Index the sorted bam file.
echo -e "Indexing bam file...\n"
samtools index ${i_sample_name}_srtd.bam
# Convert conensus fasq from bam file.
echo -e "Making fastq consensus sequence...\n"
samtools mpileup -uf ${reference} ${i_sample_name}_srtd.bam | bcftools call -c | vcfutils.pl vcf2fq > ${i_sample_name}_cns.fq
# Making consensus fasta from consensus fasq...
echo -e "Getting consensus fasta... Bases quality lower than 20 set to N...\n"
seqtk seq -aQ64 -q20 -n N ${i_sample_name}_cns.fq > ${i_sample_name}_cns.fasta
rm ${i_sample_name}_cns.fq ${i_sample_name}_srtd.bam.bai

# Call Rscript to rename the header from the reference to the current sample.
echo -e "Renaming consensus header to ${i_sample_name} \n"
Rscript ${scripts}/change_fasta_header.R ${i_sample_name}_cns.fasta ${i_sample_name} 

# Get depth information to text file and remove bam files.
echo -e "Getting depth from bam file...\n"
samtools depth -a -H ${i_sample_name}_srtd.bam -o ${i_sample_name}_depth.txt

# Make results directory for each sample and move files.
mkdir ${data_path}/${i_sample_name}_results
mv ${i_sample_name}.fastq ${i_sample_name}_srtd.bam ${i_sample_name}_depth.txt ${i_sample_name}_cns.fasta ${data_path}/${i_sample_name}_results
echo -e "Results of ${i_sample_name} were moved to ${data_path}\n"

done

# Change directory back to data path.
cd ${data_path}
# Create empty quality report csv file with headers that will be appended during the following iterations.
echo -e "Creating empty quality_report.csv file...\n"
echo -e "prefix_barcode,average_depth,perc_over_30X" > quality_report.csv

# Loop through the sample names.
for i_sample_name in ${sample_names[@]}; do
# Call the depth rscript to plot depth of each sample and append the quality report.
echo -e "Calling plot_depth.R for ${i_sample_name}...\n" 
Rscript ${scripts}/plot_depth_ONT.R ${i_sample_name}_results/${i_sample_name}_depth.txt `pwd`/quality_report.csv
mv ${i_sample_name}_depth.pdf ${i_sample_name}_low_depth_positions.csv ${i_sample_name}_results


cat ${i_sample_name}_results/${i_sample_name}_cns.fasta >> consensus_genomes.fasta 

done

# Rscript to plot coverage. It will also append the quality report with the coverage of all samples.
echo -e "Calling plot_coverage.R...\n" 
Rscript ${scripts}/plot_coverage_ONT.R `pwd`/consensus_genomes.fasta `pwd`/quality_report.csv
