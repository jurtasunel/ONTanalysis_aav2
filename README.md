# ONTanalysis_aav2

This pipeline maps ONT reads from tile-amplicons to a reference with minimap2. The outputs are consensus.fasta files and depth.pdf files for each sample, a quality_report.csv and a coverage.png file of all samples.

"ONTanalysis_aav2.sh" is the parent script to run, the pipeline goes as follows:
- Get all sample names and barcodes for each sample and remove adapters with porechop.
- Merge all adapter-removed fastqs of each sample into a single sample.fastq.
- Align sample.fastq with minimap2, then make sam > bam > srt_bam > consensus.fasta.
- Call Rscript "change_fasta_header.R" to rename consensus fasta header to be the sample name (by default is the reference genome header name).
- Call Rscript "plot_depth_ONT.R" to produce depth.pdf plots for each sample.
- Merge all samples consensus into a single "consensus.fasta" and call Rscript plot_coverage_ONT.R to produce coverage.pdf plot of all samples.
- Make a quality_report.csv

  ### On the data path should ONLY be directories with the Sample names, and within each sample directory should ONLY be barcodes directories containing fastq.gz files
