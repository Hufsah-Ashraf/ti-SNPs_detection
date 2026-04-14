# tiSNPs Detection

## Pipeline to analyse inversions for recurrence by identifying SNPs discrepant with a single event origin.

This workflow uses [Snakemake](https://bitbucket.org/snakemake/snakemake) to execute all steps in order. Before running the pipeline, the file paths should be added to Snakefile.json:

  1. path_to_bams: the folder with alignmnet (.bam) files for each sample. The folder structure should be "path_to_bams/{sample}/{one bam file per cell for this sample}"
  2. path_to_snps: path to the vcf containig bi-allelic SNPs
  3. path_to_finals: files with background cell state for all samples via [mosaicatcher](https://github.com/friendsofstrandseq/mosaicatcher), (example provided in sample_files/). Folder structure should be "path_to_finals/{sample}/{sample}_counts.gz"
  4. inversions: bed file wiht inversions to be tested for recurrence (example provided in sample_files/) 'invs.bed'
  5. output_folder: path to the folder where all the output files should be written
  
A successful run of the pipeline would generate a file 'one_table_per_snp.txt' which maintains a record of haplotype configuration counts and some other statistics, for each within inversion SNP, a file 'sig_nosig_invs.txt' with recurrence label 'signal' for each inversion depending on the aggregated SNP evidence across the whole locus and a file 'sig_invs.txt' with a subset of the file 'sig_nosig_invs.txt', containing only the inversions where a recurrence signal was observed. 
