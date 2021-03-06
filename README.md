# tiSNPs Detection

## Pipeline to analyse inversions for recurrence by identifying SNPs discrepant with a single event origin.

This workflow uses [Snakemake](https://bitbucket.org/snakemake/snakemake) to execute all steps in order. Before running the Snakefile:

  1. the alignmnet (.bam) files should be placed in '.../bam/{{sample}}/all/' and the path should be provided to 'path_to_bams'
  2. path to bi-allelic SNPs belonging to chr1-chrX should be provided in 'path_to_snps_1toX' and path to bi-allelic SNPs belonging to chrY should be provided in 'path_to_snps_Y'
  3.  files with background cell state for all samples via [mosaicatcher](https://github.com/friendsofstrandseq/mosaicatcher), (example provided in sample_files/) should be placed as '{path}/{{sample}}/100000_fixed_norm.selected_j0.1_s0.1/final.txt' and path should be added to 'path_to_finals'
  4.  Inversions belonging to chr1-chrX to be tested should be placed in a bed file (sample provided in sample_files/) 'invs_1toX.bed' and those from chrY in 'invs_Y.bed'
  
A successful run of the pipeline would generate a file 'one_table_per_snp.txt' which maintains a record of haplotype configuration counts and some other statistics, for each within inversion SNP, a file 'sig_nosig_invs.txt' with recurrence label 'signal' for each inversion depending on the aggregated SNP evidence across the whole locus and a file 'sig_invs.txt' with a subset of the file 'sig_nosig_invs.txt', containing only the inversions where a recurrence signal was observed. 
