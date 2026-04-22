# Hufsah Ashraf
# Pipeline to detect toggling-indicating SNPs inside an inverted region using Strand-Seq data and determine whether the resp. locus is recurrent or single-event


import math
from collections import defaultdict
configfile: "Snakefile.json"
path_to_bams = config["path_to_bams"]
path_to_snps = config["path_to_snps"]
path_to_finals = config["path_to_finals"]
output_folder = config["output_folder"]
inversions = config["inversions"]
SAMPLE,ID = glob_wildcards(path_to_bams+"{sample}/{id}.bam")
SAMPLES = sorted(set(SAMPLE))

    
print("Detected {} samples:".format(len(SAMPLES)))


rule all:
    input:
        expand(output_folder+'counts/{sample}/snp_strand_counts.txt', sample = SAMPLES),
        expand(output_folder+'per_sample_configs/{sample}_snp_ref_inv.txt', sample = SAMPLES),
        #expand(output_folder+'per_sample_configs_clean/{sample}_snp_ref_inv.txt', sample = SAMPLES),
        #output_folder+'one_table_per_snp.txt'


rule snp_strand_counts:
    input:
        bam = path_to_bams+'{sample}/',
        bed = inversions,
        vcf = path_to_snps
    output:
        strand_counts = output_folder+"counts/{sample}/snp_strand_counts.txt",
        other_info = output_folder+ "other/{sample}/read_info.txt"
    resources:
    	runtime_hrs = 60 ,
    	mem_total_mb = 10000
    shell:
        "python3 recurrence_first_v2.py -i {input.bam} -b {input.bed} -o {output.strand_counts} -v {input.vcf} -s {wildcards.sample} > {output.other_info} "
               
rule add_cell_states:
    input:
        states = path_to_finals+"{sample}/{sample}_counts.gz",
        snp_counts_file = output_folder+"counts/{sample}/snp_strand_counts.txt"

    output:
    	out = output_folder+'per_sample_configs/{sample}_snp_ref_inv.txt'
    resources:
    	runtime_hrs = 24,
    	mem_total_mb = 150000
    conda: "envs/r4.yaml"
    shell:
    	"""
    	Rscript cell_states_from_final.R \
    	-f {input.states} \
    	-b {input.snp_counts_file} \
    	-o {output.out}\
    	"""
rule clean_configs:
    input:
    	configs = output_folder+'per_sample_configs/{sample}_snp_ref_inv.txt'

    output:
    	out = output_folder+'per_sample_configs_clean/{sample}_snp_ref_inv.txt'
    resources:
    	runtime_hrs = 24,
    	mem_total_mb = 30000
    shell:
    	"""
    	Rscript clean_reads.R
    	"""
rule one_table:
    input:
    	configs=expand("per_sample_configs_clean/{sample}_snp_ref_inv.txt",  sample = SAMPLES),

    output:
    	out = 'one_table_per_snp.txt'
    shell:
    	"""
    	Rscript 1snp_1table_allsamples.R
    	"""




