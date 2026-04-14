#Hufsah Ashraf, 9.11.2020
#latest script for getting (per single cell) Watson and Crick read counts for each bi allelic SNP allele that lies within an inversion
import argparse
import sys
import pysam
from pathlib import Path
import glob
import pdb
from collections import defaultdict
import vcf
from pysam import VariantFile

def counts(input_bam, input_bed,count_output, input_vcf,sample):
    dictionary={}
    #change 'a' to 'w' in following line if the file doesn't already exist, 'a' appends to an already existing file
    counts_file= open(count_output, 'w')
    # Write header to output file
    counts_file.write('sample'+'\t'+'cell'+'\t'+'chrom'+'\t'+'inv_start'+'\t'+'inv_end'+'\t'+'inv_ID'+'\t'+'snp_pos'+ "\t"+'ref_allele'+ "\t" +'alt_allele'+ "\t"+'seen_allele'+ "\t" + 'orientation' + "\t" + 'count' + "\t" + 'snp_AF' +"\n")
    # Get all bam file paths
    path = Path(input_bam)
    glob_path = path.glob('*.bam')

    # Iterate over each cell / bam file
    for file in glob_path:
        # Load the according file
        file_name = str(file).strip().split("/")[-1]
        cell=file_name.strip().split(".bam")[0]
        bam_file = pysam.AlignmentFile(file, "rb")
        # Also get the bed_file
        with open(input_bed, 'r') as bed_file:
            # Iterate over each bed file entry
            next(bed_file)#skip the header
            for line in bed_file:
                line_r = line.strip().split("\t")
                chromosome= line_r[0]
                interval_start= int(line_r[1])
                interval_end= int(line_r[2])
                inv_length= interval_end-interval_start
                inversion_id=chromosome+'-'+str(interval_start)+'-INV-'+str(inv_length)
                # Fetch all bam entries that fall in this bin.
                vcf_in = VariantFile(input_vcf)  # auto-detect input format
                recs=[]
                for variant in vcf_in.fetch(chromosome, interval_start, interval_end):
                    assert len(variant.alts)==1 #to make sure all snps are bi allelic
                    #if variant.info["AF"][0] >=0.2 and variant.info["AF"][0] <=0.8: #to consider only frequent snps
                    recs.append((variant.pos,variant.ref, variant.alts[0], variant.info["AF"][0]))


                for read in bam_file.fetch(chromosome, interval_start, interval_end):
                    orientation= None
                    allele=None
                    readname= read.query_name
                    query_seq = read.query_sequence
                    query_seq_index= (read.query_alignment_start)-1
                    if read.reference_start is not None:
                        coord_ref_start = (read.reference_start)+1 # leftmost ref_coordinate where the read maping started
                    if read.reference_end is not None:
                        coord_ref_end = (read.reference_end)-1 #rightmost ref_coordinate where the read maping ended
                    readinfo=read.get_aligned_pairs(matches_only=False, with_seq=False)
                    if len(readinfo)>0:
                        print(readname, readinfo[0])
                    all_tuples=list()
                    i=0
                    snp_count=0 #no. of snps a read covers
                    index_counter=0
                    c1 = 'read.is_read2'
                    c2 = 'read.is_qcfail'
                    c3 = 'read.is_secondary'
                    c4 = 'read.is_duplicate'
                    c5 = 'read.mapq < 10'
                    c6 = 'read.pos < interval_start'
                    c7 = 'read.pos >= interval_end'
                    if any([eval(c1),eval(c2),eval(c3),eval(c4),eval(c5),eval(c6),eval(c7)]):
                        pass
                    else:
                        if read.is_reverse:
                            orientation='Watson'
                        elif not read.is_reverse:
                            orientation='Crick'
                        else:
                            orientation=None

                        if orientation is not None:
                            seen_positions = set()
                            for r in range(len(recs)):
                                rec_pos,rec_ref,rec_alt,AF= recs[r]
                                assert rec_pos not in seen_positions
                                seen_positions.add(rec_pos)
                                # a variant should be processed only if it lies somehere within the read mapping region
                                if rec_pos>=coord_ref_start and rec_pos<=coord_ref_end:
                                    index_all_tuples=index_counter
                                    for t in range(len(readinfo)):
                                        snp_count+=1
                                        tup=readinfo[t]
                                        index_counter+=1
                                        if tup[1] is not None:
                                            all_tuples.append(tup[0])
                                            i+=1
                                            if tup[1]==rec_pos:
                                                posit= all_tuples[i-2]
                                                if posit is not None:
                                                    allele= query_seq[posit]
                                                    print('allele was obtained from position', posit , 'in the query and all tuples was accessed at ', i-2 )
                                                else:
                                                    allele=None
                                                print(cell,readname,coord_ref_start, coord_ref_end,inversion_id,len(recs),rec_pos,rec_ref, rec_alt, allele, orientation)
                                                if allele is not None and (allele==rec_ref or allele==rec_alt):
                                                    if (sample, cell,chromosome, interval_start, interval_end, inversion_id,rec_pos,rec_ref,rec_alt,allele,orientation) in dictionary:
                                                        dictionary[(sample, cell,chromosome, interval_start, interval_end, inversion_id,rec_pos,rec_ref,rec_alt,allele,orientation,AF)]+=1
                                                    else:
                                                        dictionary[(sample, cell,chromosome, interval_start, interval_end, inversion_id,rec_pos,rec_ref,rec_alt,allele,orientation,AF)]=1

                                                else:
                                                    print(rec_pos,rec_ref, rec_alt, read.query_name, 'neither of the two alleles but carries ', allele)
                                                break
                                #if it is before the read start, move to the next variant
                                elif rec_pos<coord_ref_start:
                                    pass
                                #if it lies after the read end, move to the next read
                                elif rec_pos>coord_ref_start:
                                    print(cell, readname, 'covered ',snp_count, 'snps')
                                    break
    for line in dictionary:

        counts_file.write(str(line[0])+ "\t"+ str(line[1])+ "\t" + str(line[2])+"\t" + str(line[3])+ "\t"+ str(line[4])+ "\t" +
        str(line[5]) + "\t"+ str(line[6])+"\t"+str(line[7])+ "\t"+ str(line[8])+"\t"+str(line[9])+"\t"+str(line[10])+"\t" + str(dictionary[line])+"\t"+str(line[11])+"\n")



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-i", "--input_bam", help="The input bam files path", required=True)
    parser.add_argument("-b", "--input_bed", type=str, help="The bed file with inversion segments",required=True)
    parser.add_argument("-v", "--input_vcf",type=str, help="The vcf file containing only bi allelic SNVs, sorted by bcftools", required=True)
    parser.add_argument("-o", "--count_output",type=str, help="The output file with watson and crick reads for each SNP allele (one per sample)", required=True)
    parser.add_argument("-s", "--sample", help="The sample name", required=True)
    args = parser.parse_args()
    counts(args.input_bam, args.input_bed, args.count_output, args.input_vcf, args.sample)
