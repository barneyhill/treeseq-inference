#!/usr/bin/env python3
"""
Take a vcf file with ancestral information (e.g. from 1000 genomes phase 1, such as
ALL.chr22.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz
and convert it to the huge samples X sites array used in tsinfer. Store this
in the hdf5 output file, and also store the variant positions and the sample names.

vcf2tsinfer.py 1000G_chr22.vcf.gz outputarrays.hdf5

NB - a slightly odd observation: all the alleles for which the ancestral state is not present in the dataset are SNPs

read using
with h5py.File(filename, 'r') as f:
    print(f['data']['variants'][()])
"""
import sys
import collections

import numpy as np
import h5py
import pysam
#To do - use argparse module 

vcf_in = pysam.VariantFile(sys.argv[1])

allele_count = {}
rows, row = {}, 0
for sample_name in vcf_in.header.samples:
    for suffix in ('a','b'):
        rows[sample_name+suffix]=row
        row+=1
position = collections.OrderedDict()
n_vars = 400000 #fill in from https://groups.google.com/forum/#!topic/pysam-user-group/-Obv5ggp9tk
sites_by_samples = np.zeros((n_vars, len(rows)), dtype="i1")
output_freq_variants = 1e3 #output status after multiples of this many variants read

def process_variant(rec):
    """
    return True only when
    """
    allele_count[len(rec.alleles)] = allele_count.get(len(rec.alleles),0) + 1
    #restrict to cases where ancestral state contains only letters ATCG
    # i.e. where the ancestral state is certain (see ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/analysis_results/supporting/ancestral_alignments/human_ancestor_GRCh37_e59.README
    if "AA" in rec.info and all(letter in "ATCG" for letter in rec.info["AA"]):
        #only use biallelic variants
        if len(rec.alleles) == 2:
            if rec.info["AA"] not in rec.alleles:
                print("Ancestral state {} not in allele list ({}) for position {}".format(\
                    rec.info["AA"], rec.alleles, rec.pos))
                pass
            else:
                omit=False
                #check for duplicate positions, often causes e.g. by C/CAT as opposed to -/AT
                #(if the first x letters are the same, we have the wrong position
                for allele_start in range(min(len(rec.alleles[0]), len(rec.alleles[1]))):
                    if rec.alleles[0][allele_start]!=rec.alleles[0][allele_start]:
                        break
                if allele_start != 0:
                    pos = rec.pos+allele_start
                    print("Starting allele at an incremented position ({} not {}) for {} (alleles:{})".format(
                        pos, rec.pos, rec.id, rec.alleles))
                else:
                    allele_start=0
                    pos = rec.pos

                if rec.pos in position:                    
                    print("More than one set of variants at position {}. ".format(rec.pos))
                    print("Previous was {}. Omitting a subsequent duplicate (id {})".format(
                        position[rec.pos], {rec.id:rec.alleles}))
                    return False

                column = np.zeros((len(rows),), dtype="i1")
                for label, samp in rec.samples.items():
                    for i,suffix in enumerate(('a','b')):
                        #print("allele {} (pos {}, sample {}, ancestral state {}, alleles {})".format( \
                        #    rec.alleles[i], rec.pos, label+suffix, rec.info["AA"], rec.alleles))
                        if samp.alleles[i] not in rec.alleles:
                            print("@ position {}{}, sample allele {} is not in {} - could be missing data. Omitting this row".format(
                                rec.pos, "+{}".format(allele_start) if allele_start else "", samp.alleles[i], rec.alleles))
                            return False
                        column[rows[label+suffix]] = samp.alleles[i][allele_start:]!=rec.info["AA"][allele_start:]
                
                sites_by_samples[len(position)]=column
                position[pos]={rec.id:rec.alleles}
                return True
    return False

for j, variant in enumerate(vcf_in.fetch()):
    process_variant(variant)
    if j % output_freq_variants == 0:
        print("{} variants read ({} saved). Base position {} Mb (alleles per site: {})".format(j+1, len(position), variant.pos/1e6, [(k, allele_count[k]) for k in sorted(allele_count.keys())]))

with h5py.File(sys.argv[2], "w") as f:
    g = f.create_group("data")
    g.create_dataset("position", data=position.keys())
    g.create_dataset("samples", data=[s.encode() for s in sorted(rows, key= rows.get)])
    g.create_dataset("variants", data=np.transpose(sites_by_samples[:len(position)]))
print("Saved {} biallelic loci for {} samples into {}".format(len(position), len(rows), sys.argv[2]))
