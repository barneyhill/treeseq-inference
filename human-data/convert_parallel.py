import sys
import os
import multiprocessing
import pysam
import cyvcf2
import tsinfer
import numpy as np
import tqdm
import subprocess

num_no_ancestral_state = 0
num_low_confidence_ancestral_state = 0

def get_ancestral_state(ancestral_states, position):
    ret = None
    ancestral_state = ancestral_states[position]
    if ancestral_state in [".", "N", "-"]:
        pass #num_no_ancestral_state += 1
    elif ancestral_state.lower() == ancestral_state:
        pass #num_low_confidence_ancestral_state += 1
    else:
        assert ancestral_state in ["A", "C", "T", "G"]
        ret = ancestral_state
    return ret


def add_diploid_sites(vcf, samples, ancestral_states, chromosome=None, start=None, stop=None):
    """
    Read the sites in the vcf and add them to the samples object, reordering the
    alleles to put the ancestral allele first, if it is available.
    """
    pos = 0
    if chromosome is None or start is None or stop is None:
        vcf_iter = vcf
        start = pos = 0  # VCFs are 1-based, so 0 is not in the first file
    else:
        # This is a partial file: use the standard HTSlib string spec
        vcf_iter = vcf(f"{chromosome}:{start}-{stop}")
        pos = start - 1
    dupe = 0
    no_anc = 0
    added = 0
    #length = len([x for x in vcf_iter])
    for variant in vcf_iter:  # Loop over variants, each assumed at a unique site
        if added % 100 == 0: print(dupe,no_anc,added, multiprocessing.current_process())
        if variant.POS < start:
            # This could be a long variant from a previous section of the VCF that
            # extends into this file. It should have been capture previously, so we skip
            continue
        if pos == variant.POS:
            dupe += 1 #print("dupe")
            continue
            raise ValueError("Duplicate positions for variant at position", pos)
        else:
            pos = variant.POS
        #if any([not phased for _, _, phased in variant.genotypes]):
        #    print("unphased")
        #    continue
        #    raise ValueError("Unphased genotypes for variant at position", pos)
        
        #print("phased")
        alleles = [variant.REF] + variant.ALT
        ancestral = get_ancestral_state(ancestral_states, variant.POS)
        if ancestral == None: 
            no_anc +=1
            #print(multiprocessing.current_process(), ancestral)
            continue
        #print(multiprocessing.current_process(), ancestral)
        # Ancestral state must be first in the allele list.
        ordered_alleles = [ancestral] + list(set(alleles) - {ancestral})
        allele_index = {old_index: ordered_alleles.index(allele)
            for old_index, allele in enumerate(alleles)}
        # Map original allele indexes to their indexes in the new alleles list.
        genotypes = []
        for row in variant.genotypes:
            if row[2] == False: continue
            for old_index in row[0:2]:
                genotypes.append(allele_index[old_index])
        #genotypes = [allele_index[old_index]
        #    for row in variant.genotypes for old_index in row[0:2]]
        #print("PHASED")
        samples.add_site(pos, genotypes=genotypes, alleles=ordered_alleles)
        added+=1
    print("TERMINATING", dupe,no_anc,added, multiprocessing.current_process())
def read_vcf_part(params):
    ancestral_file = "homo_sapiens_ancestor_GRCh38/homo_sapiens_ancestor_20.fa" 
    fasta = pysam.FastaFile(ancestral_file)
    # NB! We put in an extra character at the start to convert to 1 based coords.
    ancestral_states = "X" + fasta.fetch(reference=fasta.references[0])
   
    vcf = cyvcf2.VCF(params[0])
    vcf.set_index(params[1])
    sequence_length = vcf.seqlens[vcf.seqnames.index("chr20")]
    print(sequence_length)
    part = params[2]
    num_parts = params[3]
    assert sequence_length > num_parts
    chunk_size = sequence_length // num_parts + 1
    start = part*chunk_size
    stop = (part+1)*chunk_size - 1  # ranges are inclusive
    stop = min(stop, sequence_length)
    
    print(start,stop,sequence_length,multiprocessing.current_process())
    
    path = f"parts/partial_file{start}-{stop}.samples"
    try:
        with tsinfer.SampleData(
            path=path,
            max_file_size = 30 * 2**30,
            sequence_length=sequence_length) as samples:
                add_diploid_sites(vcf, samples, ancestral_states, vcf.seqnames[vcf.seqnames.index("chr20")], start, stop)
    except: return None
    print("FINISHED:", samples.num_individuals, samples.num_samples, samples.num_sites)

    if samples.num_sites == 0:
        samples.close()
        return None 
    
    samples.close()  # No need to leave the file open

    return path

if __name__ == '__main__':

    # URL for the VCF
    url = "/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb/data/phased/wes_union_calls/with_parents/ukb_eur_wes_union_calls_200k_chr20.vcf.bgz"
    idx = "/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb/data/phased/wes_union_calls/with_parents/ukb_eur_wes_union_calls_200k_chr20.vcf.bgz.tbi"
    num_data_sites = int(subprocess.check_output(["bcftools", "index", "--nrecords", url]))
 
    num_threads = int(sys.argv[1])
    with multiprocessing.Pool(num_threads) as pool:
        sample_filenames = [r for r in pool.map(
            read_vcf_part, ((url, idx, i, num_threads) for i in range(num_threads)))]
    all_samples = [tsinfer.load(name) for name in sample_filenames if name != None]
    samples = all_samples[0].copy()
    samples.append_sites(*all_samples[1:])
    samples.finalise()
    print(
        "Created full samples file with",
        samples.num_sites, "sites and",
        samples.num_samples, "samples")
