#!/usr/bin/python3

'''
remove samples from vcf or vcf.gz file
'''

import gzip
import sys
import argparse
import warnings

def main():
    parser = argparse.ArgumentParser(description='Remove samples from vcf or vcf.gz file, write to new compressed file. ALT-alleles are not updated.', usage='remove_samples_from_vcf.py input.vcf sample1,sample2,sample3 output.vcf.gz')
    parser.add_argument("vcf_file", type=str, help="vcf file to remove samples from (.vcf or .vcf.gz)")
    parser.add_argument("samples", type=str, help="comma-separated list of samples to remove")
    parser.add_argument("out_file", type=str, help="outfile (.vcf.gz)")
    
    args = parser.parse_args()
    print(args)
    
    i = 0
    
    ### open vcf and outfile
    # allow zipped and unzipped vcf
    if args.vcf_file.endswith('.gz'):
        opener = gzip.open
    else:
        opener = open
    with opener(args.vcf_file, 'rt') as vcf:
        # tmp-outfile 
        with gzip.open(args.out_file, 'wt') as out:
        
            ### read input strings
            samples = args.samples.split(",")
            print(samples)

            ### skip header until #CHROM
            # skip vcf_1 header, only keep line with donor-ids
            v1 = vcf.readline()
            while v1[:6] != "#CHROM":
                out.write(v1)
                v1 = vcf.readline()
            
            ### find columns matching input strings (warn about not-in)
            sample_header = v1.rstrip().split()
            samples_not_in_header = [s for s in samples if s not in sample_header]
            if len(samples_not_in_header) > 0:
                warnings.warn("samples not in vcf: " + ", ".join(samples_not_in_header))
            if len(samples_not_in_header) == len(samples):
                sys.exit("None of the input samples is present in vcf")
            print(samples_not_in_header)
            keep_indi = [i for i in range(len(sample_header)) if sample_header[i] not in samples]
            
            ### print remove from header and print rest for all other lines
            while v1:
                v1 = v1.rstrip().split()
                out_line = [v1[i] for i in keep_indi]
                out.write("\t".join(out_line) + "\n")
                if i % 100000 == 0:
                    print(i)
                    print(v1[:12])
                i = i+1
                v1 = vcf.readline()



if __name__ == "__main__":
	main()



