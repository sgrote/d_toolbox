#!/usr/bin/python3

'''
take a samtools mpileup *.bam output from stdin
and sample a random read
write [chr | pos | base] to outfile (one file per chrom)
'''

'''
input
['1', '86329', 'N', '1', 'G', ']']
['1', '86330', 'N', '1', 'G', ']']
['1', '86331', 'N', '1', 'A', ']']
['1', '86332', 'N', '1', 'T', ']']
'''

'''
Mpileup
A pattern `\\+[0-9]+[ACGTNacgtn]+' indicates there is an insertion between this reference position and the next reference position. The length of the insertion is given by the integer in the pattern, followed by the inserted sequence. Similarly, a pattern `-[0-9]+[ACGTNacgtn]+' represents a deletion from the reference. The deleted bases will be presented as `*' in the following lines. Also at the read base column, a symbol `^' marks the start of a read. The ASCII of the character following `^' minus 33 gives the mapping quality. A symbol `$' marks the end of a read segment. 
'''


import sys
import random
#import imp
import gzip
import argparse

###########################################################

# one position
def sample_base(pos, min_cov, min_qual, n_reads):
    ## remove 0-read-coverage and Indels (not a full line)
    if int(pos[3]) < min_cov or "+" in pos[4] or "-" in pos[4]:
        #print("aussortiert")
        #print(pos)
        return None
    bases = list(pos[4].upper())
    quals = list(pos[5])
    
    # remove read-start/end encodings
    if "^" in bases or "$" in bases:
        bases = [bases[i] for i in range(len(bases)) if not (bases[i] in "$^" or bases[i-1]=="^")]
        #print("clean bases")
        #print(bases)
            
    # filter for quality, also remove "*" which has qual associated and therefor cannot be removed before
    q = [ord(qual)-33 for qual in quals]
    if any([qu < min_qual for qu in q]) or any([base == "*" for base in bases]):
        #print("filtere!")
        #print(pos)
        #print(q)
        #print([qu < min_qual for qu in q])
        bases = [bases[i] for i in range(len(bases)) if q[i] >= min_qual and bases[i]!="*"]
        if len(bases) == 0:
            return None
    
    if any([base not in ["A","C","T","G"] for base in bases]):
        print("unclean!")
        print(pos)
        print(bases)
        bases = [b for b in bases if b in ["A","C","T","G"]]
        
    # check coverage again after filtering
    if len(bases) < min_cov:
        return None
    
    # select a random base
    try:
        ran_base = random.sample(bases, n_reads)
    except ValueError as errore:
        print(bases)
        print(n_reads)
        print(min_cov)
        sys.exit(errore)
    out = pos[:2] + ran_base
    return out 


''' test
pos = ["21", "9579965", "N", "4", "ACTGE", "]A\]]"]
sample_base(pos, 3, 30, 3)
pos = ["21", "9579965", "N", "4", "ACTG", "???N"]
sample_base(pos, 1, 30, 1)
# only one passing quality
sample_base(pos, 1, 40, 1) == ['21', '9579965', 'G']
# remove read-range encodings
pos = ["21", "9579965", "N", "4", "^FA$", "N"]
sample_base(pos, 1, 30, 1) == ['21', '9579965', 'A']
sample_base(pos, 2, 30, 2) == None

'''

###########################################################


def main():
    parser = argparse.ArgumentParser(description="take a samtools mpileup *.bam output from stdin and sample one or more random bases. write [chr | pos | base1 | base2 ...] to outfile (one file per chrom)", usage="samtools mpileup sample.bam | bam_ran_base.py outfile_trunk")
    parser.add_argument("outfile_trunk", type=str, help="trunk of output filename (one file per chrom will be generated)")
    parser.add_argument("--min_cov", type=int, default=1, help="minimum number of filtered-passing reads covering a site. Site will be skipped if lower coverage.")
    parser.add_argument("--min_qual", type=int, default=30, help="minimum quality of base, bases with lower qual will be disregarded.")
    parser.add_argument("--n_reads", type=int, default=1, help="number of bases to be sampled (without replacement).")
    
    args = parser.parse_args()
    
    # check coverage filter
    if args.n_reads > args.min_cov:
        args.min_cov = args.n_reads
    
    # mpileup from stdin
    if sys.stdin.isatty():
        sys.exit("Error: Needs mpileup from stdin, e.g. 'samtools mpileup sample.bam | bam_ran_base.py outfile_trunk'")
    mpileup = sys.stdin
    
    # go through file
    pos = mpileup.readline().split()
    i = 1
    while True:
        chrom = pos[0]
        print("next chrom!")
        outfile = args.outfile_trunk + "_chr" + chrom + ".tab.gz"
        with gzip.open(outfile, "wt") as out:
            outline = sample_base(pos, args.min_cov, args.min_qual, args.n_reads)
            if outline:
                outline = "\t".join(outline)+"\n"
                out.write(outline)
            for pos in mpileup:
                pos = pos.split()
                ## next chrom?
                if pos[0] != chrom:
                    print("Next chrom!")
                    break
                outline = sample_base(pos, args.min_cov, args.min_qual, args.n_reads)
                if outline:
                    outline = "\t".join(outline)+"\n"
                    out.write(outline)
                    i += 1
                    if i % 100000 == 0: 
                        print(pos)
                        print(outline.rstrip())
            else:
                break


if __name__ == "__main__":
    main()



