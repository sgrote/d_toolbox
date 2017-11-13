#!/usr/bin/python

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
def sample_base(pos, min_reads, min_qual):
	## remove 0-read-coverage and Indels (not a full line)
	if int(pos[3]) < min_reads or "+" in pos[4] or "-" in pos[4]:
		#print "aussortiert"
		#print pos
		return(None)
	bases = list(pos[4].upper())
	quals = list(pos[5])

	# remove read-start/end encodings
	if "^" in bases or "$" in bases:
		bases = [bases[i] for i in range(len(bases)) if not (bases[i] in "$^" or bases[i-1]=="^")]
		#print "clean bases"
		#print bases
			
	# filter for quality, also remove "*" which has qual associated and therefor cannot be removed before
	q = [ord(qual)-33 for qual in quals]
	if any([qu < min_qual for qu in q]) or any([base == "*" for base in bases]):
		#print "filtere!"
		#print pos
		#print q
		bases = [bases[i] for i in range(len(bases)) if q[i] >= min_qual and bases[i]!="*"]
		#print bases
		if len(bases) == 0:
			return(None)
		
	if any([base not in "ACTG" for base in bases]):
		print "unclean!"
		print pos
		print bases
	
	# select a random base
 	ran_base = random.choice(bases)
 	#if len(set(bases)) > 1: 
		#print pos
		#print ran_base
	out = pos[:2] + [ran_base]
	return(out)



###########################################################


def main():
    parser = argparse.ArgumentParser(description="take a samtools mpileup *.bam output from stdin and sample a random read write [chr | pos | base] to outfile (one file per chrom)", usage="samtools mpileup sample.bam | bam_ran_base.py outfile_trunk")
    parser.add_argument("outfile_trunk", type=str, help="trunk of output filename (one file per chrom will be generated)")
    parser.add_argument("--min_reads", type=int, default=1, help="minimum number of reads, site will be skipped if lower coverage.")
    parser.add_argument("--min_qual", type=int, default=30, help="minimum quality of base, bases with lower qual will be disregarded.")
    
    args = parser.parse_args()
    
    # mpileup from stdin
    if sys.stdin.isatty():
        sys.exit("Error: Needs mpileup from stdin, e.g. 'samtools mpileup sample.bam | bam_ran_base.py outfile_trunk'")
    mpileup = sys.stdin
    
    # go through file
    pos = mpileup.readline().split()
    i = 1
    while True:
        chrom = pos[0]
        print "next chrom!"
        #if chrom=="3": break
        outfile = args.outfile_trunk + "_chr" + chrom + ".tab.gz"
        print pos
        with gzip.open(outfile, "wb") as out:
            outline = sample_base(pos, args.min_reads, args.min_qual)
            print outline
            if outline:
                outline = "\t".join(outline)+"\n"
                out.write(outline)
            for pos in mpileup:
                pos = pos.split()
                ## next chrom?
                if pos[0] != chrom:
                    print "Next chrom!"
                    break
                outline = sample_base(pos, args.min_reads, args.min_qual)
                #print outline
                if outline:
                    outline = "\t".join(outline)+"\n"
                    out.write(outline)
                    i += 1
                    if i % 100000 == 0: 
                        print pos
                        print outline.rstrip()
            else:
                break


if __name__ == "__main__":
	main()



