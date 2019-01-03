#!/usr/bin/python3

''' 
Add "0.0", "1.0" or "." from [chr,pos,base] file to allele-freq file (all biallelic sites)
only for sites present in allele freqs file (no new singletons)
sites where new sample has a third state are removed

'''
	
import sys
import argparse
import gzip



def is_transver(ref, alt):
	'''
	True for transversion
	False for transistion
	'''
	bases = ref + alt
	transi = ("A" in bases and "G" in bases) or ("C" in bases and "T" in bases)
	transver = not transi
	return transver
	
''' test
is_transver("A","T") == True
is_transver("T","A") == True
is_transver("T","C") == False
is_transver("A","G") == False

'''


def merge_base(al, base):
	'''
	in: allele freq line as list, base
	out: merged line or None if triallelic
	'''
	ref = al[3]
	alt = al[4]
	if base == ref:
		out = al + ["0.0"]
	elif base == alt:
		out = al + ["1.0"]
	else:
		out = None
	return out

''' test
merge_base(["21","9828",".","T","G","0.0","1.0"], "T") == ["21","9828",".","T","G","0.0","1.0","0.0"]
merge_base(["21","9828",".","T","G","0.0","1.0"], "G") == ["21","9828",".","T","G","0.0","1.0","1.0"]
merge_base(["21","9828",".","T",".","0.0","1.0"], "T") == ["21","9828",".","T",".","0.0","1.0","0.0"]
merge_base(["21","9828",".","T",".","0.0","1.0"], "G") == None
merge_base(["21","9828",".","T","G","0.0","1.0"], "A") == None

'''
	
	


def main():
	parser = argparse.ArgumentParser(description='Merge an allele-freq file (as created by allele_freqs.py) from STDIN with a file containing [chr | pos | base] (compressed or uncompressed). Positions not in allele-freqs file are skipped, i.e. no new singletons are inserted. Sites that the new sample would make triallelic are removed. Writes to STDOUT.', usage='zcat file1.gz | merge_freqs_base.py base_file.gz newSample --require2')
	parser.add_argument("base_file", help="File with chrom, pos, base; 1-based.")
	parser.add_argument("base_name", help="Name of new column in the allele-freqs header.")
	parser.add_argument("-t", "--transver", action="store_true", help="Transversions only.")
	parser.add_argument("--require2", action="store_true", help="Restrict to positions present in base_file")
	
	args = parser.parse_args()

	# check instream
	if sys.stdin.isatty():
		sys.exit("Error: Needs allele-freqs from stdin, e.g. 'zcat file1.gz | merge_freqs_base.py base_file.gz newSample'")
	al_in = sys.stdin
	
	trans_county = 0
	tri_county = 0
	
	# print header, add new sample to line with donor-ids
	al = al_in.readline()
	sys.stdout.write(al.rstrip()  + "\t" + args.base_name + "\n")

	# open base-file
	if args.base_file.endswith('.gz'):
	    opener = gzip.open
	else:
	    opener = open
	with opener(args.base_file, 'rt') as bases:
		al = al_in.readline().split()
		base = bases.readline().split()
		if al[0] != base[0]:
			sys.exit("Error: Chromosomes do not match: al-chrom={}, base-chrom={}.".format(al[0], base[0]))
		## number of donors
		n_al = len(al)-9
		### iterate over lines until al ends (pos not in al are ignored)
		while not(len(al) == 0):
			try:
				# filter out transitions
				if args.transver and not is_transver(al[3], al[4]):			
					trans_county += 1
					al = al_in.readline().split()
					continue
				
				## merge
				out = None
				## pos not in base, including the case that base file has reached end
				# -> if base not required fill "." read next al line
				if (len(base) == 0) or (int(al[1]) < int(base[1])):
					if not args.require2:
						out = al + ["."]
					al = al_in.readline().split()
				## pos not al --> skip: since REF at this pos is not known
				elif int(base[1]) < int(al[1]):
					base = bases.readline().split()
				## positions are present in both files
				elif al[1] == base[1]:
					out = merge_base(al, base[2])
					# here out is only None if base not in [ref,alt]
					if not out:
						tri_county += 1
					# read next lines from both files
					al = al_in.readline().split()
					base = bases.readline().split()
					
				# print merged line if any
				if out:
					sys.stdout.write("\t".join(out) + "\n")
				
			except (IndexError, ValueError) as errore:
				sys.stderr.write(str(errore) + "\n")
				sys.stderr.write("\t".join(al) + "\n")
				sys.stderr.write("\t".join (base) + "\n")
				sys.exit()
	sys.stderr.write("Nr. of skipped transitions: " + str(trans_county) + "\n")
	sys.stderr.write("Nr. of skipped triallelic sites: " + str(tri_county) + "\n")


if __name__ == "__main__":
	main()

''' test
FREQ=/mnt/sequencedb/gendivdata/4_processed/ALT_freqs_from_merged_VCF/af_high_sgdp_g1000_apes/freq_var_chr21.tab.gz
BASE=/mnt/scratch/steffi/D/random_bases/Denis11/denis11_deam_chr21.tab.gz
MERGE=/mnt/expressions/steffi/D/d_toolbox/merge_freqs_base.py
zcat $FREQ | cut -f-8 | $MERGE $BASE Denis11Deam | less -S
zcat $FREQ | cut -f-8 | $MERGE $BASE Denis11Deam --require2 | less -S
zcat $FREQ | cut -f-8 | $MERGE $BASE Denis11Deam --require2 --transver | less -S

'''


