#!/usr/bin/python3

'''
restrict to sites in range of allele-frequencies given index of allele-freq column
(can be applied to several file-formats)
take file from stdin
write to stdout
'''

import sys
import argparse


def main():
	parser = argparse.ArgumentParser(description='restrict to sites in range of allele-frequencies given index of allele-freq column (can be applied to several file-formats). Reads from stdin and writes to stdout.', usage='cat freq_file | remove_allele_freqs.py -i 3 --min 0.3 --max 0.7 > freq_file_filtered')
	# mandatory
	parser.add_argument("-i","--freq_index", type=int, default=5, help="Column number of allele freqs (0-based).")
	parser.add_argument("--min", type=float, default=0.3, help="Minimum allele freq (inclusive).")
	parser.add_argument("--max", type=float, default=0.7, help="Maximum allele freq (inclusive).")

	args = parser.parse_args()

	# infile from stdin
	if sys.stdin.isatty():
		sys.exit("Error: Needs file from stdin.")
	infile = sys.stdin

	for line in infile:
		try:
			fields = line.rstrip().split()
			# frequency
			freq = float(fields[args.freq_index])
			# reduce to freq-range 
			if freq >= args.min and freq <= args.max:
				sys.stdout.write(line)
		except (IndexError, ValueError) as errore:
			sys.stderr.write(str(errore) + "\n")
			sys.stderr.write(line + "\n")
			sys.exit()

if __name__ == "__main__":
	main()





