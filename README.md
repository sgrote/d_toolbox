## tools for D-statistics and file preparation

Here is a quick overview of the functionalities.  
For more details see the help pages, e.g.  
`./d_stats.py -h`  
`./d_genomewide.R -h`

### Common D-statistics

##### `d_stats.py`

Count ABBA and BABA for _D(pop1, pop2, pop3, pop4)_ based on a 
VCF file or allele frequencies as computed by `allele_freqs.py`.  

##### `d_genomewide.R`

Go through all autosome subdirectories and combine ABBA and BABA counts from `d_stats.py`  
to compute  _D(pop1, pop2, pop3, pop4)_.
Perform blockwise jackknife to compute Z-scores.

##### `plot_results.R`

Plot a heatmap of _D(X, Y, pop3, pop4)_ for every _pop3-pop4_ combination.  

##### `plot_results_bar.R`

Plot _D(pop1, pop2, X, pop4)_ for every _pop1-pop2-pop4_ combination as barplot.  

##### `plot_results_hbar.R`

Plot _D(X, Y, pop3, pop4)_ for every _pop3-pop4_ combination as horizontal barplot.

  

### Allele-frequency stratified D-statistics

See `slides/D_stats_freqbins_patterns.pdf` for some notes on the interpretation of stratified D-statistics

##### `d_genomewide_freqbins.R`

compute _D(pop1, pop2, pop3, pop4)_ per B-allele-frequency-bin in _pop3_.  
Introgression from _pop1_ or _pop2_ into _pop3_ should be most visible at the expected 
frequency of the introgressed B-allele in _pop3_.   
Requires prior counting of ABBA and BABA patterns with `d_stats.py` with `--sites 'full'` and `d_genomewide.R`.   

##### `plot_results_freqbins.R`

Plot _D(pop1, pop2, X, pop4)_ per B-allele-frequency-bin for every _pop1-pop2-pop4_ combination as lines or dots.  



### Other D-stats related tools

##### `abba_baba_per_base_combi.R`

compute ABBA and BABA counts for _D(pop1, pop2, pop3, pop4)_ per ancestral/derived (A/B) allele pair.  
Requires prior counting of ABBA and BABA patterns with `d_stats.py` with `--sites 'full'` and `d_genomewide.R`.   

##### `plot_abba_baba_per_base_combi.R`

plot ABBA and BABA counts for _D(pop1, pop2, pop3, pop4)_ per ancestral/derived (A/B) allele pair.     

##### `count_AAAA_BAAA_ABAA.py`

counts AAAA, BAAA and ABAA patterns for all combinations of _pop1-pop4_ (transitions and transversions separately).  
Useful for quickly checking excess of B-alleles in _pop1_ or _pop2_, potential errors.  




### VCF tools

##### `allele_freqs.py`

Calculate allele frequencies per population for biallelic sites.  
The output is similar to a VCF-file and can be used as alternative input for `d_stats.py` (instead of a VCF).  
`d_stats.py` is much faster with precomputed allele frequencies.  

##### `merge_vcf_vcf.py`

Merge two VCF files, one from stdin the other as first command-line argument (.vcf or .vcf.gz).  
Do filtering on the fly, write to stdout for multiple merging steps.  


##### `merge_vcf_base.py`

Merge a VCF file from stdin with a file containing [chr | pos | base | (base)],  
like created by `bam_ran_base.py`.  
Bases get converted to pseudo-genotypes.  

##### `base_to_gt.py`

Take a VCF file from stdin that has bases in LAST columns instead of genotype,  
convert bases to pseudo-genotype.  
Use case: Kays `BamSNPAddMaf` was used to add ape-bases to VCF.

##### `vcf_filter.py`

Filter one VCF file from stdin. Acts on whole lines and/or individual genotypes.  

##### `filter_vcf_with_bed.py`

Filter VCF-like file from stdin with bed-file.  
Also works with allele frequency files created by `allele_freqs.py`.  

##### `subsample_vcf.py`

Subsample streamed in vcf-like file randomly given the proportion of passing sites.  
E.g. `zcat file.vcf.gz | subsample_vcf.py 0.5` will remove around half of the positions from `file.vcf.gz`.

##### `remove_group_singletons.py`

Remove lines from VCF that are only variable due to alleles in a group of samples.  
Mask genotypes that are private to the group in lines that are kept due to other variation.  
Use case: Remove genotypes from VCF that have mutations only observed in low-coverage genomes (potential errors).  

##### `count_triallelic_by_sample.py`

Count sites in VCF that are triallelic only due to a group of specified samples.  
Use case: How many biallelic sites got lost, because some other sample was added and introduced a third allele? 

##### `haploidize_male_X.py`

Take VCF from stdin and an infofile with sex per sample,  
and convert diploid genotypes for males outside the PAR region to haploid genotypes.  
E.g. `'0/0' -> '0'`, `'1/1' -> '1'`, `'1/0' -> '.'`



### Other

##### `bam_ran_base.py`

Take a _samtools mpileup *.bam_ output from stdin and sample a random read.  
Write [chr | pos | base] to outfile (one file per chrom).  

##### `remove_transitions.py`

Restrict to biallelic transversions given index of REF and ALT column (can be applied to several file-formats).

##### `remove_allele_freqs.py`

Restrict to sites in range of allele-frequencies given index of allele-freq column
(can be applied to several file-formats).  

##### `pw_diverge.py`

Calculate pairwise divergences between populations given a VCF file.





