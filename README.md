# GenomePermute
Permutation analysis to investigate the significance of overlap between genomic features

This perl script does permutation analysis to calculate significance of enrichment or depletion of overlap between genomic regions given in two bed files.

Usage example: perl ./permute_human.pl GRCh37_ens75_genes.bed GRCh37_G4_1-3.bed

This will check overlap between gene bodies annotated in release 75 of Ensembl database (http://www.ensembl.org) and predicted G4 quadruplex structures in GRCh37 reference genome that have following motif:  G3+ N1-3 G3+ N1-3 G3+ N1-3 G3+ )
The example input files, chromosome lengths and gap locations for GRCh37 assembly are part of this repository

Expected output:
=== FINAL CONCLUSION: Enriched, p<0.001 ===

Program also generates an output file with number of overlaps dicovered for every permutation. This file can be used for generation of boxplot or violin plots that show value of non-permuted overlap value and distribution of permuted values.
