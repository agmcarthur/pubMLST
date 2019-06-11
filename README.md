## PubMLST Script

Entering *pubmlst* without arguments will give a help screen:

```bash
pubmlst
	PubMLST will analyze all .fasta files in the current directory.
	First argument is name of analysis, used in output file names.
	Second argument is the path to the directory containing PubMLST allele and type files
```

Example analysis:

```bash
pubmlst myEfaecalis ~/db/Enterococcus_faecalis_PubMLST
```

Updated or additional reference files can be obtained at https://pubmlst.org/databases/. There is no guarantee the copies here are up to date. Obtain the Allele sequences file and MLST profiles file, e.g. https://pubmlst.org/bigsdb?db=pubmlst_bhenselae_seqdef

The output is (1) a .mlst file containing the individual MLST calls (tab-delimited, view in EXCEL), (2) a log file giving the frequency of each MLST, (3) a file for the relabel.pl script, and (4) a relaxed PHYLIP alignment file suitable for RAxML phylogenetic analysis.

*pubmlst* can only detect known alleles in the Allele sequences file, it cannot call new allele sequences (yet). However, it can call novel MLSTs (i.e. new combinations of known alleles).