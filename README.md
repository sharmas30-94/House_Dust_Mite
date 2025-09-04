## House Dust Mite detection from unmapped Whole Genome Sequencing reads
This repository documents a pipeline to detect house dust mite (HDM) sequences-specifically Dermatophagoides farinae and Dermatophagoides pteronyssinus from unmapped human WGS reads.

## Overview
We aligned 1,217 WGS libraries to the human T2T-CHM13 v2.0 reference and extracted the unmapped reads. These unmapped reads were then aligned against HDM reference assemblies:

- Dermatophagoides farinae (GCA_020809275.1)
- Dermatophagoides pteronyssinus (GCA_003076615.3)

From each alignment, we retained read alignments that passed stringent mapping and length filters and then validated species identity by BLAST against the NCBI nt database with low-complexity masking enabled.

## Running Blast on Biowulf Server
Easyblast is an easy interface to Blast on Biowulf. It is a wrapper script that will prompt you for all required parameters, set up your jobs appropriately and submit them to the cluster. You will need to have all your query sequences in multiple files in a single directory (multiple sequences per file is fine). Easyblast will run the latest version of Blast+. The version will be printed at the beginning of the Easyblast run, and will also appear in the summary file in the output directory, and the actual Blast output files.
