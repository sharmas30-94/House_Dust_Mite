## House Dust Mite detection from unmapped Whole Genome Sequencing reads
This repository documents a pipeline to detect house dust mite (HDM) sequences-specifically Dermatophagoides farinae and Dermatophagoides pteronyssinus from unmapped human WGS reads.

## Overview
We aligned 1,217 WGS libraries to the human T2T-CHM13 v2.0 reference and extracted the unmapped reads. These unmapped reads were then aligned against HDM reference assemblies:
- Dermatophagoides farinae (GCA_020809275.1)
- Dermatophagoides pteronyssinus (GCA_003076615.3)
From each alignment, we retained read alignments that passed stringent mapping and length filters and then validated species identity by BLAST against the NCBI nt database with low-complexity masking enabled.
