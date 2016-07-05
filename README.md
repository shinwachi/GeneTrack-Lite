# GeneTrack-Lite
Smoothing and peak calling module used in early version of GeneTrack, with no dependency on other modules for GeneTrack.

## Why?
Created mainly to avoid going through the gene-track web interface
by allowing command-line scripting of peak calling.
## How?
Ripped out from old genetrack code (MIT Open Source License), this module produces the 
"gaussian kernel" smoothed values and picks peaks, just like the web-interface counterpart.

## Citation
For citation, please refer to:
Albert I, Wachi S, Jiang C, Pugh BF, Genetrack--a genomic data processing and visualization framework.
Bioinformatics 2008 May 15;24(10):1305-6
PMID: 18388141
