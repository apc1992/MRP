# MRP

<img src="img/logo.png" width="400">

Match a series of peptides into a custom reference proteome and filter the identical matches

# Description and workflow

This tool is aimed to filter out any peptide that aligns perfectly to any entry in a custom protein reference. It takes as input a protein fasta file and generates all possible derived n-mer peptides for each protein through a slidding window scan with a pre-defined minimum and maximum length. It can directly take as input the peptides list in a text file with a column named 'peptides'. 

The custom protein reference must be in fasta format, which will be indexed prior to aligning. Each peptide will be aligned into the reference and the perfect matches will be dropped from the list. 

The tool was originally designed for filtering MHC Class I peptides that match the human proteome, in order to prioritize foreign peptides for immunogenicity assays. 

<img src="img/flow_diagram.png" width="800"> 

# Input

T
