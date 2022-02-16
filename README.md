# SARS-CoV-2_Gene_Finder
Gene Prediction Process for the SARS-CoV-2 Virus 

## Introduction
Computational approaches for analyzing biological data have revolutionized many subfields of biology.
This impact has been so large that a new term has been coined to describe this type of interdisciplinary fusion of biology, 
computation, math, and statistics: bioinformatics.

In this project, I built a gene predictor of the SARS-CoV-2 virus using Python to identify potential protein-coding genes in the SARS-CoV-2 virus, 
and use the [protein-BLAST search engine](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastp&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome)
to determine which of these genes may encode key proteins of the virus, gaining insight of the process that goes into understanding COVID-19.

Within the `sars-cov-2.ipynb` notebook is the function implementation used to extract the genes in the virus and completed unit tests on the various functions in
`test_gene_finder.py`

## Credits 

This project was developed within my Software Design class during my first year at the [Olin College of Engineering](https://github.com/olin).
