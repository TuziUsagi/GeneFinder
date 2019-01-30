# Readme

## Required package

### numpy
### biopython
### pandas

## File description

### exonChainMain.py
This file implement the exon chaining method for finding a gene in a unannotated genome. The putative exons are found with sequence alignment. 

### exonChainFunctions.py
This file contains all helper functions for the exonChainMain.py. THe align_local function runs Smith-Waterman algorithm for finding a region of the genome containing multiple putative exons. 

### exonChainRun.py
This is an example for running the exon chaining algorithm. 

## Results
The putative exons found by local alignment:
![putative exons](https://raw.githubusercontent.com/TuziUsagi/PyPSC/master/exampleFigures/PSC_example.png)

The gene finding result:
![Gene](https://raw.githubusercontent.com/TuziUsagi/PyPSC/master/exampleFigures/baseline_example.png)
