# vIPer


[TOC1.pdf](https://github.com/3BioCompBio/vIPer/files/10472583/TOC1.pdf)




Here there is the new python code relative to the paper "Calculating vertical ionization potentials of DNA stacks" (M.Rooman and F. Pucci, submitted).

The usage is very simple, just clone the repository and run python vip.py, input the sequence of the nucleotide stack of any length and in output you will receive it vertical ionization potential. 

For single bases, doublets, triplets and quadruplets the vIP are calculated via MP2 quantum chemistry approach. Fore longer stacks, we predict the vIP value using an iterative model that is a function of the 4plet values (for more details see the manuscript).  
