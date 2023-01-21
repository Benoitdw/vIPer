# vIPer

Here there is the new python code relative to the paper "Estimating the vertical ionization potential of single-stranded DNA molecules" (M.Rooman and F. Pucci, submitted).

The usage is very simple, just clone the repository and run python vip.py, input the sequence of the nucleotide stack of any length and in output you will receive it vertical ionization potential. 

For single bases, doublets, triplets and quadruplets the vIP are calculated via MP2 quantum chemistry approach. Fore longer stacks, we predict the vIP value using an iterative model that is a function of the 4plet values (for more details see the manuscript).  
