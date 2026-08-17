# Explicit $S_n$ representations
## Explicit transformation matrices for $S_n$ irreducible representations
Algorithmic determination of 
the explicit matrices for the irreducible representations of $S_6$ using the 
`Sn_matrices.py` script. <br> 
All the functions needed for this are collected in 
`young_functions.py`. <br>
The default is set to $n=6$, but it works for arbitrary $n$.

## Determine the cosets of $S_n$ elements
The program `cosets_Sn.py` automatically determines the cosets of 
$S_n$ with respect to the dihedral group $D_4$ using `SymPy`. <br>
Also here, the default is set to $n=6$, but the program works for 
arbitrary $n$ with $n$ even.

## Explicitly construct the representations of the multiplets
The Mathematica file `explicit_multiplet_construction.nb` collects all the 
cosets and constructs the explicit representations for the $S_6$ 
(anti)multiplets.
