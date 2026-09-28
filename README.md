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
cosets into plaquette variables $\pi$ 
and constructs the explicit representations for the $S_6$ 
(anti)multiplets. <br>
In the folder `Sn/S6/matrices/` are the collected matrices $a_i$ needed 
for the explicit construction of a given multiplet $M_j = \sum_i a_i X_i$, where 
$X_i$ are a combination of the plaquettes $\pi$. 
The permutations to be applied to the ordering of the 
plaquettes $\pi$ on the Chessboard Planes (CPs) 
for a specific multiplet $M_j$ are collected in `Change_M_j.mx` 
and the corresponding powers of the 
rotation matrix in `R_j.mx`. <br> 
For the 16-plet, one additionally needs 
transformations applied to the plaquettes, which are collected in 
`replacements_1_2_lc.mx`, `replacements_5_6_lc.mx`, `replacements_10_11_lc.mx`, 
`replacements_13_14_lc.mx`. 
