# Explicit $S_n$ representations
This repository contains programs to automatically construct the 
explicit form for the multiplets of $S_6$. Some parts of it can also 
be generalized to generic $S_n$. <br>
***General workflow*** <br>
Before any Python scripts can be run, it is advisable to run `python setup_venv.py` 
in order to set up a virtual environment with all the package 
requirements. This step, however, is optional. After that, one 
continues with the following workflow.
1. Determine which $S_n$ you want to consider.
2. Run `python3 Sn_matrices.py` with with the chosen $n$ to automatically generates 
and saves 
the transposition and cyclic permutation representation matrices for every 
multiplet of the chosen $S_n$. The default here is set to $n=6$, however, this 
runs for any $n$ (given enough computing power). 
3. Once these are computed, one can determine the Plaquettes of elements 
under the dihedral subgroup $D_4$ by running `python3 determine_cosets_Sn.py`. 
This, however, currently only works for even $n$, since odd $n$ do not 
display the $D_4$ subgroup structure.
4. 

## Explicit transformation matrices for $S_n$ irreducible representations
Algorithmic determination of 
the explicit representation matrices for transpositions 
$(ij)$ and cyclic permutations $(123... n)$ for each 
multiplet $\mathcal{M}$ of $S_6$ using the 
`Sn_matrices.py` script. <br> 
All the functions needed for this are collected in 
`src/young_functions.py`. <br>
The default is set to $n=6$, but it works for arbitrary $n$. <br>
The resulting matrices are stored as `matrix_mu_ij.txt`in the 
folders `S6/matrices/mu` for a given multiplet $\mathcal{M}$. 
The `mu` denotes the associated partition of a multiplet $\mathcal{M}$. 

## Determine the cosets of $S_n$ elements
The program `determine_cosets_Sn.py` automatically determines the cosets of 
$S_n$ with respect to the dihedral group $D_4$ using `SymPy`. <br>
Also here, the default is set to $n=6$, but the program works for 
arbitrary $n$ with $n$ even. The odd $n$ do not display the coset structure 
under $D_4$. <br>
To check if the determined cosets are an exact cover for the permutation 
group $S_n$, one can use the program `check_cosets_solution.py`. This uses 
the Dancing Links (DLX) algorithm of Donald Knuth [see here](https://arxiv.org/abs/cs/0011047) 
to determine the solution to the exact cover problem.

## Explicitly construct the representations of the multiplets
The Mathematica file `Multiplet_construction.nb` collects all the 
cosets into plaquette variables $\pi$ 
and constructs the explicit representations for the $S_6$ 
(anti)multiplets as outlined in the corresponding publication. 