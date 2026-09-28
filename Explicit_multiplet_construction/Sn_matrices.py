# --------------------------------------------------------------------------- #
#                     Necessary libraries                                     #
# --------------------------------------------------------------------------- #
import numpy as np
import math
from numpy.linalg import multi_dot
import sympy as sp
import os

# --------------------------------------------------------------------------- #
#                    Load helper functions                                    #
# --------------------------------------------------------------------------- #
from young_functions import *       

# --------------------------------------------------------------------------- #
#               Create folders to store resulting matrices                    #
# --------------------------------------------------------------------------- #
def create_folder(dir_name):
    '''
        This function automatically creates all necessary subfolders to 
        store the explicit matrices for a given irreducible 
        representation/partition
    '''
    try:
        os.mkdir(dir_name)
        print(f"Directory '{dir_name}' created successfully.")
    except FileExistsError:
        print(f"Directory '{dir_name}' already exists.")
    except PermissionError:
        print(f"Permission denied: Unable to create '{dir_name}'.")
    except Exception as e:
        print(f"An error occurred: {e}")


# --------------------------------------------------------------------------- #
#           Main function to construct the representation matrices            #
# --------------------------------------------------------------------------- #
def determine_everything(n):
    '''
        This function uses all helper functions to determine all representation 
        matrices for a given permutation group S_n.
    '''

    # --- create folders to store the found matrices --- #
    path = "S"+str(n)
    create_folder(path)
    pathmat = path+"/matrices"
    create_folder(pathmat)

    # --- Start with simplest partition: n times the 1 --- #
    comb = list(np.ones(n))
    comblist = [comb]

    # --- Determine all possible partitions --- #
    N = 2
    ic = 0
    while True:
        determine_combination(comblist[ic],comblist,N,n)
        if comblist[-1][0] == n:
            print("All combinations found!")
            break
        ic += 1

    # --- Convert partitions from string to integers and sort the list --- #
    for i in range(0,len(comblist)):
        for j in range(0,len(comblist[i])):
            comblist[i][j] = int(comblist[i][j])

    comblist.sort()

    # --- Create subfolders for the partitions to store the matrices --- #
    for i in range(0,len(comblist)):
        name = ""
        for j in range(0,len(comblist[i])):
            name += str(comblist[i][j])
            
        create_folder(pathmat+"/"+name)

    # --- Determine sizes of conjugacy classes --- #
    c_size = np.zeros((len(comblist)))
    for i in range(0,len(c_size)):
        c_size[i] = compute_c_size(comblist[i],n)
    print(c_size)

    # --- Determine the Young tableaus for all the partitions and get their dimension --- #
    Sn_matrices = []
    for i in range(0,len(comblist)):
        Sn_matrices.append(young_tableaux_from_dim(comblist[i],n))
    Sn_dims = []
    for i in range(0,len(Sn_matrices)):
        dim = calc_dim(Sn_matrices[i])
        Sn_dims.append(int(dim))
    print("Dims: ",Sn_dims)

    # --- Check if the total dimension of the group S_n is recovered --- #
    dimcheck = 0.0
    for i in Sn_dims:
        dimcheck += i**2
    print("Dimension check: ", dimcheck)

    # --- Determine the matrices for all the partitions -- #
    for i in range(0,len(comblist)):
        determine_matrix(comblist[i],comblist,Sn_dims,Sn_matrices)
        
    
# --------------------------------------------------------------------------- #
#          Construct all representation matrices for a given S_n              #
# --------------------------------------------------------------------------- #
n = 6
determine_everything(n)