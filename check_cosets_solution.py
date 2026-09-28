import numpy as np
import itertools
from functools import partial
import matplotlib.pyplot as plt
from matplotlib.patches import RegularPolygon
import sys
import multiprocessing as mp
from multiprocessing import freeze_support
sys.setrecursionlimit(1000) 
import sympy as sp

# --- load the helper file with the function definitions --- #
from src.iteration_functions import *


n = 6
shape_cycles = ["CCC","T"]
determine_everything(n,shape_cycles,"")
check_sympy_solution(n)


# check if the produced index lists are the same
with open('S{}/associated_permutation.json'.format(n), 'r') as file:
    sol_s6_normal = json.load(file)

with open('S{}/indices_sympy.json'.format(n), 'r') as file:
    sol_s6_sympy = json.load(file)


true_array = [False]*len(sol_s6_normal['permutations'])
for i in range(0,len(sol_s6_normal['permutations'])):
    if sorted(sol_s6_normal['permutations'][i]) == sorted(sol_s6_sympy['permutations'][i]):
        true_array[i] = True
    else:
        print("No match for index ",i)
        print("Normal: ",sol_s6_normal['permutations'][i])
        print("Sympy:  ",sol_s6_sympy['permutations'][i]) 

true_check = all(x == True for x in true_array)
print("All permutations match: ",true_check)



with open('S{:.0f}/ngons.json'.format(n), 'r') as file:
    ngons_normal = json.load(file)

with open('S{:.0f}/ngons_sympy.json'.format(n), 'r') as file:
    ngons_sympy = json.load(file)

true_array = [False]*len(ngons_normal['ngons'])
for i in range(0,len(ngons_normal['ngons'])):
    if sorted(ngons_normal['ngons'][i]) == sorted(ngons_sympy['ngons'][i]):
        true_array[i] = True
    else:
        print("No match for index ",i)
        print("Normal: ",ngons_normal['ngons'][i])
        print("Sympy:  ",ngons_sympy['ngons'][i]) 

true_check = all(x == True for x in true_array)
print("All ngons match: ",true_check)
