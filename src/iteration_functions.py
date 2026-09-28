import itertools
from itertools import permutations
from functools import partial
import sys
from src.helper_functions import *
import multiprocessing
from sympy.combinatorics import Permutation
from multiprocessing import Pool, cpu_count, Process, Manager
import json
import time
import numpy as np

# --------------------------------------------------------------------------- #
#                              Timing utility                                 #
# --------------------------------------------------------------------------- #

# This function records the current wall-clock time and process CPU time.
def get_time():
    '''
        This function records the current wall-clock time and process CPU time.
    '''
    t = time.time()
    t_cpu = time.process_time()
    return t, t_cpu

# This function calculates the elapsed wall-clock time and process CPU time.
def calc_time_taken(_t, _t_cpu):
    '''
        This function calculates the elapsed wall-clock time and process CPU time.
    '''
    t_a = time.time()
    t_a_cpu = time.process_time()

    dt = t_a - _t
    dt_cpu = t_a_cpu - _t_cpu
    return dt, dt_cpu

sys.setrecursionlimit(1000) 

# --------------------------------------------------------------------------- #
#                          File-writing utilities                             #
# --------------------------------------------------------------------------- #

# This function writes a collection of combinations to a text file.
def write_file(combs, name):
    '''
        This function writes a collection of combinations to a text file.
    '''

    # --- Open the output file for writing. --- #
    hexfile = open(name,"w")
    for i in range(0,len(combs)):
        hexfile.write("[")
        for j in range(0,len(combs[i])):
            if j < len(combs[i])-1:
                hexfile.write("{:.0f}, ".format(combs[i][j]))
            else:
                hexfile.write("{:.0f}".format(combs[i][j])) 
        hexfile.write("]\n")
    hexfile.close()

# --------------------------------------------------------------------------- #
#                Permutation and cyclic-operation functions                   #
# --------------------------------------------------------------------------- #

# This function cyclically reorders a sequence so that the element at the given index comes first.
def reorder_from_idx(idx, a):
    '''
        This function cyclically reorders a sequence so that the element at the given index comes first.
    '''
    return a[idx:] + a[:idx]

# This function creates partial functions representing all cyclic reorderings of a sequence.
def cyclic_perm(a):
    '''
        This function creates partial functions representing all cyclic reorderings of a sequence.
    '''
    return [partial(reorder_from_idx, i) for i in range(len(a))]

# This function swaps the elements 1 and 2 in a permutation.
def P12(a):
    '''
        This function swaps the elements 1 and 2 in a permutation.
    '''
    at = list(a[:])
    ind1 = at.index(1)
    ind2 = at.index(2)

    at[ind1],at[ind2] = at[ind2],at[ind1]
    at = tuple(at)
    return at


# This function swaps two specified elements i and j in a permutation.
def Pij(a,i,j):
    '''
        This function swaps two specified elements i and j in a permutation.
    '''
    at = list(a[:])
    ind1 = at.index(i)
    ind2 = at.index(j)
    at[ind1],at[ind2] = at[ind2],at[ind1]
    at = tuple(at)
    return at

# This function applies the fixed sequence of swaps used for the six-element cyclic operation.
def P123456(a):
    '''
        This function applies the fixed sequence of swaps used for the six-element cyclic operation.
    '''
    aorg = a
    for i in range(0,5):
        aorg = Pij(aorg,5-i,6-i)
    return aorg

# This function applies the cyclic permutation operation to a sequence of arbitrary length.
def Pcycl(a):
    '''
        This function applies the cyclic permutation operation to a sequence of arbitrary length.
    '''
    aorg = a
    dima = len(a)
    dimam1 = len(a)-1
    for i in range(0,dimam1):
        aorg = Pij(aorg,dimam1-i,dima-i)
    return aorg

# --------------------------------------------------------------------------- #
#                       Path-determination functions                          #
# --------------------------------------------------------------------------- #

# This function determines a path using the six-element cyclic operation together with the 1-2 swap.
def determine_path(perm,cycle,permlist):
    '''
        This function determines a path using the six-element cyclic operation together with the 1-2 swap.
    '''
    perms_of_met_hex = []

    pnew = perm
    counter = 0
    ellist = []

    # --- Apply the requested number of cyclic operations before each transposition. --- #
    while True:

    # --- Apply the requested number of cyclic operations before each transposition. --- #
        for i in range(0,cycle):
            pnew = P123456(pnew)

    # --- Store the permutation indices encountered along the path. --- #
        ellist.append(permlist.index(pnew))
        pnew = P12(pnew)

    # --- Store the permutation indices encountered along the path. --- #
        ellist.append(permlist.index(pnew))

    # --- Generate all cyclic positions associated with the current path element. --- #
        perms_of_met_hex.append([pnew])
        pinter = pnew[:]
        for i in range(0,5):
            pinter = P123456(pinter)
            perms_of_met_hex[counter].append(pinter)

        counter += 1
        if pnew == perm:
            break
    
    return perms_of_met_hex, ellist

# This function determines a path using the general cyclic operation together with the 1-2 swap.
def determine_general_path(perm,cycle,permlist):
    '''
        This function determines a path using the general cyclic operation together with the 1-2 swap.
    '''
    perms_of_met_shape = []

    pnew = perm
    counter = 0
    ellist = []

    # --- Apply the requested number of general cyclic operations before each transposition. --- #
    while True:

    # --- Apply the requested number of general cyclic operations before each transposition. --- #
        for i in range(0,cycle):
            pnew = Pcycl(pnew)

    # --- Store the permutation indices encountered along the path. --- #
        ellist.append(permlist.index(pnew))
        pnew = P12(pnew)

    # --- Store the permutation indices encountered along the path. --- #
        ellist.append(permlist.index(pnew))

    # --- Generate all cyclic positions associated with the current path element. --- #
        perms_of_met_shape.append([pnew])
        pinter = pnew[:]
        for i in range(0,len(pinter)):
            pinter = Pcycl(pinter)
            perms_of_met_shape[counter].append(pinter)

        counter += 1
        #test if we are back at the original element
        if pnew == perm:
            break
    print("Took ",counter," steps to get back to the same element")
    
    return perms_of_met_shape, ellist


# This function determines a path from an arbitrary sequence of cyclic (C) and transposition (T) operations.
def determine_arbitrary_path(perm,cycle,_permlist):
    '''
        This function determines a path from an arbitrary sequence of cyclic (C) and transposition (T) operations.
    '''
    perms_of_met_shape = []

    pnew = perm
    counter = 0
    ellist = []

    # --- Apply each operation in the supplied path description. --- #
    while True:

    # --- Apply each operation in the supplied path description. --- #
        for i in range(0,len(cycle)):

    # --- Apply each operation in the supplied path description. --- #
            for k in range(0,len(cycle[i])):
                if cycle[i][k] == "C":
                    pnew = Pcycl(pnew)
                elif cycle[i][k] == "T":
                    pnew = P12(pnew)

    # --- Store the permutation index reached after each group of operations. --- #
            ellist.append(_permlist.index(pnew))

    # --- Generate the cyclic positions associated with the current path element. --- #
        perms_of_met_shape.append([pnew])
        pinter = pnew[:]
        for i in range(0,len(pinter)):
            pinter = Pcycl(pinter)
            perms_of_met_shape[counter].append(pinter)

        counter += 1
        #test if we are back at the original element
        if pnew == perm:
            break
    
    return perms_of_met_shape, ellist, counter

# --------------------------------------------------------------------------- #
#                       Shape-determination functions                         #
# --------------------------------------------------------------------------- #

# This function determines the cyclic shape containing a permutation and adds it to the shape list if it is new.
def find_shape(perm,cycle,shape):
    '''
        This function determines the cyclic shape containing a permutation and adds it to the shape list if it is new.
    '''
    nodes = [perm]
    testperm = perm
    for i in range(0,cycle-1):
        testperm = Pcycl(testperm)
        nodes.append(testperm)
    
    nodes.sort()
    if nodes not in shape:
        shape.append(nodes)

# This function determines and returns the sorted cyclic shape containing a permutation.
def find_shape_parallel(perm,cycle):
    '''
        This function determines and returns the sorted cyclic shape containing a permutation.
    '''
    nodes = [perm]
    testperm = perm
    for i in range(0,cycle-1):
        testperm = Pcycl(testperm)
        nodes.append(testperm)
    
    nodes.sort()
    return nodes

# This function determines all distinct shapes represented by the supplied permutation list.
def determine_shape_list(_perm_list,_n):
    '''
        This function determines all distinct shapes represented by the supplied permutation list.
    '''
    shape = []

    for i in range(0,len(_perm_list)):
        find_shape(_perm_list[i],_n,shape)
    print("Found {} shapes with {} corners from the {} permutations!".format(len(shape),_n,len(_perm_list)))
    return shape

# This function determines a hexagon from a permutation and adds it to the supplied list if it is new.
def find_hexagons(perm,cycle,hexagons):
    '''
        This function determines a hexagon from a permutation and adds it to the supplied list if it is new.
    '''
    nodes = [perm]
    testperm = perm
    for i in range(0,cycle-1):
        testperm = P123456(testperm)
        nodes.append(testperm)

    nodes.sort()
    if nodes in hexagons:
        print("This hexagon has already been added!")
    else:
        hexagons.append(nodes)


# This function determines and returns the six-element cyclic structure associated with a permutation.
def determine_single_hex(perm,cycle):
    '''
        This function determines and returns the six-element cyclic structure associated with a permutation.
    '''
    nodes = [perm]
    testperm = perm
    for i in range(0,cycle-1):
        testperm = P123456(testperm)
        nodes.append(testperm)
    return nodes

# This function finds the index of the hexagon containing a given permutation.
def get_hex(perm,hexagons):
    '''
        This function finds the index of the hexagon containing a given permutation.
    '''
    for i in range(0,len(hexagons)):
        if perm in hexagons[i]:
            return i
        else:
            continue


# --------------------------------------------------------------------------- #
#                     Permutation-generation functions                        #
# --------------------------------------------------------------------------- #

# This function generates all permutations of the integers from 1 through n_dim.
def generate_Sn_permutations(n_dim):
    '''
        This function generates all permutations of the integers from 1 through n_dim.
    '''

    # --- Create the list of values whose permutations are required. --- #
    sn_list = []
    for i in range(1,n_dim+1):
        sn_list.append(i)

    # --- Generate all permutations of the values. --- #
    sn_perm_list = list(itertools.permutations(sn_list))

    del sn_list

    return sn_perm_list


# --------------------------------------------------------------------------- #
#                      Shape and path index functions                         #
# --------------------------------------------------------------------------- #

# This function determines the unique index lists associated with the corners of the supplied shapes.
def determine_shape_corner_index_list(shapelist, cycl, i_permlist):
    '''
        This function determines the unique index lists associated with the corners of the supplied shapes.
    '''
    i_compl_indlist_sorted = []
    i_compl_indlist = []
    for j in range(0,len(shapelist)):
        for i in range(0,len(shapelist[0])):
            perms_of_met_hex, indlist_ = determine_general_path(shapelist[j][i],cycl,i_permlist)
            perms_of_met_hex, indlist_sort = determine_general_path(shapelist[j][i],cycl,i_permlist)
            indlist_sort.sort()
            if indlist_sort in i_compl_indlist_sorted:
                continue
            else:
                i_compl_indlist_sorted.append(indlist_sort)
                i_compl_indlist.append(indlist_)
    return i_compl_indlist

# This function checks whether an index from one list occurs in any of the supplied sorted index lists.
def test_indices(_indlist,_sorted_indlist):
    '''
        This function checks whether an index from one list occurs in any of the supplied sorted index lists.
    '''
    flag = False
    for i in _indlist:
        count = 0
        for j in range(0,len(_sorted_indlist)):
            if i in _sorted_indlist[j]:
                flag = True
                return flag


# --------------------------------------------------------------------------- #
#                           Exact-cover functions                             #
# --------------------------------------------------------------------------- #

# This function searches recursively for an exact cover of the supplied universe.
def find_exact_cover(universe, subsets):
    '''
        This function searches recursively for an exact cover of the supplied universe.
    '''
    if not universe:
        return []

    # get next number in range 0-n
    element = next(iter(universe))

    # build set of all arrays containing the number element
    options = [s for s in subsets if element in s]

    # backtracking search algorithm
    for option in options:
        new_universe = universe - set(option)
        new_subsets = [s for s in subsets if not (set(s) & set(option))]

        # recursive calling of the algorithm
        result = find_exact_cover(new_universe, new_subsets)

        if result is not None:
            return [option] + result

    return None    


class Node:

    def __init__(self, col=None, row_data=None):
        self.left = self
        self.right = self
        self.up = self
        self.down = self
        self.col = col
        self.row_data = row_data

class Header(Node):

    def __init__(self, name):
        super().__init__(self)
        self.size = 0
        self.name = name


# --------------------------------------------------------------------------- #
#                       Dancing-links (DLX) functions                         #
# --------------------------------------------------------------------------- #

# This function constructs the dancing-links matrix used for the exact-cover problem.
def get_dlx_matrix(universe_elements, subsets):
    '''
        This function constructs the dancing-links matrix used for the exact-cover problem.
    '''

    # --- Create one column header for every element in the universe. --- #
    universe = {j: Header(j) for j in universe_elements}

    # --- Create the root of the circular column-header list. --- #
    root = Header("root")
    
    root.right = root
    root.left = root
    
    # --- Link all column headers into the root's circular doubly linked list. --- #
    for j in universe_elements:
        j_header = universe[j]
        j_header.left = root.left
        j_header.right = root
        root.left.right = j_header
        root.left = j_header

    # --- Add one row of nodes for every subset. --- #
    for subset in subsets:
        row_nodes = []
        for element in subset:
            if element not in universe:
                continue
            
            col_header = universe[element]
            node = Node(col_header, row_data=subset)
            row_nodes.append(node)
            
    # --- Link nodes vertically into their corresponding columns. --- #
            node.up = col_header.up
            node.down = col_header
            col_header.up.down = node
            col_header.up = node
            col_header.size += 1

        if row_nodes:

    # --- Link nodes horizontally within each subset. --- #
            for i in range(len(row_nodes) - 1):
                current_node = row_nodes[i]
                next_node = row_nodes[i + 1]
                current_node.right = next_node
                next_node.left = current_node
            
            row_nodes[0].left = row_nodes[-1]
            row_nodes[-1].right = row_nodes[0]
            
    return root


# This function covers a column and removes the associated conflicting nodes.
def cover(c):
    '''
        This function covers a column and removes the associated conflicting nodes.
    '''
    c.right.left = c.left
    c.left.right = c.right
    
    i = c.down
    while i != c:
        j = i.right
        while j != i:
            j.down.up = j.up
            j.up.down = j.down
            j.col.size -= 1
            j = j.right
        i = i.down

# This function restores a previously covered column and its associated nodes.
def uncover(c):
    '''
        This function restores a previously covered column and its associated nodes.
    '''
    i = c.up
    while i != c:
        j = i.left
        while j != i:
            j.col.size += 1
            j.down.up = j
            j.up.down = j
            j = j.left
        i = i.up
    
    c.right.left = c
    c.left.right = c


# This function solves the exact-cover problem using an iterative implementation of the DLX algorithm.
def solve_dlx_iterative(root):
  """
  An iterative implementation of the DLX algorithm.
  This version uses an explicit stack to manage state and avoids RecursionError.
  """

    # --- Store selected solution rows and the explicit search stack. --- #
  solution_nodes = []
  stack = []
  
  c = root.right # Start with the first column

  loop_counter = 0

    # --- Continue searching until either an exact cover is found or all possibilities are exhausted. --- #
  while True:

    # --- Continue searching until either an exact cover is found or all possibilities are exhausted. --- #
      # Step 1: Check for solution
      if c == root:
          return [node.row_data for node in solution_nodes]
      
    # --- Choose the column with the fewest remaining nodes. --- #
      # Step 2: Choose column with MRV heuristic
      min_size = sys.maxsize
      best_c = None
      current_c = root.right
      while current_c != root:
          if current_c.size < min_size:
              min_size = current_c.size
              best_c = current_c
          current_c = current_c.right
      c = best_c

    # --- If the chosen column is empty, backtrack to the previous choice. --- #
      # Step 3: Handle a column with no available rows (backtrack)
      if c.down == c:
          if not stack:
            elements_recovered = [node.row_data for node in solution_nodes]
            print(f"Recovered only {len(elements_recovered)} elements. Excited with {loop_counter} counts.")
            return None # No solution
          c = stack.pop()
          uncover(c)
          continue
      
    # --- Cover the chosen column and select its first available row. --- #
      # Step 4: Cover column and try the first row
      cover(c)
      r = c.down
      solution_nodes.append(r)
      stack.append(c)
      
    # --- Cover all other columns represented by the selected row. --- #
      # Step 5: Cover other columns in the chosen row
      j = r.right
      while j != r:
          cover(j.col)
          j = j.right
      
      c = root.right # Continue search with the first column
      loop_counter += 1
      
  return None  

# --------------------------------------------------------------------------- #
#                   Arbitrary-path exact-cover functions                      #
# --------------------------------------------------------------------------- #

# This function determines path index sets for all permutations and searches for an exact cover.
def determine_shape_corner_index_list_arb_path(shapelist, cycl, i_permlist):
    '''
        This function determines path index sets for all permutations and searches for an exact cover.
    '''

    sol_flag = False
    i_compl_permlist = []

    print("Determining path for all permutations.")

    for ii in range(0,len(i_permlist)):
        perms_of_met_hex, _indlist_, _counter = determine_arbitrary_path(i_permlist[ii],cycl,i_permlist)
        i_compl_permlist.append(_indlist_)
    

    print("Found all {} paths. Takes {} steps with {} indices to get back".format(len(i_compl_permlist),_counter,len(_indlist_)) )

    # --- The path index sets form the subsets of the exact-cover problem. --- #
    universe = set(range(len(i_permlist)))

    all_subsets = i_compl_permlist


    # --- The universe contains one element for every permutation. --- #
    print("Finding solutions")
    dlx_matrix_time = time.process_time()

    # --- Construct the DLX matrix. --- #
    matrix_root = get_dlx_matrix(universe, all_subsets)
    print(f"Time taken for get_dlx_matrix {(time.process_time() - dlx_matrix_time):.4f} s")

    dlx_solver_time = time.process_time()

    # --- Solve the exact-cover problem. --- #
    solution = solve_dlx_iterative(matrix_root)
    print(f"Time taken for iterative solver {(time.process_time() - dlx_solver_time):.4f} s")
    

    CRED = '\033[91m'
    CEND = '\033[0m'

    # solution = find_exact_cover(universe,all_subsets)
    if solution:
        print(f"\033[92mFound a solution with\033[0m {len(solution)} \033[92mlists.\033[0m")
        sol_flag = True

    else:
        print(CRED+"No solution found. -> Please provide a different path!"+CEND)
        raise Exception(CRED+"No solution found. -> Please provide a different path!"+CEND)


    return solution, sol_flag

# --------------------------------------------------------------------------- #
#                         SymPy solution functions                            #
# --------------------------------------------------------------------------- #

# This function reads a saved SymPy solution and converts its permutations to indices in the supplied permutation list.
def read_the_sympy_solution(_n,_permlist):
    '''
        This function reads a saved SymPy solution and converts its permutations to indices in the supplied permutation list.
    '''
    print("Reading the solution... ",end="")
    time_read = time.process_time()
    with open('sol_s{:.0f}.json'.format(_n), 'r') as file:
        data = json.load(file)

    tuple_to_index = {tuple(val): idx for idx, val in enumerate(_permlist)}
    
    # Now just look up in O(1)
    all_subsets_net = []
    for subset in data['solution']:
        inter_subsets = [tuple_to_index[tuple(item)] for item in subset]
        all_subsets_net.append(inter_subsets)
    print(f"Time taken to read the solution {(time.process_time() - time_read):.4f} s")

    for i in range(0,len(data['solution'])):
        for j in range(0,len(data['solution'][i])):
            data['solution'][i][j] = tuple(data['solution'][i][j])

    return all_subsets_net, data['solution']


# This function uses the SymPy-derived subsets to construct and solve the corresponding exact-cover problem.
def determine_shape_corner_index_list_arb_path_sympy(_n, _all_subsets_net):
    '''
        This function uses the SymPy-derived subsets to construct and solve the corresponding exact-cover problem.
    '''
    sol_flag = False
    n_fac = math.factorial(_n)
    universe = set(range(n_fac))

    print("Finding solutions...")

    dlx_matrix_time = time.process_time()

    matrix_root = get_dlx_matrix(universe, _all_subsets_net)

    print(f"Time taken for get_dlx_matrix {(time.process_time() - dlx_matrix_time):.4f} s")

    dlx_solver_time = time.process_time()

    solution = solve_dlx_iterative(matrix_root)

    print(f"Time taken for dlx solver {(time.process_time() - dlx_solver_time):.4f} s")

    CRED = '\033[91m'
    CGREEN = '\033[92m'
    CEND = '\033[0m'

    if solution:
        print(f"\033[92mFound a solution with\033[0m {len(solution)} \033[92mlists.\033[0m")
        sol_flag = True

    else:
        raise Exception(CRED+"No solution found. -> Please provide a different path!"+CEND)
        


    return solution, sol_flag


# --------------------------------------------------------------------------- #
#                      Permutation subclass functions                       #
# --------------------------------------------------------------------------- #

# This function tests whether a permutation is related to a target permutation by the T operation.
def test_T_relation(a_list,b,_perml):
    '''
        This function tests whether a permutation is related to a target permutation by the T operation.
    '''
    perm_shapes = []
    for i in a_list:
        perm_shapes.append(_perml[i])
    for il in perm_shapes:
        if P12(il) == _perml[b]:
            return True

# This function separates associated permutations into subclasses according to their T relation.
def determine_subclasses(_assoc_perm,_permlist):
    '''
        This function separates associated permutations into subclasses according to their T relation.
    '''
    subclasses_perm = []
    subclasses_permT = []
    for i in range(0,len(_assoc_perm)):
        apl = [_assoc_perm[i][0]]
        aplT = []
        for k in range(1,len(_assoc_perm[i][:])):
            T_related = test_T_relation(apl,_assoc_perm[i][k],_permlist)
            if T_related == True:
                aplT.append(_assoc_perm[i][k])
            else:
                apl.append(_assoc_perm[i][k])

        apl_perms = []
        aplT_perms = []
        apl_perms_12 = []
        for j in range(0,len(apl)):
            apl_perms.append(_permlist[apl[j]])
            aplT_perms.append(_permlist[aplT[j]])
        subclasses_perm.append(apl_perms)
        subclasses_permT.append(aplT_perms)

    return subclasses_perm, subclasses_permT


# This function tests the T relation directly between permutation tuples.
def test_T_relation_sympy(a_list,b,):
    '''
        This function tests the T relation directly between permutation tuples.
    '''
    for il in a_list:
        if P12(il) == b:
            return True
        

# This function separates the permutations in a SymPy solution into T-related subclasses.
def determine_subclasses_sympy_solution(_data):
    '''
        This function separates the permutations in a SymPy solution into T-related subclasses.
    '''
    subclasses_perm = []
    subclasses_permT = []
    for i in range(0,len(_data)):
        apl = [_data[i][0]]
        aplT = []
        for k in range(1,len(_data[i][:])):
            T_related = test_T_relation_sympy(apl,_data[i][k])
            if T_related == True:
                aplT.append(_data[i][k])
            else:
                apl.append(_data[i][k])
        
        subclasses_perm.append(apl)
        subclasses_permT.append(aplT)


    return subclasses_perm, subclasses_permT


# --------------------------------------------------------------------------- #
#                  Shape association and output functions                     #
# --------------------------------------------------------------------------- #

# This function determines which shapes form the n-gons associated with the supplied permutation index groups.
def determine_ngons(_n,_shapes,_indlist,_permlist):
    '''
        This function determines which shapes form the n-gons associated with the supplied permutation index groups.
    '''
    ngon_of_shapes = []

    #precompute reverse lookup of permutations to shape
    perm_to_shape = {}
    for k, shape in enumerate(_shapes):
        for val in shape:
            perm_to_shape[val] = k


    pos_map = [{val: idx for idx, val in enumerate(shape)} for shape in _shapes]


    for ind_group in _indlist:
        inter_shape = []
        for ind in ind_group:

            inter_perm = _permlist[ind]
            k = perm_to_shape.get(inter_perm)
            if k is not None:
                inter_shape.append(k)
            # for k, shape in enumerate(_shapes):
            #     if inter_perm in pos_map[k]:
            #         inter_shape.append(k)

        inter_shape = list(dict.fromkeys(inter_shape))  # Remove duplicates while preserving order
        ngon_of_shapes.append(inter_shape)

    ngon_dict = {"ngons": ngon_of_shapes}
    with open("S{:.0f}/ngons_sympy.json".format(_n), 'w') as f:
        json.dump(ngon_dict, f)


# This function determines permutation positions and the shapes associated with each path.
def position_and_associated_shapes(shapes, _indlist, _permlist):
    '''
        This function determines permutation positions and the shapes associated with each path.
    '''
    pos_on_shape = []
    assoc_perm = []
    ngon_of_shapes = []
    

    pos_map = [{val: idx for idx, val in enumerate(shape)} for shape in shapes]


    for ind_group in _indlist:
        inter_pos = []
        inter_assoc_perm = []
        inter_shape = []
        for ind in ind_group:
            
            inter_perm = _permlist[ind]
            inter_assoc_perm.append(ind)
            for k, shape in enumerate(shapes):
                if inter_perm in pos_map[k]:
                    inter_shape.append(k)
                    inter_pos.append(pos_map[k][inter_perm])

        pos_on_shape.append(inter_pos)
        assoc_perm.append(inter_assoc_perm)
        inter_shape = list(dict.fromkeys(inter_shape))  # Remove duplicates while preserving order
        ngon_of_shapes.append(inter_shape)

    return assoc_perm, ngon_of_shapes


# --------------------------------------------------------------------------- #
#                          Hexagon path functions                             #
# --------------------------------------------------------------------------- #

# This function determines the unique path index lists associated with hexagon corners.
def determine_hexagon_corner_index_list(hexlist,i_permlist):
    '''
        This function determines the unique path index lists associated with hexagon corners.
    '''
    i_compl_indlist_sorted = []
    i_compl_indlist = []
    for j in range(0,len(hexlist)):
        for i in range(0,len(hexlist[0])):
            perms_of_met_hex, indlist_ = determine_path(hexlist[j][i],3,i_permlist)
            perms_of_met_hex, indlist_sort = determine_path(hexlist[j][i],3,i_permlist)
            indlist_sort.sort()
            
            if indlist_sort in i_compl_indlist_sorted:
                continue
            else:
                i_compl_indlist_sorted.append(indlist_sort)
                i_compl_indlist.append(indlist_)
                
    return i_compl_indlist

# --------------------------------------------------------------------------- #
#                             Output functions                                #
# --------------------------------------------------------------------------- #

# This function writes permutation groups to a text file in the required notation.
def write_perms(_permlist,filename):
    '''
        This function writes permutation groups to a text file in the required notation.
    '''
    file = open(filename,"w")
    file.write("# No. of vectors: {}\n".format(len(_permlist)))
    for i in range(0,len(_permlist)):
        file.write("")
        for j in range(0,int(len(_permlist[i]))):
            file.write("f_{")
            for k in range(0,len(_permlist[i][j])):
                if k < len(_permlist[i][j])-1:
                    file.write("{:.0f}".format(_permlist[i][j][k]))
                else:
                    file.write("{:.0f}".format(_permlist[i][j][k])) 
            if j < len(_permlist[i])-1:
                file.write("}\\ ")
            else:
                file.write("}")
        file.write("\n")
    file.close()

# This function writes permutation groups using alphabetic vector labels.
def write_vectors(_permlist,filename):
    '''
        This function writes permutation groups using alphabetic vector labels.
    '''
    alphabet = ["a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r"]
    file = open(filename,"w")
    for i in range(0,len(_permlist)):
        for j in range(0,int(len(_permlist[i]))):
            file.write("{")
            for k in range(0,len(_permlist[i][j])):
                if k < len(_permlist[i][j])-1:
                    file.write("{},".format(alphabet[_permlist[i][j][k]-1]))
                else:
                    file.write("{}".format(alphabet[_permlist[i][j][k]-1])) 
            if j < len(_permlist[i])-1:
                file.write("} ")
            else:
                file.write("}")
        file.write("\n")
    file.close()


# --------------------------------------------------------------------------- #
#                        High-level search functions                          #
# --------------------------------------------------------------------------- #

# This function generates the permutation list and the corresponding list of shapes.
def predetermine_perm_shape(_n):
    '''
        This function generates the permutation list and the corresponding list of shapes.
    '''
    _perm_list = generate_Sn_permutations(_n)

    _shape = determine_shape_list(_perm_list,_n)

    return _perm_list, _shape


# This function determines the remaining path solution for a supplied shape and cycle.
def determine_rest(_shape,_cycle,_perm_list):
    '''
        This function determines the remaining path solution for a supplied shape and cycle.
    '''
    _compl_indlist, _sol_flag = determine_shape_corner_index_list_arb_path(_shape,_cycle,_perm_list)
    if _sol_flag == False:
        return False

# This function tests all supplied combinations and collects the successful solutions.
def test_for_solution(_n,_combinations,_permlist_):
    '''
        This function tests all supplied combinations and collects the successful solutions.
    '''
    _sol_list = []
    print("{} combinations found".format(len(_combinations)))
    for i in range(0,len(_combinations)):
        obtained_solution = determine_rest(_n,_combinations[i],_permlist_)
        if obtained_solution == False:
            print("No solution found for path {}, continuing!".format(_combinations[i]))
            continue
        else:
            print("Solution found for path: ",_combinations[i])
            _sol_list.append(obtained_solution)
            continue

# --------------------------------------------------------------------------- #
#                   Path-combination generation functions                     #
# --------------------------------------------------------------------------- #

# This function generates and filters all candidate odd path combinations.
def generate_all_odd_combinations(_permlist_):
    '''
        This function generates and filters all candidate odd path combinations.
    '''
    all_odd_combinations = get_all_combinations_as_array()
    path_of_length_2 = []
    rest_of_the_paths = []
    for i_comb in all_odd_combinations:
        met_shapes, ellist, _counter = determine_arbitrary_path(_permlist_[0],i_comb,_permlist_)
        if _counter == 2:
            path_of_length_2.append(all_odd_combinations.index(i_comb))
        else:
            rest_of_the_paths.append(all_odd_combinations.index(i_comb))

    rest_of_the_paths.reverse()
    for i in rest_of_the_paths:
        all_odd_combinations.pop(i)

    return all_odd_combinations

# This function generates and filters all candidate even path combinations.
def generate_all_even_combinations(_permlist_):
    '''
        This function generates and filters all candidate even path combinations.
    '''
    all_even_combinations = get_all_even_combinations_as_array(4)
    path_of_length_4 = []
    rest_of_the_paths = []
    for i_comb in all_even_combinations:
        met_shapes, ellist, _counter = determine_arbitrary_path(_permlist_[0],i_comb,_permlist_)
        if _counter == 4:
            path_of_length_4.append(all_even_combinations.index(i_comb))
        else:
            rest_of_the_paths.append(all_even_combinations.index(i_comb))

    rest_of_the_paths.reverse()
    for i in rest_of_the_paths:
        all_even_combinations.pop(i)

    return all_even_combinations

# This function generates candidate paths of the requested type and tests them for solutions.
def brute_force_path(_n,_type):
    '''
        This function generates candidate paths of the requested type and tests them for solutions.
    '''
    _permlist, _shapes = predetermine_perm_shape(_n)

    
    if _type == "odd":
        _combs = generate_all_odd_combinations(_permlist)
    else:
        _combs = generate_all_even_combinations(_permlist)

    print("Paths to test:")
    for icomb in _combs:
        print(icomb)
    test_for_solution(_n,_combs,_permlist)    

# This function generates permutations and shapes and determines the corresponding path solution without writing output files.
def determine_everything_wo_write(_n, _cycle):
    '''
        This function generates permutations and shapes and determines the corresponding path solution without writing output files.
    '''
    _perm_list = generate_Sn_permutations(_n)

    _shape = determine_shape_list(_perm_list,_n)

    _compl_indlist, _sol_flag = determine_shape_corner_index_list_arb_path(_shape,_cycle,_perm_list)
    if _sol_flag == False:
        return False


# --------------------------------------------------------------------------- #
#                   Representative-based shape generation                     #
# --------------------------------------------------------------------------- #

# This function generates the orbit of a permutation under repeated cyclic operations.
def orbit_of(perm, n):
    '''
        This function generates the orbit of a permutation under repeated cyclic operations.
    '''
    nodes = []
    cur = perm
    for _ in range(n):
        nodes.append(cur)
        cur = Pcycl(cur)
    return nodes  # immutable representation

# This function generates representative permutations together with their cyclic orbits.
def gen_shapes(n):
    '''
        This function generates representative permutations together with their cyclic orbits.
    '''
    # domain of values
    vals = list(range(n+1))
    # enumerate permutations for positions 1..n-1 (fix value 0 at position 0)
    for tail in permutations(vals[2:]):  # (n-1)! iterations
        perm = (1,) + tail
        yield perm, orbit_of(perm, n)

# This function generates the shape orbits, optionally converting the generator to a list.
def determine_shape_list_by_reps(n, as_list=False):
    '''
        This function generates the shape orbits, optionally converting the generator to a list.
    '''

    gen = (orbit for _, orbit in gen_shapes(n))
    if as_list:
        return list(gen)
    return gen  # user can iterate, avoiding storing everything at once


# --------------------------------------------------------------------------- #
#                     High-level calculation functions                        #
# --------------------------------------------------------------------------- #

# This function generates the permutations and shapes, determines the path solution, and writes the associated results to files.
def determine_everything(_n, _cycle,cycle_word):
    '''
        This function generates the permutations and shapes, determines the path solution, and writes the associated results to files.
    '''
    _perm_list = generate_Sn_permutations(_n)

    _shape = determine_shape_list(_perm_list,_n)

    _compl_indlist, _sol_flag = determine_shape_corner_index_list_arb_path(_shape,_cycle,_perm_list)
    
    if _sol_flag == False:
        return False
    _assoc_perm, _shape_ngon = position_and_associated_shapes(_shape,_compl_indlist,_perm_list)

    subset_net_dict = {"permutations": _compl_indlist}
    with open(f"S{_n:.0f}/associated_permutation"+cycle_word+".json", 'w') as f:
        json.dump(subset_net_dict, f)

    subset_ngon_dict = {"ngons": _shape_ngon}
    with open(f"S{_n:.0f}/ngons"+cycle_word+".json", 'w') as f:
        json.dump(subset_ngon_dict, f)

    _subclass_shape, _subclass_shape_T = determine_subclasses(_assoc_perm,_perm_list)
    write_perms(_subclass_shape,f"S{_n:.0f}/subclasses"+cycle_word+".txt")
    write_perms(_subclass_shape_T,f"S{_n:.0f}/subclasses"+cycle_word+"_T.txt")

    write_vectors(_subclass_shape,f"S{_n:.0f}/vectors"+cycle_word+".txt")

# check solution of algebraic program using sympy 

# --------------------------------------------------------------------------- #
#                       SymPy verification functions                          #
# --------------------------------------------------------------------------- #

# This function checks a solution obtained from SymPy and writes the derived permutation, shape, and subclass data.
def check_sympy_solution(_n):
    '''
        This function checks a solution obtained from SymPy and writes the derived permutation, shape, and subclass data.
    '''
    print('\033[94m'+"Using the results from sympy for S_{}".format(_n)+'\033[0m')

    # --- Generate the permutation list. --- #
    print("Generating permutations... ",end="\n")
    time_perm = time.process_time()
    _perm_list = generate_Sn_permutations(_n)
    

    print(f"Time taken to generate permutations {(time.process_time()-time_perm):.4f} s")

    # --- Generate the shape list using representative permutations. --- #
    print("Generating shapes... ",end="\n")
    time_perm = time.process_time()
    _shape = determine_shape_list_by_reps(_n, as_list=True)


    # for large numbers -> save the data
    # shapes_dict = {"shapes": _shape}
    # with open("S{:.0f}/shapes.json".format(_n), 'w') as f:
    #     json.dump(shapes_dict, f)

    print(f"Time taken to generate shapes {(time.process_time()-time_perm):.4f} s")

    # --- Read the saved SymPy solution. --- #
    all_subsets_net, data_array = read_the_sympy_solution(_n,_perm_list)
    
    # --- Save the associated permutation index sets. --- #
    subset_net_dict = {"permutations": all_subsets_net}
    with open("S{:.0f}/associated_permutation_sympy.json".format(_n), 'w') as f:
        json.dump(subset_net_dict, f)

    # --- Determine the n-gons corresponding to the SymPy solution. --- #
    time_ngon = time.process_time()
    print("Determining n-gons... ",end="\n")
    determine_ngons(_n,_shape,all_subsets_net,_perm_list)
    print(f"Time taken to determine n-gons: {(time.process_time()-time_ngon):.4f} s")

    # --- Determine the path index sets using the SymPy solution. --- #
    _compl_indlist, _sol_flag = determine_shape_corner_index_list_arb_path_sympy(_n,all_subsets_net)
    if _sol_flag == False:
        return False

    # --- Save the associated permutation index sets. --- #
    subset_net_dict = {"permutations": all_subsets_net}
    with open("S{:.0f}/indices_sympy.json".format(_n), 'w') as f:
        json.dump(subset_net_dict, f)

    # --- Determine and write the permutation subclasses. --- #
    print("Determining subclasses...\n")
    _subclass_shape, _subclass_shape_T = determine_subclasses_sympy_solution(data_array)
    write_perms(_subclass_shape,"S{:.0f}/subclasses_sympy.txt".format(_n))
    write_perms(_subclass_shape_T,"S{:.0f}/subclasses_sympy_T.txt".format(_n))

    write_vectors(_subclass_shape,"S{:.0f}/vectors_sympy.txt".format(_n))

