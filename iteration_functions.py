import itertools
from itertools import permutations
from functools import partial
import sys
from helper_functions import *
import multiprocessing
from sympy.combinatorics import Permutation
from multiprocessing import Pool, cpu_count, Process, Manager
import json
import time
import numpy as np

def get_time():
    t = time.time()
    t_cpu = time.process_time()
    return t, t_cpu

def calc_time_taken(_t, _t_cpu):
    t_a = time.time()
    t_a_cpu = time.process_time()

    dt = t_a - _t
    dt_cpu = t_a_cpu - _t_cpu
    return dt, dt_cpu

sys.setrecursionlimit(1000) 

def write_file(combs, name):
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

def reorder_from_idx(idx, a):
    return a[idx:] + a[:idx]

def cyclic_perm(a):
    return [partial(reorder_from_idx, i) for i in range(len(a))]

def P12(a):
    at = list(a[:])
    ind1 = at.index(1)
    ind2 = at.index(2)

    at[ind1],at[ind2] = at[ind2],at[ind1]
    at = tuple(at)
    return at


def Pij(a,i,j):
    at = list(a[:])
    ind1 = at.index(i)
    ind2 = at.index(j)
    at[ind1],at[ind2] = at[ind2],at[ind1]
    at = tuple(at)
    return at

def P123456(a):
    aorg = a
    for i in range(0,5):
        aorg = Pij(aorg,5-i,6-i)
    return aorg

def Pcycl(a):
    aorg = a
    dima = len(a)
    dimam1 = len(a)-1
    for i in range(0,dimam1):
        aorg = Pij(aorg,dimam1-i,dima-i)
    return aorg

def determine_path(perm,cycle,permlist):
    perms_of_met_hex = []

    pnew = perm
    counter = 0
    ellist = []
    while True:
        for i in range(0,cycle):
            pnew = P123456(pnew)
        ellist.append(permlist.index(pnew))
        pnew = P12(pnew)
        ellist.append(permlist.index(pnew))

        perms_of_met_hex.append([pnew])
        pinter = pnew[:]
        for i in range(0,5):
            pinter = P123456(pinter)
            perms_of_met_hex[counter].append(pinter)

        counter += 1
        #test if we are back at the original element
        if pnew == perm:
            break
    
    return perms_of_met_hex, ellist

def determine_general_path(perm,cycle,permlist):
    perms_of_met_shape = []

    pnew = perm
    counter = 0
    ellist = []
    while True:
        for i in range(0,cycle):
            pnew = Pcycl(pnew)
        ellist.append(permlist.index(pnew))
        pnew = P12(pnew)
        ellist.append(permlist.index(pnew))

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


def determine_arbitrary_path(perm,cycle,_permlist):
    perms_of_met_shape = []

    pnew = perm
    counter = 0
    ellist = []
    while True:
        for i in range(0,len(cycle)):
            for k in range(0,len(cycle[i])):
                if cycle[i][k] == "C":
                    pnew = Pcycl(pnew)
                elif cycle[i][k] == "T":
                    pnew = P12(pnew)

            ellist.append(_permlist.index(pnew))

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

def find_shape(perm,cycle,shape):
    nodes = [perm]
    testperm = perm
    for i in range(0,cycle-1):
        testperm = Pcycl(testperm)
        nodes.append(testperm)
    
    nodes.sort()
    if nodes not in shape:
        shape.append(nodes)

def find_shape_parallel(perm,cycle):
    nodes = [perm]
    testperm = perm
    for i in range(0,cycle-1):
        testperm = Pcycl(testperm)
        nodes.append(testperm)
    
    nodes.sort()
    return nodes

def determine_shape_list(_perm_list,_n):
    shape = []

    for i in range(0,len(_perm_list)):
        find_shape(_perm_list[i],_n,shape)
    print("Found {} shapes with {} corners from the {} permutations!".format(len(shape),_n,len(_perm_list)))
    return shape

def find_hexagons(perm,cycle,hexagons):
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


def determine_single_hex(perm,cycle):
    nodes = [perm]
    testperm = perm
    for i in range(0,cycle-1):
        testperm = P123456(testperm)
        nodes.append(testperm)
    return nodes

def get_hex(perm,hexagons):
    for i in range(0,len(hexagons)):
        if perm in hexagons[i]:
            return i
        else:
            continue



# generate permutations of
def generate_Sn_permutations(n_dim):
    sn_list = []
    for i in range(1,n_dim+1):
        sn_list.append(i)
    sn_perm_list = list(itertools.permutations(sn_list))

    del sn_list

    return sn_perm_list


# determine indices of shape corners in path
def determine_shape_corner_index_list(shapelist, cycl, i_permlist):
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

def test_indices(_indlist,_sorted_indlist):
    flag = False
    for i in _indlist:
        count = 0
        for j in range(0,len(_sorted_indlist)):
            if i in _sorted_indlist[j]:
                flag = True
                return flag


def find_exact_cover(universe, subsets):
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

# --- DLX Matrix Builder (unchanged) ---

def get_dlx_matrix(universe_elements, subsets):
    universe = {j: Header(j) for j in universe_elements}
    root = Header("root")
    
    root.right = root
    root.left = root
    
    for j in universe_elements:
        j_header = universe[j]
        j_header.left = root.left
        j_header.right = root
        root.left.right = j_header
        root.left = j_header

    for subset in subsets:
        row_nodes = []
        for element in subset:
            if element not in universe:
                continue
            
            col_header = universe[element]
            node = Node(col_header, row_data=subset)
            row_nodes.append(node)
            
            node.up = col_header.up
            node.down = col_header
            col_header.up.down = node
            col_header.up = node
            col_header.size += 1

        if row_nodes:
            for i in range(len(row_nodes) - 1):
                current_node = row_nodes[i]
                next_node = row_nodes[i + 1]
                current_node.right = next_node
                next_node.left = current_node
            
            row_nodes[0].left = row_nodes[-1]
            row_nodes[-1].right = row_nodes[0]
            
    return root

# --- Core DLX Operations (unchanged) ---

def cover(c):
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

def uncover(c):
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


def solve_dlx_iterative(root):
  """
  An iterative implementation of the DLX algorithm.
  This version uses an explicit stack to manage state and avoids RecursionError.
  """
  solution_nodes = []
  stack = []
  
  c = root.right # Start with the first column

  loop_counter = 0
  while True:
      # Step 1: Check for solution
      if c == root:
          return [node.row_data for node in solution_nodes]
      
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

      # Step 3: Handle a column with no available rows (backtrack)
      if c.down == c:
          if not stack:
            elements_recovered = [node.row_data for node in solution_nodes]
            print(f"Recovered only {len(elements_recovered)} elements. Excited with {loop_counter} counts.")
            return None # No solution
          c = stack.pop()
          uncover(c)
          continue
      
      # Step 4: Cover column and try the first row
      cover(c)
      r = c.down
      solution_nodes.append(r)
      stack.append(c)
      
      # Step 5: Cover other columns in the chosen row
      j = r.right
      while j != r:
          cover(j.col)
          j = j.right
      
      c = root.right # Continue search with the first column
      loop_counter += 1
      
  return None  

  
# determine indices of shape corners in arbitrary path
def determine_shape_corner_index_list_arb_path(shapelist, cycl, i_permlist):

    sol_flag = False
    i_compl_permlist = []

    print("Determining path for all permutations.")

    for ii in range(0,len(i_permlist)):
        perms_of_met_hex, _indlist_, _counter = determine_arbitrary_path(i_permlist[ii],cycl,i_permlist)
        i_compl_permlist.append(_indlist_)
    

    print("Found all {} paths. Takes {} steps with {} indices to get back".format(len(i_compl_permlist),_counter,len(_indlist_)) )

    # reduce number of sets
    # unique_sets = {frozenset(l) for l in i_compl_permlist}

    # # Convert the unique frozensets back into a list of lists.
    # final_lists = [list(s) for s in unique_sets]
    # print("Reduced the number of lists from {} to {}".format(len(i_compl_permlist),len(final_lists)))
    # print(final_lists)

    universe = set(range(len(i_permlist)))

    all_subsets = i_compl_permlist


    print("Finding solutions")
    dlx_matrix_time = time.process_time()
    matrix_root = get_dlx_matrix(universe, all_subsets)
    print(f"Time taken for get_dlx_matrix {(time.process_time() - dlx_matrix_time):.4f} s")

    dlx_solver_time = time.process_time()
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

def read_the_sympy_solution(_n,_permlist):
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

# determine indices of shape corners in arbitrary path -> check sympy path
def determine_shape_corner_index_list_arb_path_sympy(_n, _all_subsets_net):
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

    # solution = find_exact_cover(universe,all_subsets)
    if solution:
        print(f"\033[92mFound a solution with\033[0m {len(solution)} \033[92mlists.\033[0m")
        sol_flag = True

    else:
        raise Exception(CRED+"No solution found. -> Please provide a different path!"+CEND)
        


    return solution, sol_flag


def test_T_relation(a_list,b,_perml):
    perm_shapes = []
    for i in a_list:
        perm_shapes.append(_perml[i])
    for il in perm_shapes:
        if P12(il) == _perml[b]:
            return True

# determine vectors of permutations
def determine_subclasses(_assoc_perm,_permlist):
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


def test_T_relation_sympy(a_list,b,):
    for il in a_list:
        if P12(il) == b:
            return True
        
# determine vectors of permutations
def determine_subclasses_sympy_solution(_data):
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


def determine_ngons(_n,_shapes,_indlist,_permlist):
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


def position_and_associated_shapes(shapes, _indlist, _permlist):
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


# find the indices of the hexagon corners involved in the path
def determine_hexagon_corner_index_list(hexlist,i_permlist):
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

def write_perms(_permlist,filename):
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

def write_vectors(_permlist,filename):
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


def predetermine_perm_shape(_n):
    _perm_list = generate_Sn_permutations(_n)

    _shape = determine_shape_list(_perm_list,_n)

    return _perm_list, _shape


def determine_rest(_shape,_cycle,_perm_list):
    _compl_indlist, _sol_flag = determine_shape_corner_index_list_arb_path(_shape,_cycle,_perm_list)
    if _sol_flag == False:
        return False

def test_for_solution(_n,_combinations,_permlist_):
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

def generate_all_odd_combinations(_permlist_):
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

def generate_all_even_combinations(_permlist_):
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

def brute_force_path(_n,_type):
    _permlist, _shapes = predetermine_perm_shape(_n)

    
    if _type == "odd":
        _combs = generate_all_odd_combinations(_permlist)
    else:
        _combs = generate_all_even_combinations(_permlist)

    print("Paths to test:")
    for icomb in _combs:
        print(icomb)
    test_for_solution(_n,_combs,_permlist)    

def determine_everything_wo_write(_n, _cycle):
    _perm_list = generate_Sn_permutations(_n)

    _shape = determine_shape_list(_perm_list,_n)

    _compl_indlist, _sol_flag = determine_shape_corner_index_list_arb_path(_shape,_cycle,_perm_list)
    if _sol_flag == False:
        return False


def orbit_of(perm, n):
    nodes = []
    cur = perm
    for _ in range(n):
        nodes.append(cur)
        cur = Pcycl(cur)
    return nodes  # immutable representation

def gen_shapes(n):
    # domain of values
    vals = list(range(n+1))
    # enumerate permutations for positions 1..n-1 (fix value 0 at position 0)
    for tail in permutations(vals[2:]):  # (n-1)! iterations
        perm = (1,) + tail
        yield perm, orbit_of(perm, n)

def determine_shape_list_by_reps(n, as_list=False):

    gen = (orbit for _, orbit in gen_shapes(n))
    if as_list:
        return list(gen)
    return gen  # user can iterate, avoiding storing everything at once


def determine_everything(_n, _cycle,cycle_word):
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
def check_sympy_solution(_n):
    print('\033[94m'+"Using the results from sympy for S_{}".format(_n)+'\033[0m')

    print("Generating permutations... ",end="\n")
    time_perm = time.process_time()
    _perm_list = generate_Sn_permutations(_n)
    

    print(f"Time taken to generate permutations {(time.process_time()-time_perm):.4f} s")

    print("Generating shapes... ",end="\n")
    time_perm = time.process_time()
    _shape = determine_shape_list_by_reps(_n, as_list=True)


    # for large numbers -> save the data
    # shapes_dict = {"shapes": _shape}
    # with open("S{:.0f}/shapes.json".format(_n), 'w') as f:
    #     json.dump(shapes_dict, f)

    print(f"Time taken to generate shapes {(time.process_time()-time_perm):.4f} s")

    all_subsets_net, data_array = read_the_sympy_solution(_n,_perm_list)
    
    subset_net_dict = {"permutations": all_subsets_net}
    with open("S{:.0f}/associated_permutation_sympy.json".format(_n), 'w') as f:
        json.dump(subset_net_dict, f)

    time_ngon = time.process_time()
    print("Determining n-gons... ",end="\n")
    determine_ngons(_n,_shape,all_subsets_net,_perm_list)
    print(f"Time taken to determine n-gons: {(time.process_time()-time_ngon):.4f} s")

    _compl_indlist, _sol_flag = determine_shape_corner_index_list_arb_path_sympy(_n,all_subsets_net)
    if _sol_flag == False:
        return False

    subset_net_dict = {"permutations": all_subsets_net}
    with open("S{:.0f}/indices_sympy.json".format(_n), 'w') as f:
        json.dump(subset_net_dict, f)

    print("Determining subclasses...\n")
    _subclass_shape, _subclass_shape_T = determine_subclasses_sympy_solution(data_array)
    write_perms(_subclass_shape,"S{:.0f}/subclasses_sympy.txt".format(_n))
    write_perms(_subclass_shape_T,"S{:.0f}/subclasses_sympy_T.txt".format(_n))

    write_vectors(_subclass_shape,"S{:.0f}/vectors_sympy.txt".format(_n))

# -------------------------------------------------------------------------------
# functions to generate arbitraty paths to brute force the path finding
def generate_c_power_t(power):
    """Generate C^power T notation"""
    cstr = ""
    for i in range(0,power):
        cstr += "C"
    return f"T"+cstr

def generate_c_power(power):
    """Generate C^power T notation"""
    cstr = ""
    for i in range(0,power):
        cstr += "C"
    return cstr

def generate_even_unique_combinations(order):
    """Generate unique combinations for each pattern type"""
    
    # Define the possible values for exponents
    exponents = []
    for i in range(1,order+1):
        exponents.append(i)
    

    results = {
        'pattern_1': [],  # [C^a''', T]
        'pattern_2': [],  # [C^a'', T, C^a''', T]
        'pattern_3': [],  # [C^a', T, C^a'', T, C^a''', T]
        'pattern_4': [],   # [C^a, T, C^a', T, C^a'', T, C^a''', T]
    }
    
    # Pattern 1: [C^a''' T, T] - only needs 1 exponent (a''')
    for a_triple_prime in [1,2,3,4,5]:
        pattern_1 = [generate_c_power(a_triple_prime),"T"]
        results['pattern_1'].append(pattern_1)

    for i_aa in range(1,len(exponents)):
        for i_aaa in range(0,i_aa+1): 
            a_double_prime = exponents[i_aa]
            a_triple_prime = exponents[i_aaa]
            pattern_2 = [generate_c_power(a_double_prime), "T",
                        generate_c_power(a_triple_prime), "T"]
            results['pattern_2'].append(pattern_2)
    
    # Pattern 3: [C^a', T, C^a'', T, C^a''' T, T] - needs 3 exponents (a', a'', a''')
    for i_a in range(0,len(exponents)):
        for i_aa in range(0,i_a+1):
            for i_aaa in range(0,i_aa+1): 
                a_prime = exponents[i_a]
                a_double_prime = exponents[i_aa]
                a_triple_prime = exponents[i_aaa]
                pattern_3 = [generate_c_power(a_prime), "T",
                            generate_c_power(a_double_prime), "T",
                            generate_c_power(a_triple_prime), "T"]
                results['pattern_3'].append(pattern_3)
    
    # Pattern 4: [C^a, T, C^a', T, C^a'', T, C^a''' T, T] - needs 4 exponents (a, a', a'', a''')
    for i in range(0,len(exponents)):
        for i_a in range(0,i+1):
            for i_aa in range(0,i_a+1):
                for i_aaa in range(0,i_aa+1): 
                    a = exponents[i]
                    a_prime = exponents[i_a]
                    a_double_prime = exponents[i_aa]
                    a_triple_prime = exponents[i_aaa]
                    pattern_4 = [generate_c_power(a), "T",
                                generate_c_power(a_prime), "T",
                                generate_c_power(a_double_prime), "T",
                                generate_c_power(a_triple_prime), "T"]
                    results['pattern_4'].append(pattern_4)
    
    return results

def generate_unique_combinations():
    """Generate unique combinations for each pattern type"""
    
    # Define the possible values for exponents
    exponents = [1, 2, 3, 4]
    
    results = {
        'pattern_1': [],  # [C^a''' T, T]
        'pattern_2': [],  # [C^a'', T, C^a''' T, T]
        'pattern_3': [],  # [C^a', T, C^a'', T, C^a''' T, T]
        'pattern_4': [],   # [C^a, T, C^a', T, C^a'', T, C^a''' T, T]
        'pattern_5': [],   # [C^a, T, C^a', T, C^a'', T, C^a''' T, T]
        'pattern_6': [],   # [C^a, T, C^a', T, C^a'', T, C^a''' T, T]
        'pattern_7': [],   # [C^a, T, C^a', T, C^a'', T, C^a''' T, T]
        'pattern_8': [],   # [C^a, T, C^a', T, C^a'', T, C^a''' T, T]
        'pattern_9': [],   # [C^a, T, C^a', T, C^a'', T, C^a''' T, T]
        'pattern_10': []
    }
    
    # Pattern 1: [C^a''' T, T] - only needs 1 exponent (a''')
    for a_triple_prime in [1,2]:
        pattern_1 = [generate_c_power_t(a_triple_prime), "T"]
        results['pattern_1'].append(pattern_1)
    
    # Pattern 2: [C^a'', T, C^a''' T, T] - needs 2 exponents (a'', a''')
    for i_a in range(0,len(exponents)):
        for i_aa in range(0,i_a):
            a_double_prime = exponents[i_a]
            a_triple_prime = exponents[i_aa]
            pattern_2 = [generate_c_power(a_double_prime), "T", 
                        generate_c_power_t(a_triple_prime), "T"]
            results['pattern_2'].append(pattern_2)
    
    # Pattern 3: [C^a', T, C^a'', T, C^a''' T, T] - needs 3 exponents (a', a'', a''')
    for i_a in range(0,len(exponents)):
        for i_aa in range(0,i_a+1):
            for i_aaa in range(0,i_aa): 
                a_prime = exponents[i_a]
                a_double_prime = exponents[i_aa]
                a_triple_prime = exponents[i_aaa]
                pattern_3 = [generate_c_power(a_prime), "T",
                            generate_c_power(a_double_prime), "T",
                            generate_c_power_t(a_triple_prime), "T"]
                results['pattern_3'].append(pattern_3)
    
    # Pattern 4: [C^a, T, C^a', T, C^a'', T, C^a''' T, T] - needs 4 exponents (a, a', a'', a''')
    for i in range(0,len(exponents)):
        for i_a in range(0,i+1):
            for i_aa in range(0,i_a+1):
                for i_aaa in range(0,i_aa): 
                    a = exponents[i]
                    a_prime = exponents[i_a]
                    a_double_prime = exponents[i_aa]
                    a_triple_prime = exponents[i_aaa]
                    pattern_4 = [generate_c_power(a), "T",
                                generate_c_power(a_prime), "T",
                                generate_c_power(a_double_prime), "T",
                                generate_c_power_t(a_triple_prime), "T"]
                    results['pattern_4'].append(pattern_4)

    for i in range(0,len(exponents)):
        for i_a in range(0,i+1):
            for i_aa in range(0,i_a):
                for i_aaa in range(0,i_a+1): 
                    a = exponents[i]
                    a_prime = exponents[i_a]
                    a_double_prime = exponents[i_aa]
                    a_triple_prime = exponents[i_aaa]
                    pattern_5 = [generate_c_power(a), "T",
                                generate_c_power(a_prime), "T",
                                generate_c_power_t(a_double_prime), "T",
                                generate_c_power(a_triple_prime), "T"]
                    results['pattern_5'].append(pattern_5)

    for i in range(0,len(exponents)):
        for i_a in range(0,i+1):
            for i_aa in range(0,i_a):
                for i_aaa in range(0,i_a+1): 
                    for i_aaaa in range(0,i_aaa):
                        a = exponents[i]
                        a_prime = exponents[i_a]
                        a_double_prime = exponents[i_aa]
                        a_triple_prime = exponents[i_aaa]
                        a_quad_prime = exponents[i_aaaa]
                        pattern_6 = [generate_c_power(a), "T",
                                    generate_c_power(a_prime), "T",
                                    generate_c_power_t(a_double_prime), "T",
                                    generate_c_power(a_triple_prime), "T",
                                    generate_c_power_t(a_quad_prime), "T"]
                        results['pattern_6'].append(pattern_6)

    for i in range(0,len(exponents)):
        for i_aa in range(0,i):
            for i_a in range(0,i+1):
                for i_aaa in range(0,i_a+1): 
                    for i_aaaa in range(0,i_aaa):
                        a = exponents[i]
                        a_prime = exponents[i_a]
                        a_double_prime = exponents[i_aa]
                        a_triple_prime = exponents[i_aaa]
                        a_quad_prime = exponents[i_aaaa]
                        pattern_7 = [generate_c_power(a), "T",
                                    generate_c_power_t(a_double_prime), "T",
                                    generate_c_power(a_prime), "T",
                                    generate_c_power(a_triple_prime), "T",
                                    generate_c_power_t(a_quad_prime), "T"]
                        results['pattern_7'].append(pattern_7)
    
    for i in range(0,len(exponents)):
        for i_a in range(0,i+1):
            for i_aa in range(0,i_a+1):
                for i_aaa in range(0,i_a+1): 
                    for i_aaaa in range(0,i_aaa):
                        a = exponents[i]
                        a_prime = exponents[i_a]
                        a_double_prime = exponents[i_aa]
                        a_triple_prime = exponents[i_aaa]
                        a_quad_prime = exponents[i_aaaa]
                        pattern_8 = [generate_c_power(a), "T",
                                    generate_c_power(a_prime), "T",
                                    generate_c_power(a_double_prime), "T",
                                    generate_c_power(a_triple_prime), "T",
                                    generate_c_power_t(a_quad_prime), "T"]
                        results['pattern_8'].append(pattern_8)


    for i in range(0,len(exponents)):
        for i_aa in range(0,1):
            for i_a in range(0,i+1):
                for i_aaaaa in range(0,1):
                    for i_aaa in range(0,i_a+1): 
                        for i_aaaa in range(0,1):
                            a = exponents[i]
                            a_prime = exponents[i_a]
                            a_double_prime = exponents[i_aa]
                            a_triple_prime = exponents[i_aaa]
                            a_quad_prime = exponents[i_aaaa]
                            a_pent_prime = exponents[i_aaaaa]
                            pattern_9 = [generate_c_power(a), "T",
                                        generate_c_power_t(a_double_prime), "T",
                                        generate_c_power(a_prime), "T",
                                        generate_c_power_t(a_pent_prime), "T",
                                        generate_c_power(a_triple_prime), "T",
                                        generate_c_power_t(a_quad_prime), "T"]
                            results['pattern_9'].append(pattern_9)


    
    return results

def print_results(results):
    """Print all generated patterns in a readable format"""
    
    pattern_names = {
        'pattern_1': '[C^a\'\'\' T, T]',
        'pattern_2': '[C^a\'\', T, C^a\'\'\' T, T]',
        'pattern_3': '[C^a\', T, C^a\'\', T, C^a\'\'\' T, T]',
        'pattern_4': '[C^a, T, C^a\', T, C^a\'\', T, C^a\'\'\' T, T]'
    }
    
    for pattern_key, pattern_name in pattern_names.items():
        print(f"\n{pattern_name}:")
        print("=" * 50)
        
        for i, combination in enumerate(results[pattern_key], 1):
            print(f"{i:2d}: {combination}")

def get_all_combinations_as_array():
    """Return all combinations as a single nested array"""
    results = generate_unique_combinations()
    all_combinations = []
    
    # Add all patterns to a single list
    for pattern_key in ['pattern_1', 'pattern_2', 'pattern_3', 'pattern_4', 'pattern_5', 'pattern_6', 'pattern_7', 'pattern_8', 'pattern_9']:
        all_combinations.extend(results[pattern_key])
    
    return all_combinations


def get_all_even_combinations_as_array(order):
    """Return all combinations as a single nested array"""
    results = generate_even_unique_combinations(order)
    all_combinations = []

    pattern_list = []
    for i in range(1,order+1):
        pattern_list.append('pattern_'+str(i))
    
    # Add all patterns to a single list
    for pattern_key in pattern_list:
        all_combinations.extend(results[pattern_key])
    
    return all_combinations
