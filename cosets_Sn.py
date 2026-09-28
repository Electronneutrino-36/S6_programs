from sympy.combinatorics import Permutation, PermutationGroup
import itertools, math
import json
import time 
from collections import deque


def check_if_disjoint(_cosets):
    """
        Function to check whether the obtained cosets are disjoint.
    """
    CRED = '\033[91m'
    CGREEN = '\033[92m'
    CEND = '\033[0m'
    explicit_disjoint_check = []
    for c in range(0,len(_cosets)):
        for ci in range(c,len(_cosets)):
            if c == ci:
                continue
            else:
                explicit_disjoint_check.append(set(_cosets[c]).isdisjoint(set(_cosets[ci])))
    
    if all(set(explicit_disjoint_check)) == True:
        print(CGREEN+"All sets disjoint!"+CEND)
        return True
    else:
        print(CRED+"Found an overlap, word is not correct!"+CEND)
        return False

def tuple_from_perm(p):
    """Convert a SymPy Permutation p to a tuple representation (0-based)."""
    n = p.size
    return tuple(p(i) for i in range(n))

def permute_values(lst, perm):
    """Provide a dictionary mapping to easily access the permutations."""
    mapping = {i: perm(i) for i in range(perm.size)}
    return [mapping[v] for v in lst]

def compose_tuple_with_perm(p_tuple, h):
    """
       Return tuple for right multiplication p*h.
       p_tuple: tuple representation of a permutation (0-based).
       h: SymPy Permutation.
    """
    n = len(p_tuple)
    h = Permutation(h.array_form, size=n)
    mapping = {i: h(i) for i in range(h.size)}
    return tuple(mapping[v] for v in p_tuple)

def set_output_color(_state):
    '''Simple function to set the colour of the terminal output depending on the state.'''
    if _state == True:
        return '\033[92m'
    else:
        return '\033[91m'

def to_one_based(p_tuple):
    """Convert a tuple (0..n-1) -> (1..n)."""
    return tuple(x+1 for x in p_tuple)

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

def coset_reps(n, skip_check=False, _reorder=True):
    '''
        Function that determines the coset representations for a given S_n (with n even)
        - skip_check: is a flag to toggle the automatic check whether all the cosets 
        are disjoint (default is False)
        - reorder: does a sorting of the obtained cosets
    '''
    if skip_check == True:
        print("\033[91mSkipping element disjoint check \033[0m")

    # ------------------------------------------------------------------------#
    #                    Build the generators: C and T                        #
    # ------------------------------------------------------------------------#
    C = Permutation(list(range(1, n)) + [0])  # (0 1 ... n-1)
    T = Permutation(1, 0)                     # swap 0 and 1
    CT = C * T
    TC = T * C

    if n == 3:                                      # path for S_3
        W = T
    elif n != 3 and (n+1)%2 == 0:                   # path that works for S_5
        W = C**(n//2) * T * C**(n//2) * T * CT * T
    elif n%2 == 0:                                  # path for S_n and n even
        W = C**(n//2) * T

    # ------------------------------------------------------------------------#
    #                            Subgroup H                                   #
    # ------------------------------------------------------------------------#
    H = PermutationGroup([W,T])
    print(H)

    order_H = H.order()
    index = math.factorial(n) // order_H
    print(f"n={n}, |H|={order_H}, index={index}, length={len(H)}")

    # ------------------------------------------------------------------------#
    #       Enumerate H elements explicitly (OK if H is small)                #
    #       use the Schreier-Sims algorithm                                   #
    # ------------------------------------------------------------------------#
    H_elems = list(H.generate_schreier_sims())
    H_elems_array = list(H.generate_schreier_sims(af=True))


    print("Collecting coset representatives...")
    t_coset, t_coset_cpu = get_time()

    # ------------------------------------------------------------------------#
    #                   Collect coset representatives                         #
    # ------------------------------------------------------------------------#
    rep_list = 0
    covered = set()
    cosets = []
    for p in itertools.permutations(range(n)):   # all S_n as tuples
        if p in covered:
            continue

        rep_list += 1
        coset_members = []
        h_count = 0
        for h in H_elems:            
            ph = compose_tuple_with_perm(p, h)
            covered.add(ph)

            if _reorder == False:
                coset_members.append(ph)
            else:
                if h_count == 0:
                    coset_members.append(p)
                else:
                    continue

        cosets.append(coset_members)
        if rep_list >= index:
            break

        print(f"Progress: {(rep_list+1)/index*100:.4f}%    ",end="\r")

    del covered                                         # delete covered array to save memory for large n
    
    print("Number of permutations: ", math.factorial(n))
    print("Number of shapes: ", math.factorial(n)//n)
    print(f"Found {rep_list} coset representatives.")

    ec_t_coset, ec_t_coset_cpu = calc_time_taken(t_coset, t_coset_cpu)
    print(f"\033[94m-> Collecting cosets took {ec_t_coset:.4f} sec (wall), {ec_t_coset_cpu:.4f} sec (CPU)\033[0m")


    # ------------------------------------------------------------------------#
    #        Sanity check: explicitly check if all cosets are disjoint        #
    # ------------------------------------------------------------------------#
    if skip_check == False:
        print("Check if \033[94mcoset elements\033[0m are disjoint...")
        tcid, tcid_cpu_time = get_time()

        check_if_disjoint(cosets)

        ec_tcid, ec_tcid_cpu = calc_time_taken(tcid,tcid_cpu_time)
        print(f"\033[94m-> Checking if elements disjoint took {ec_tcid:.4f} sec (wall), {ec_tcid_cpu:.4f} sec (CPU)\033[0m")
        total_covered = sum(len(c) for c in cosets)

        disjoint = (len(set().union(*[set(c) for c in cosets])) == total_covered)
        full = (total_covered == math.factorial(n))

        CEND = '\033[0m'
        print("Collected cosets: ", len(cosets)," \n"
            "total covered: ", total_covered, " \n"
            "disjoint? "+set_output_color(disjoint)+"{}".format(disjoint)+CEND+" \n"
            "full? "+set_output_color(full)+"{}".format(full)+CEND)

    if _reorder == True:
        print("Reordering cosets representatives...")
        t_ord, t_ord_cpu_time = get_time()
        print(cosets[0])
        cosets_ordered = []

        cosets = deque(cosets)

        icc = 0
        l_cosets = len(cosets)
        while True:
            i_contents_ord = [cosets[0][0]]
            iccc = cosets[0][0]
            while True:
                iccc_new = compose_tuple_with_perm(iccc, C**(n//2))
                i_contents_ord.append(iccc_new)
    
                iccc_new = compose_tuple_with_perm(iccc_new, T)
                if iccc_new == cosets[0][0]:
                    cosets.popleft()
                    icc += 1
                    break
                else:
                    i_contents_ord.append(iccc_new)
                iccc = iccc_new

            cosets_ordered.append(i_contents_ord)

            print(f"Progress: {(icc+1)/l_cosets*100:.4f}%     ",end="\r")
            if icc >= l_cosets:
                break

        ec_t_ord, ec_t_ord_cpu = calc_time_taken(t_ord,t_ord_cpu_time)
        print(f"\033[94m-> Reordering the elements took {ec_t_ord:.4f} sec (wall), {ec_t_ord_cpu:.4f} sec (CPU)\033[0m")

    else:
        cosets_ordered = cosets[:]

    cosets_transformed = cosets_ordered[:]
    for i in range(0,len(cosets_transformed)):
        for j in range(0,len(cosets_transformed[i])):

            cosets_transformed[i][j] = to_one_based(cosets_transformed[i][j])

    out = {'n': n, '|H|': order_H, 'coset reps:': index, 'solution': cosets_transformed}
    with open("sol_s{:.0f}".format(n)+".json", 'w') as f:
        json.dump(out, f)
    return rep_list



def permutation_replace(lst, perm):
    """
    Value-based action of perm on list lst.
    Each element v in lst (an integer label) is sent to perm(v-1)+1 
    because SymPy uses 0-based indexing internally.
    """
    if perm.size < len(lst):
        perm = Permutation(perm.array_form, size=len(lst))
    return [perm(v) for v in lst]


def check_equality_of_data(data1, data2, keyword):
    with open(data1, 'r') as file:
        data_arr_1 = json.load(file)

    with open(data2, 'r') as file:
        data_arr_2 = json.load(file)

    true_array = [False]*len(data_arr_1[keyword])
    for i in range(0,len(data_arr_1[keyword])):
        if sorted(data_arr_1[keyword][i]) == sorted(data_arr_2[keyword][i]):
            true_array[i] = True
        else:
            print("No match for index ",i)
            print("Normal: ",data_arr_1[keyword][i])
            print("Sympy:  ",data_arr_2[keyword][i]) 

    true_check = all(x == True for x in true_array)
    print("Data matches: ",true_check)


# main file
n_expl = 6
reps = coset_reps(n_expl, True, _reorder=True)


#check_equality_of_data("sol_s6.json","sol_s6_test.json","solution")
