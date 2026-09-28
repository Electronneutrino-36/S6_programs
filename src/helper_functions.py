# ---------- Helpful algebraic utilities ----------
import math
import itertools

def make_generators(n):
    """Return C (n-cycle as list) and T=(1 2) as list."""
    C = [((i % n) + 1) for i in range(1, n+1)]
    T = list(range(1, n+1))
    T[0], T[1] = T[1], T[0]
    return C, T

def compose(p, q):
    """Return permutation p∘q as lists (1-based images)."""
    return [p[q[i]-1] for i in range(len(p))]

def pow_perm(p, k):
    """Return p^k (k>=0)."""
    n = len(p)
    if k % n == 0:
        return list(range(1, n+1))
    res = list(range(1, n+1))
    k = k % n
    for _ in range(k):
        res = compose(p, res)
    return res

def perm_to_cycles(p):
    """Return cycle decomposition as list of tuples (excluding 1-cycles)."""
    n = len(p)
    seen = [False]*n
    cycles = []
    for i in range(n):
        if not seen[i] and p[i] != i+1:
            cur = i+1
            cyc = []
            while not seen[cur-1]:
                seen[cur-1] = True
                cyc.append(cur)
                cur = p[cur-1]
            cycles.append(tuple(cyc))
    return cycles

def order_of_perm(p):
    """Return order of permutation p (lcm of cycle lengths)."""
    cycles = perm_to_cycles(p)
    if not cycles:
        return 1
    lens = [len(c) for c in cycles]
    l = 1
    for a in lens:
        l = l * a // math.gcd(l, a)
    return l

# ---------- Convert textual word to permutation ----------
def parse_word_to_list(word_tokens, n):
    """
    Convert textual tokens like ["C^2","T","C^5"] or ["C","C","T","C","C","T"] into list-of-pairs form
    Accepts two styles: 'C^k' or a string of Cs like 'CC' (you already use those forms).
    Returns list of ('C',k) and ('T',1).
    """
    out = []
    for tok in word_tokens:
        if tok == "T":
            out.append(('T',1))
        elif tok.startswith("C^"):
            k = int(tok[2:])
            out.append(('C',k))
        elif all(ch == 'C' for ch in tok):
            out.append(('C', len(tok)))
        else:
            # fallback: try to parse like "CC" or "C^2"
            if 'C' in tok:
                out.append(('C', tok.count('C')))
            else:
                raise ValueError("Unknown token: " + str(tok))
    # simplify consecutive C-powers (mod n)
    C, T = make_generators(n)
    simp = []
    i = 0
    while i < len(out):
        if out[i][0] == 'C':
            total = out[i][1] % n
            j = i+1
            while j < len(out) and out[j][0] == 'C':
                total = (total + out[j][1]) % n
                j += 1
            if total % n != 0:
                simp.append(('C', total % n))
            i = j
        else:
            simp.append(out[i])
            i += 1
    return simp

def word_to_perm_from_tokens(word_tokens, n):
    """Return permutation g as list for given token list and n."""
    C, T = make_generators(n)
    parsed = parse_word_to_list(word_tokens, n)
    curr = list(range(1, n+1))
    for letter, exp in parsed:
        if letter == 'C':
            perm = pow_perm(C, exp)
        else:
            perm = T
        curr = compose(perm, curr)
    return curr, parsed

# ---------- Cheap algebraic checks ----------
def is_involution(p):
    return compose(p, p) == list(range(1, len(p)+1))

def is_single_transposition(p):
    cycles = perm_to_cycles(p)
    return len(cycles) == 1 and len(cycles[0]) == 2

def any_power_is_single_transposition(p):
    ordp = order_of_perm(p)
    for k in range(1, ordp+1):
        pk = pow_perm(p, k)
        if is_single_transposition(pk):
            return True
    return False

# ---------- Orbit enumeration and vector pairing (expensive) ----------
def enumerate_all_perms(n):
    return [list(p) for p in itertools.permutations(range(1, n+1))]

def enumerate_orbits_and_vectors(g, T):
    """
    g: permutation list (right-multiplication action)
    T: permutation list representing (1 2)
    Returns:
      - orbits: list of orbits (each orbit is list of perm-tuples)
      - vectors: list of vectors (each is union of 1 or 2 orbits, as indices into perms)
      - stats: dict with counts and size distribution
    """
    n = len(g)
    perms = enumerate_all_perms(n)
    tup_to_idx = {tuple(p): i for i, p in enumerate(perms)}
    visited = [False] * len(perms)
    orbits = []
    for i, p in enumerate(perms):
        if visited[i]:
            continue
        cur = p.copy()
        orbit = []
        while True:
            orbit.append(tuple(cur))
            visited[tup_to_idx[tuple(cur)]] = True
            cur = compose(g, cur)
            if visited[tup_to_idx[tuple(cur)]]:
                break
        orbits.append(orbit)

    # map each permutation tuple to its orbit index
    perm_to_orb = {}
    for oi, orb in enumerate(orbits):
        for tup in orb:
            perm_to_orb[tup] = oi

    # pair orbits by T
    paired = set()
    vectors = []
    for oi, orb in enumerate(orbits):
        if oi in paired:
            continue
        rep = orb[0]
        rep_after = tuple(compose(T, list(rep)))
        oj = perm_to_orb[rep_after]
        if oj == oi:
            # fixed orbit -> vector size = len(orb)
            vectors.append(list(orb))
            paired.add(oi)
        else:
            # union
            vectors.append(list(orb) + list(orbits[oj]))
            paired.add(oi); paired.add(oj)

    # stats
    vsizes = [len(v) for v in vectors]
    stats = {
        'num_orbits': len(orbits),
        'num_vectors': len(vectors),
        'vector_size_counts': {s: vsizes.count(s) for s in set(vsizes)},
        'total_covered': sum(vsizes),
    }
    return orbits, vectors, stats
