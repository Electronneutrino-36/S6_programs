# --------------------------------------------------------------------------- #
#                     Necessary libraries                                     #
# --------------------------------------------------------------------------- #
import numpy as np
import math
from numpy.linalg import multi_dot
import sympy as sp
import pickle

def calc_dim(dia):
    '''
        Function that takes a S_n young diagram and calculates its dimension. 
        This is done by calculating the horizontal hook length (hlen) and 
        vertical hook length (vlen) for each box in the diagram, and then 
        calculating the total hook length by summing hl_helper=hlen+vlen and 
        then hooklength = prod_i hl_helper(i)
    '''
    ni = len(dia[:,0])
    hlen = np.zeros((ni,ni))

    # ---------------------------------------------------------------- #
    #        Calculate horizontal hook length for each box             #
    # ---------------------------------------------------------------- #
    for i in range(0,ni):
        for k in range(0,ni):
            isum = 0
            for j in range(k,ni):
                isum += dia[i,j]
            hlen[i,k] = (isum)

    # ---------------------------------------------------------------- #
    #           Calculate vertical hook length for each box            #
    # ---------------------------------------------------------------- #
    vlen = np.zeros((ni,ni))
    for j in range(0,ni):
        for k in range(0,ni):
            vsum = 0
            for i in range(k+1,ni):
                vsum += dia[i,j]
            # print(vsum)
            vlen[k,j] = (vsum)

    # ---------------------------------------------------------------- #
    #   Introduce hooklength helper function for following product     #
    # ---------------------------------------------------------------- #
    hooklengthm = np.array(hlen)+np.array(vlen)

    Sndim = math.factorial(ni)
    hooklen = 1.0

    # ---------------------------------------------------------------- #
    #       Calculate total hooklength h for the given diagram         #
    # ---------------------------------------------------------------- #    
    for i in range(0,ni):
        for j in range(0,ni):
            if hooklengthm[i,j] != 0.0:
                hooklen *= hooklengthm[i,j]

    # ---------------------------------------------------------------- #
    #            Calculate dimension d of the diagram d=n!/h           #
    # ---------------------------------------------------------------- #    
    d = Sndim/hooklen
    return d


def calc_axial_dist(dia,i,j):
    '''
        This function determines the axial distances of two numbers i and j 
        in a given diagram dia. It searches for the entry of the numbers in the 
        diagram of given shape and stores them as "position vectors" \vec{indi} 
        and \vec{indj} and. Then it calculates the vector difference 
        \vec{diff} = \vec{indi}-\vec{indj} and obtains the axial distance as 
        given below.
    '''

    # ---------------------------------------------------------------- #
    #   Find position of numbers i and j in given diagram shape        #
    # ---------------------------------------------------------------- #    
    indih = np.where(dia == i)
    indjh = np.where(dia == j)

    # ---------------------------------------------------------------- #
    #    Store position of numbers i and j in "position vectors"       #
    # ---------------------------------------------------------------- #    
    indi = np.array([indih[0][0],indih[1][0]])
    indj = np.array([indjh[0][0],indjh[1][0]])

    # ---------------------------------------------------------------- #
    #                   Calculate distance vector                      #
    # ---------------------------------------------------------------- #
    diff = indi - indj

    # ---------------------------------------------------------------- #
    #                     Calculate axial distance                     #
    # ---------------------------------------------------------------- #    
    axdist = 0
    for i in range(0,len(diff)):
        axdist += (-1)**(i+1)*diff[i]
    return axdist

def index_exchange(M,j):
    '''
        Simple function to exchange the indices j with the next j+1 in matrix M
    '''
    Mtest = M[:]
    Mtest[j],Mtest[j+1] = Mtest[j+1],Mtest[j]
    return Mtest

def calc_coeff(dia):
    '''
        This function calculates the coefficients for the two unit vectors 
        e1 and e2 of a given diagram (irrep).
    '''

    n = len(dia[:,0])
    rho = []

    # ---------------------------------------------------------------- #
    #         Calculate axial distance rho for given diagram           #
    # ---------------------------------------------------------------- #
    for i in range(0,n-1):
        rhon = calc_axial_dist(dia,i+2,i+1)
        rho.append(rhon)

    # ---------------------------------------------------------------- #
    #                       Define unit vectors                        #
    # ---------------------------------------------------------------- #
    e1 = np.array([1,0])
    e2 = np.array([0,1])

    # ---------------------------------------------------------------- #    
    # Determine coefficients in the two directions of the unit vectors #
    # ---------------------------------------------------------------- #    
    coefs = [] 
    for i in range(0,len(rho)):
        Dtimese1 = 1/(rho[i])*e1 + sp.sqrt(rho[i]**2-1)/np.abs((rho[i]))*e2
        coefs.append(Dtimese1)
    return coefs, rho

def calc_coeff_general(M,M_base):
    '''
        This function is a generalization of the calc_coef from above. 
        It takes a Yamanouchi symbol M and the set of Yamanouchi symbols 
        M_base. From M it determines to corresponding Young tableaux and 
        from that 
        it defines the unit vector e1 in the direction that the matrix 
        M takes in the vector space spanned by M_base. It then tries to 
        exchange the indices of M to find another element of M_base to 
        span a two dimensional subspace of M_base.
    '''
    dim = len(M_base)

    # ---------------------------------------------------------------- #
    #  Determine Young tableaux corresponding to Yamanouchi symbol M   #
    # ---------------------------------------------------------------- #
    dia = determine_y_tableaux(M)

    n = len(dia[:,0])
    rho = []

    # ---------------------------------------------------------------- #
    #         Calculate the axial distance for that diagram            #
    # ---------------------------------------------------------------- #
    for i in range(0,n-1):
        rhon = calc_axial_dist(dia,i+2,i+1)
        rho.append(rhon)
    
    # ---------------------------------------------------------------- #
    #               Set unit vector in direction of M                  #
    # ---------------------------------------------------------------- #
    e1 = np.zeros(dim)
    inde1 = M_base.index(M)
    e1[inde1] = 1

    # ---------------------------------------------------------------- #    
    #     Determine second direction e2 and calculate coefficients     #
    # ---------------------------------------------------------------- #    
    coefs = [] 
    for i in range(0,len(rho)):
        e2 = np.zeros(dim)
        if np.abs(rho[i]) != 1:
            M_ex = index_exchange(M,i)
            if M_ex in M_base:
                inde2 = M_base.index(M_ex)
            else:
                continue
        else:
            inde2 = inde1
        e2[inde2] = 1

        # ---------------------------------------------------------------- # 
        #                    Calculate coefficients                        #
        # ---------------------------------------------------------------- # 
        Dtimese1 = 1/(rho[i])*e1 + sp.sqrt(rho[i]**2-1)/np.abs((rho[i]))*e2
        coefs.append(Dtimese1)
    return coefs, rho


def turn_to_matrix(coeffsmat,ind,matrix):
    '''
        Simple function to turn the coefficient matrix to a matrix. Obsolete...
    '''
    counti1neq0 = 0
    for i in range(0,len(coeffsmat)):
        arr = coeffsmat[i][ind]
        if arr[1] != 0:
            counti1neq0 += 1
        iar = 0
        if arr[1] == 0:
            matrix[i,i] = arr[0]
        else:
            if counti1neq0 == 1:
                matrix[i,i] = arr[0]
                matrix[i,i+1] = arr[1]
            elif counti1neq0 == 2:
                matrix[i,i-1] = arr[1]
                matrix[i,i] = arr[0]  
            
    # print(matrix)

def determine_y_tableaux(M):
    '''
        Function to determine Young tableaux from Yamanouchi symbol M
    '''
    coord = []
    n = len(M)
    counter = np.zeros(n)

    for i in range(0,len(M)):
        counter[M[i]-1] += 1
        if counter[M[i]-1] == 1:
            coord_i = [M[i]-1,0]
        else:
            coord_i = [M[i]-1,int(counter[M[i]-1]-1)]
        coord.append(coord_i)

    # ---------------------------------------------------------------- # 
    #       Create nxn matrix which has entries 1 in the               #
    #       shape of the Young tableaux for M                          #
    # ---------------------------------------------------------------- # 
    Sy_test = np.zeros((n,n))
    for i in range(0,len(Sy_test)):
        Sy_test[coord[i][0],coord[i][1]] = i+1

    return Sy_test

def determine_M_basis(dim,M_basis):
    '''
        Function to automatically determine the Yamanouchi symbol basis given 
        the dimension (dim) of the irrep.
    '''

    # ---------------------------------------------------------------- # 
    #      Determine Young tableaux and calculate its coefficients     #
    # ---------------------------------------------------------------- # 
    for i in range(0,dim):
        rhoi = calc_coeff(determine_y_tableaux(M_basis[i]))[1]
        for j in range(0,len(rhoi)):
            Mtest = M_basis[i][:]
            if np.abs(rhoi[j]) != 1:
                Mtest[j],Mtest[j+1] = Mtest[j+1],Mtest[j]
                if Mtest in M_basis:
                    continue
                else:
                    M_basis.append(Mtest)

def compile_matrix(coef,n):
    '''
        This function stores the individual matrices for a partition n from the 
        coefficient array.
    '''
    Mat = []
    for i in range(0,n-1):
        Mat.append([])
    for i in range(0,len(coef)):
        for j in range(0,len(Mat)):
            Mat[j].append(list(coef[i][j]))
    return Mat

def write_matrix(matrix,name):
    '''
        Function to write the found matrix to a file with filename "name"
    '''
    matrixfile = open(name,"w")
    for i in range(0,len(matrix)):
        for j in range(0,len(matrix[i])):
            matrixstr = str(matrix[i,j])
            if matrixstr.find("sqrt"):
                matrixstr = matrixstr.replace("sqrt(","Sqrt[")
                matrixstr = matrixstr.replace(")","]")
            if j <= len(matrix[i])-2:
                matrixfile.write("{}  ".format(matrixstr))
            else:
                matrixfile.write("{}".format(matrixstr))
        if i <= len(matrix)-2:
            matrixfile.write("\n")
    matrixfile.close()

    # ---------------------------------------------------------------- # 
    #                 Also write it as a pickle file                   #
    # ---------------------------------------------------------------- # 
    with open(name[:-3]+"pickle", 'wb') as outf:
        outf.write(pickle.dumps(matrix))


def ordertest(A):
    '''
        Function that tests the ordering of a vector A.
    '''

    return all(A[i] >= A[i+1] for i in range(len(A)-1))

def determine_combination(comb,comblist,N,n):
    '''
        Function that automatically determines all the possible Yamanouchi 
        symbols by starting from an inital Yamanouchi symbol comb[0] (usually 
        n times the number 1) and 
        summing up the constituents of the Yamanouchi symbols that 
        come before. The Yamanouchi symbols are stored in comblist. n is the 
        n of S_n and N=2 tells us that always two elements of the Yamanouchi 
        symbols should be summed together. This is in principle variable.
    '''

    # ---------------------------------------------------------------- # 
    #        Usually start with the n times 1 Yamanouchi symbol        #
    # ---------------------------------------------------------------- # 
    for i in range(0,len(comb)-(N-1)):
        nsum = 0
        for j in range(0,N):
            nsum += comb[i+j]
        ncomb = []
        for k in range(0,i):
            ncomb.append(comb[k])
        ncomb.append(nsum)
        for k in range(i+(N),len(comb)):
            ncomb.append(comb[k])
        if ordertest(ncomb) == True:
            if sum(ncomb) == n:
                if ncomb in comblist:
                    continue
                else:
                    comblist.append(ncomb)


def young_tableaux_from_dim(dim,n):
    '''
        Function to determine the Young diagram from Yamanouchi symbol dim.
    '''
    Sy = np.zeros((n,n))
    for i in range(0,len(dim)):
        for j in range(0,dim[i]):
            Sy[i,j] = 1
    return Sy

def determine_initial_Yamanouchi_symbol(S):
    '''
        This function determines the inital Yamanouchi symbol given a Young tableaux S.
    '''
    nl = len(S[:,0])
    numbers = list(np.linspace(1,nl,nl))
    for i in range(0,nl):
        for j in range(0,nl):
            if S[i,j] == 1:
                S[i,j] = numbers[0]
                numbers.pop(0)
     
    M_init = []
    for i in range(0,nl):
        indnum = np.where(S == i+1)
        M_init.append(indnum[0][0]+1)
    return M_init

def determine_matrix(mu, mu_list,Sndims,Snmatrices):
    '''
        Given a partition mu, the list of partitions mu_list 
        the dimensions of the irreps Sndims and the Young corresponding 
        Young tableaus Snmatrices, this function determines the full 
        basis of Yamanouchi symbols for the partition mu and constructs 
        the matrices.
    '''
    nl = len(Snmatrices[0][:,0])

    # ---------------------------------------------------------------- # 
    #         Get index of partition in the list of partitions         #
    # ---------------------------------------------------------------- # 
    ind = mu_list.index(mu)
    name = ""
    for i in range(0,len(mu)):
        name += str(mu[i])

    # ---------------------------------------------------------------- # 
    #    Get dimension of the Young diagram for the obtained index     #
    # ---------------------------------------------------------------- # 
    dim = Sndims[ind]

    # ---------------------------------------------------------------- # 
    #    Determine inital Yamanouchi symol for matrix construction     #
    # ---------------------------------------------------------------- # 
    Mi = determine_initial_Yamanouchi_symbol(Snmatrices[ind])

    # ---------------------------------------------------------------- # 
    #             Determine basis of Yamanouchi symbols                #
    # ---------------------------------------------------------------- # 
    M_basis = [Mi]
    determine_M_basis(dim,M_basis)
    print("Dimension: ", dim, " mu: ", mu)

    # ---------------------------------------------------------------- # 
    #         Sort the basis to avoid the matrices not matching        #
    # ---------------------------------------------------------------- # 
    M_basis.sort()

    print("M: ",M_basis)

    # ---------------------------------------------------------------- # 
    #                     Write the basis to file                      #
    # ---------------------------------------------------------------- # 
    basisfile = open("S"+str(nl)+"/matrices/"+name+"/basis.dat","w")
    basisfile.write("Basis\n: {}".format(M_basis))
    basisfile.close()

    # ---------------------------------------------------------------- # 
    #             Write the dimension of the basis to file             #
    # ---------------------------------------------------------------- # 
    dimensionfile = open("S"+str(nl)+"/matrices/"+name+"/dimension.dat","w")
    dimensionfile.write("Dimension: {}".format(dim))
    dimensionfile.close()

    # ---------------------------------------------------------------- # 
    #                   Calculate coefficents and axial                #
    #                 distances for the constructed basis              #
    # ---------------------------------------------------------------- # 
    coef = []
    rho = []
    for i in range(0,len(M_basis)):
        calccoef = calc_coeff_general(M_basis[i],M_basis)
        coef.append(calccoef[0])
        rho.append(calccoef[1])

    # ---------------------------------------------------------------- # 
    #    Construct explicit matrix from coef array for partition mu    #
    # ---------------------------------------------------------------- # 
    mat = compile_matrix(coef,sum(mu))
    product = multi_dot(mat)

    # ---------------------------------------------------------------- # 
    #                   Write the matrices to files                    #
    # ---------------------------------------------------------------- # 
    trafos = []
    for i in range(1,nl):
        trafos.append(str(i)+str(i+1))
    for i in range(0,len(mat)):
        write_matrix(np.array(mat[i]),"S"+str(nl)+"/matrices/"+name+"/matrix_"+name+"_"+trafos[i]+".txt")
    write_matrix(product,"S"+str(nl)+"/matrices/"+name+"/matrix_"+name+".txt")

def compute_c_size(_part, _n):
    '''
        Compute size of comjugacy class belonging to partition _part
    '''
    name = ""
    for i in range(0,len(_part)):
        name += str(_part[i])
    _prod = 1
    _red_part = list(dict.fromkeys(_part))
    for i in _red_part:
        _c = _part.count(i)
        _prod *= (i**_c * math.factorial(_c))


    _csize = math.factorial(_n)/_prod

    ciszefile = open("S"+str(_n)+"/matrices/"+name+"/csize.dat","w")
    ciszefile.write("{}".format(_csize))
    ciszefile.close()
    return(_csize)