from math import log, exp, pi, lgamma
from fpylll import IntegerMatrix, LLL, GSO
from vector import *
from random import choice, randint


def convert(path):
    with open(path, 'r') as file:
        lst = file.readlines()
    lst = lst[:-1]
    output = []
    for ligne in lst:
        ligne = ligne.replace('[', '')
        ligne = ligne.replace(']', '')
        ligne = ligne.replace('\n', '')
        ligne = ligne.split(' ')
        ligne = [int(elt) for elt in ligne]
        output.append(ligne)
    return output


def gh(B):
    n = B.ncols
    Bgso = GSO.Mat(B)
    _ = Bgso.update_gso()
    start_row = 0
    stop_row = -1
    log_det_B = Bgso.get_log_det(start_row, stop_row)
    a = (lgamma(n/2.+1)+log_det_B/2)/n
    return exp(a)/(pi**0.5)


def Sample_Init(B):
    n = B.ncols
    S = [Vector(n) for _ in range(n)]
    for i in range(n):
        for j in range(n):
            S[i][j] = B[i, j]
    return S


def Sample_Inflate(S):
    n = len(S[0].vec)
    v_temp = S[0]
    for i in range(n-1):
        for j in range(i+1, n):
            v_temp = S[i]-S[j]
            if v_temp < S[i] or v_temp < S[j]:
                S.append(v_temp)
    return S
        
        
def Sample_Reduce(S):
    n = len(S[0].vec)
    v_zero = Vector(n)
    v_start, v_temp = S[0], S[0]
    for _ in range(n):
        v_start = choice(S)
        for _ in range(n//2):
            v_choice = choice(S)
            v_temp = v_start - v_choice
            if (v_zero < v_temp < v_start or v_zero < v_temp < v_choice) and v_temp not in S:
                S.append(v_temp)
    return S


def Inflate_L4(B):
    S = Sample_Init(B)
    n = len(S[0].vec)
    v_min = S[0]
    S = Sample_Inflate(S)
    S = sorted(S, key = lambda vec : vec.sq_length)
    while S[0] < v_min:
        v_min = S[0]
        m = len(S)
        B = IntegerMatrix(m, n)
        for i in range(m):
            for j in range(n):
                B[i, j] = int(S[i].vec[j])
        B = LLL.reduction(B)
        S = [Vector(n) for _ in range(n)]
        for i in range(n):
            for j in range(n):
                S[i][j] = B[i+m-n, j]
        S = Sample_Inflate(S)
        S = sorted(S, key = lambda vec : vec.sq_length)
    return S


def Sample_L4(B):
    S = Sample_Init(B)
    n = len(S[0].vec)
    v_min = S[0]
    S = Sample_Reduce(S)
    S = sorted(S, key = lambda vec : vec.sq_length)
    while S[0] < v_min:
        v_min = S[0]
        m = len(S)
        B = IntegerMatrix(m, n)
        for i in range(m):
            for j in range(n):
                B[i, j] = int(S[i].vec[j])
        B = LLL.reduction(B)
        S = [Vector(n) for _ in range(n)]
        for i in range(n):
            for j in range(n):
                S[i][j] = B[i+m-n, j]
        S = Sample_Reduce(S)
        S = sorted(S, key = lambda vec : vec.sq_length)
    return S
        


###############################################
################ Randomization ################
###############################################

def gen_unimodular(n, density):
    sortie = [[0 for _ in range(n)] for _ in range(n)]
    liste_coefficients = [-1, 1]+[0 for _ in range(density)]
    for i in range(n):
        sortie[i][i] = choice([-1, 1])
        for j in range(i+1, n):
           sortie[i][j] = choice(liste_coefficients)
    # Permutations
    for i in range(n):
        j = randint(0, n-1)
        sortie[i], sortie[j] = sortie[j], sortie[i]
    return sortie


def mult(S, U):
    n = len(S[0].vec)
    B = IntegerMatrix(n, n)
    for i in range(n):
        for j in range(n):
            coeff = 0
            for k in range(n):
                coeff = coeff + S[k].vec[i]*U[j][k]
            B[j, i] = int(coeff)
    return B


def basis_reduction(S, density = 148):
    n = len(S[0].vec)
    U = gen_unimodular(n, density)
    return mult(S, U)
