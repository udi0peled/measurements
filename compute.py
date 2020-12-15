#!/usr/bin/python3

import itertools as it
import sys
import traceback


Q = 256
ELL = Q
EPS = 2*Q
ELLP = 5*Q
PEDERSEN = 1024
N = 2048
N2 = 2*N

VER_BATCH_FACTOR = 10

def E(calc_dict, calc_amount, calc_type, exp, mod, C=0, P=0, V=0):
    
    if C:
        exp //= 2
        mod //= 2
        res = E(calc_dict, 2*calc_amount, calc_type, exp, mod, C=0, P=P, V=V)
        return
    
    res = 0

    if P >= 0:
        base_time_ms = 7.929542
        if (exp, mod) == (2048,4096) : res = 1.000000 * base_time_ms
        if (exp, mod) == (2048,1024) : res = 0.081596 * base_time_ms
        if (exp, mod) == (1024,2048) : res = 0.141994 * base_time_ms
        if (exp, mod) == (1024,4096) : res = 0.510803 * base_time_ms
        if (exp, mod) == (1024,1024) : res = 0.042022 * base_time_ms
        if (exp, mod) == (4096,4096) : res = 1.973207 * base_time_ms
        if (exp, mod) == (512,4096) : res = 0.265187 * base_time_ms
        if (exp, mod) == (512,2048) : res = 0.074424 * base_time_ms
        if (exp, mod) == (512,1024) : res = 0.021511 * base_time_ms
        if (exp, mod) == (256,4096) : 
            #traceback.print_stack() 
            res = 0.140861 * base_time_ms
        if (exp, mod) == (256,2048) : res = 0.037966 * base_time_ms
        if (exp, mod) == (256,1024) : res = 0.011588 * base_time_ms
        if (exp, mod) == (128,4096) : res = 0.074498 * base_time_ms
        if (exp, mod) == (128,2048) : res = 0.021300 * base_time_ms
        if (exp, mod) == (128,1024) : res = 0.006852 * base_time_ms
        if (exp, mod) == (128,512) : res = 0.007956 * base_time_ms
        if (exp, mod) == (384,4096) : res = 0.202235 * base_time_ms
        if (exp, mod) == (384,2048) : res = 0.057764 * base_time_ms
        if (exp, mod) == (384,1024) : res = 0.016174 * base_time_ms
        if (exp, mod) == (384,512) : res = 0.024014 * base_time_ms
        if (exp, mod) == (896,512) : res = 0.052850 * base_time_ms
        if (exp, mod) == (384,512) : res = 0.023436 * base_time_ms
        if (exp, mod) == (768,4096) :
            #traceback.print_stack() 
            res = 0.394831 * base_time_ms
        if (exp, mod) == (768,2048) : res = 0.108155 * base_time_ms
        if (exp, mod) == (768,1024) : res = 0.040745 * base_time_ms
        if (exp, mod) == (768,1024) : res = 0.032987 * base_time_ms
        if (exp, mod) == (1792,1024) : res = 0.070891 * base_time_ms
        if (exp, mod) == (1280,1024) : res = 0.052689 * base_time_ms
        if (exp, mod) == (1536,4096) : res = 0.811380 * base_time_ms

        if res == 0: print(f'error {exp}, {mod}, {C}, {P}, {V}')

    if V == 1:
        res /= VER_BATCH_FACTOR
        
    if P == 1 and (exp == 2048 and mod == 4096) or (exp == 1024 and mod == 2048):
        res /= 10

    curr_amount, time = calc_dict[calc_type].get((exp,mod, C, P, V), (0, 0))    
    calc_dict[calc_type][(exp,mod, C, P, V)] = (curr_amount + calc_amount, time + calc_amount*res)

def Ped(calc_dict, calc_amount, is_private, exp1_bits, exp2_bits, exp1_rand=0, exp2_rand=0):
    pool1 = False #int((exp1_rand>0) and (pool>0))
    pool2 = False #int((exp2_rand>0) and (pool>0))
    
    to_crt = int((is_private>0) and (crt>0))
    to_ver = int((is_private>0) and (ver>0))
    
#    if is_private == 1 and exp1_rand == 0 and exp2_rand == 0:
#            E(calc_dict, calc_amount, "pedersen", max(exp1_bits,exp2_bits), PEDERSEN, C=to_crt, V=to_ver)

    E(calc_dict, calc_amount, "pedersen", exp1_bits, PEDERSEN, C=to_crt, P=pool1, V=to_ver)
    E(calc_dict, calc_amount, "pedersen", exp2_bits, PEDERSEN, C=to_crt, P=pool2, V=to_ver)

def G(calc_dict, calc_amount, P=0,V=0):
    res = 0.476187
    if P == 1: res = 0
    if V == 1: res /= VER_BATCH_FACTOR

    curr_amount, time = calc_dict["group"].get((256, 0, P, V), (0, 0))
    calc_dict["group"][(256, 0, P, V)] = (curr_amount + calc_amount, time + calc_amount*res)

# Factorization of Paillier and Pedersen CRT, batch: Pedersen \lambda, Group Verification

def paillier(b, d, a):
    if b == 0:
        E(d, a, "paillier", N,N2,P=pool)
    else:
        E(d, a, "paillier", N,N2, C=crt)
    
def sch(b, d, a):
    if b == 0:
        G(d, a, P=pool)
    else:
        G(d, a)
        G(d, a, V=ver)

def ddh(b, d, a):
    if b == 0:
        G(d, 2*a, P=pool)
        G(d, 2*a)
    else:
        G(d, 5*a)
        G(d, 2*a, V=ver)

def log(b, d, a):
    if b == 0:
        Ped(d, a, 0,ELL,ELL+PEDERSEN,0,1)
        Ped(d, a, 0,ELL+EPS,ELL+EPS+PEDERSEN,1,1)
        E(d, a, "paillier", Q,N) 
        E(d, a, "paillier", N,N2,C=crt, P=pool)
        G(d, a, P=pool)
    else:
        E(d, a, "paillier", N,N2,V=ver)
        E(d, a, "paillier", Q,N2)  ###256
        Ped(d, a, 1, ELL+EPS, ELL+EPS+PEDERSEN)
        E(d, a, "pedersen", Q,PEDERSEN, C=crt)
        G(d, a)
        G(d, a, V=ver)

def Rddh(b, d, a):
    if b == 0:
        Ped(d, a, 0,ELL,ELL+PEDERSEN,0,1)
        Ped(d, a, 0, ELL+EPS, ELL+EPS+PEDERSEN,1,1)
        E(d, a, "paillier", Q,N, C=crt)
        E(d, a, "paillier", N,N2,C=crt, P=pool)
        G(d, 3*a, P=pool)
    else:
        E(d, a, "paillier", N,N2,V=ver)
        E(d, a, "paillier", Q,N2)    ###256
        E(d, a, "pedersen", Q,PEDERSEN,C=crt)
        Ped(d, a, 1, ELL+EPS, ELL+EPS+PEDERSEN)
        G(d, 3*a)
        G(d, 2*a, V=ver)

def affg(b, d, a):
    if b == 0:
        E(d, a, "paillier", ELL+EPS,N2)  ###768
        E(d, a, "paillier", N,N2,C=crt,P=pool)
        E(d, a, "paillier", N,N2,P=pool)
        Ped(d, a, 0,ELL+EPS,ELL+EPS+PEDERSEN,1,1)
        Ped(d, a, 0,ELLP+EPS,ELL+EPS+PEDERSEN,1,1)
        Ped(d, a, 0,ELL,ELL+PEDERSEN,0,1)
        Ped(d, a, 0,ELLP,ELL+PEDERSEN,0,1)
        E(d, a, "pedersen", Q,PEDERSEN,C=crt)
        E(d, a, "pedersen", Q,PEDERSEN)
        G(d, a, P=pool)
    else:
        E(d, a, "paillier", ELL+EPS,N2,C=crt)
        E(d, a, "paillier", Q, N2,C=crt)
        E(d, a, "paillier", N, N2,C=crt,V=ver)
        E(d, a, "paillier", N, N2,V=ver)
        E(d, a, "paillier", Q, N2) ###256
        Ped(d, a, 1, ELL+EPS, ELL+EPS+PEDERSEN)
        Ped(d, a, 1, ELLP+EPS, ELL+EPS+PEDERSEN)
        E(d, 2*a, "pedersen", Q, PEDERSEN)
        G(d, a)
        G(d, a, V=ver)
        

def ECDSA_presign(d, n):
    paillier(0, d, 2)
    G(d, 7)
    ddh(0, d, 1)
    Rddh(0, d, n-1)
    Rddh(1, d, n-1)
    paillier(0, d, 4*(n-1))
    E(d, 2*(n-1), "paillier", Q, N2) ###256
    affg(0, d, 2*(n-1))
    affg(1, d, 2*(n-1))
    log(0, d, n-1)
    log(1, d, n-1)
    paillier(1, d, 2*(n-1))
    ddh(1, d, n-1)

def set_values(n):
    calc_dict = {"paillier": {}, "pedersen": {}, "group": {}}
    ECDSA_presign(calc_dict, n)

    return calc_dict

def print_dict(d):
    sum = 0
    for t, v in d.items():
        print(t)
        curr_sum = 0
        for expmod, res in v.items():
            print(f'{expmod}: {res}')
            curr_sum += res[1]
        print(f'{curr_sum}\n')
        sum += curr_sum
    print(f'Total: {sum}')


crt = 0
pool = 0
ver = 0

if len(sys.argv) >= 2: crt  = int(sys.argv[1])
if len(sys.argv) >= 3: pool = int(sys.argv[2])
if len(sys.argv) >= 4: ver  = int(sys.argv[3])

print(f'crt exponentiation(C) = {crt}, r^N pool sampling(P) = {pool}, batch verification(V) = {ver}\n')
print("(exp, mod C, P, V): (amount, total time)")
total_group = 0
print_dict(set_values(3))

