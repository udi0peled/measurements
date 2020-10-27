#!/usr/bin/python3

Q = 256
ELL = Q
EPS = 2*Q
ELLP = 5*Q
PEDERSEN = 1024
N = 2048
N2 = 2*N

pool = 0
ver = 0
crt = 0

def show_params():
    print(f'r^N pool sampling = {pool}, batch verification = {ver}, crt exponentiation = {crt}')

def E(exp, mod, C=0, P=0, V=0):
    
    if P > 0: return 0
    
    if V > 0: return  E(exp,mod,C,P,V=0)/10
    
    if C:
        exp //= 2
        mod //= 2
    
    res = 0
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
    if (exp, mod) == (256,4096) : res = 0.140861 * base_time_ms
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
    if (exp, mod) == (768,4096) : res = 0.394831 * base_time_ms
    if (exp, mod) == (768,2048) : res = 0.108155 * base_time_ms
    if (exp, mod) == (768,1024) : res = 0.040745 * base_time_ms
    if (exp, mod) == (768,1024) : res = 0.032987 * base_time_ms
    if (exp, mod) == (1792,1024) : res = 0.070891 * base_time_ms
    if (exp, mod) == (1280,1024) : res = 0.052689 * base_time_ms
    if (exp, mod) == (1536,4096) : res = 0.811380 * base_time_ms

    if res == 0: print(f'error {exp}, {mod}, {C}, {P}, {V}')
    
    if C: res *= 2
    return res

def Ped(is_private, exp1_bits, exp2_bits, exp1_rand=0, exp2_rand=0):
    pool1 = (exp1_rand>0) and (pool>0)
    pool2 = (exp2_rand>0) and (pool>0)
    
    to_crt = (is_private>0) and (crt>0)
    to_ver = (is_private>0) and (ver>0)
    
    if is_private and exp1_rand == 0 and exp2_rand ==0:
            res = E(max(exp1_bits,exp2_bits), PEDERSEN, C=to_crt, V=to_ver)

    res = E(exp1_bits, PEDERSEN, C=to_crt, P=pool1, V=to_ver) \
        + E(exp2_bits, PEDERSEN, C=to_crt, P=pool2, V=to_ver)
        
    return res

def G(P=0,V=0):
    #if b > 0: return 0
    return 0.476187

# Factorization of Paillier and Pedersen CRT, batch: Pedersen \lambda, Group Verification

def set_values(n):
    paillier = (E(N,N2,P=pool), E(N,N2))
    
    sch = (G(pool), 2*G())
    
    ddh = (2*G(pool) + 2*G(), 7*G())
    
    log = (Ped(0,ELL,ELL+PEDERSEN,0,1) + Ped(0,ELL+EPS,ELL+EPS+PEDERSEN,1,1) + E(Q,N) + E(N,N2,C=crt, P=pool) + G(pool), \
            E(N,N2,V=ver) + E(Q,N2) + 2*G() + Ped(1, ELL+EPS, ELL+EPS+PEDERSEN) + E(Q,PEDERSEN, C=crt) )
    
    Rddh = (Ped(0,ELL,ELL+PEDERSEN,0,1) + Ped(0, ELL+EPS, ELL+EPS+PEDERSEN,1,1) + E(Q,N, C=crt) + E(N,N2,C=crt, P=pool) + 3*G(pool), \
            E(N,N2,V=ver) + E(Q,N2) + 5*G() + Ped(1, ELL+EPS, ELL+EPS+PEDERSEN) + E(Q,PEDERSEN,C=crt) )
    
    affg = (E(ELL+EPS,N2) + E(N,N2,C=crt,P=pool) + E(N,N2,P=pool) + Ped(0,ELL+EPS,ELL+EPS+PEDERSEN,1,1) +  Ped(0,ELLP+EPS,ELL+EPS+PEDERSEN,1,1) + \
            Ped(0,ELL,ELL+PEDERSEN,0,1) + Ped(0,ELLP,ELL+PEDERSEN,0,1) + E(Q,PEDERSEN,C=crt) + E(Q,PEDERSEN) + G(pool), \
            E(ELL+EPS,N2,C=crt) + E(Q, N2,C=crt) + E(N,N2,C=crt,V=ver) + E(N,N2,V=ver) + E(Q,N2) + 2*G() + Ped(1, ELL+EPS, ELL+EPS+PEDERSEN) + Ped(1, ELLP+EPS, ELL+EPS+PEDERSEN) + 2*E(Q,PEDERSEN) )

    presign = 2*paillier[0] + 7*G() + ddh[0] + (n-1)*(Rddh[0] + Rddh[1] + 4*paillier[0] + 2*E(Q,N2) + 2*affg[0] + 2*affg[1] + log[0] + log[1] + 2*paillier[1] + ddh[1])
    return presign

import itertools as it

for pool in range(2):
    for ver in range(2):
        for crt in range(2):
            show_params()
            print(set_values(3))

