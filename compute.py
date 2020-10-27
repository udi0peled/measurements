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
    print(f'pool = {pool}, ver = {ver}, crt = {crt}')

def E(exp, mod, C=0, P=0, V=0):
    
    if P > 0: return 0
    
    if V > 0: return E(exp,mod,C,P,V=0)/10
    
    if C:
        exp //= 2
        mod //= 2
    
    res = 0
    base_time_ms = 5.685242
    if (exp, mod) == (2048,4096) : res = 1.000000 * base_time_ms
    if (exp, mod) == (2048,1024) : res = 0.077083 * base_time_ms
    if (exp, mod) == (1024,2048) : res = 0.136546 * base_time_ms
    if (exp, mod) == (1024,4096) : res = 0.514042 * base_time_ms
    if (exp, mod) == (1024,1024) : res = 0.039685 * base_time_ms
    if (exp, mod) == (4096,4096) : res = 1.980670 * base_time_ms
    if (exp, mod) == (512,4096) : res = 0.255801 * base_time_ms
    if (exp, mod) == (512,2048) : res = 0.076914 * base_time_ms
    if (exp, mod) == (512,1024) : res = 0.021245 * base_time_ms
    if (exp, mod) == (256,4096) : res = 0.140735 * base_time_ms
    if (exp, mod) == (256,2048) : res = 0.037493 * base_time_ms
    if (exp, mod) == (256,1024) : res = 0.011757 * base_time_ms
    if (exp, mod) == (128,4096) : res = 0.069756 * base_time_ms
    if (exp, mod) == (128,2048) : res = 0.019113 * base_time_ms
    if (exp, mod) == (128,1024) : res = 0.006222 * base_time_ms
    if (exp, mod) == (128,512) : res = 0.002515 * base_time_ms
    if (exp, mod) == (384,4096) : res = 0.196598 * base_time_ms
    if (exp, mod) == (384,2048) : res = 0.053895 * base_time_ms
    if (exp, mod) == (384,1024) : res = 0.017076 * base_time_ms
    if (exp, mod) == (768,4096) : res = 0.377192 * base_time_ms
    if (exp, mod) == (768,2048) : res = 0.102272 * base_time_ms
    if (exp, mod) == (768,1024) : res = 0.030858 * base_time_ms
    if (exp, mod) == (768,1024) : res = 0.030961 * base_time_ms
    if (exp, mod) == (1792,1024) : res = 0.075015 * base_time_ms
    if (exp, mod) == (1280,1024) : res = 0.055187 * base_time_ms
    if (exp, mod) == (1536,4096) : res = 0.761429 * base_time_ms

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

pool = 0
ver = 0
crt = 0

print(set_values(3))

pool = 1
ver = 1
crt = 1

print(set_values(3))

crt = 0
pool = 0

print(set_values(3))

