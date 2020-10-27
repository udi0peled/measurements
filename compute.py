#!/usr/bin/python3

B = 0
V = 0

Q = 256
P = 1024
N = 2048
N2 = 2*N

def set_params(background = 0, verification_batch = 0):
    B = background
    V = verification_batch

def E(exp, mod, B=0, V=0):
    if B > 0: return 0
    if V > 0: return 0

    # Paillier Cipher mod exp
    if (exp, mod) == (2048,4096) : return 7.00959
    if (exp, mod) == (256,4096)  : return 1.032578
    if (exp, mod) == (3*256,4096): return 2.748462

    # Paillier Plain mod exp
    if (exp, mod) == (256,2048) : return 0.279093
    # RingPedersen mod exp
    if (exp, mod) == (1024,1024): return 0.300123
    if (exp, mod) == (256,1024): return 0.084757

    return 0

def G(b=0):
    if b > 0: return 0
    return 0.476187

# Factorization of Paillier and Pedersen CRT, batch: Pedersen \lambda, Group Verification

def set_values(N, N2, P, Q, n):
    b = (E(N,N2, B), E(N,N2))
    c = (G(B), 2*G())
    d = (2*G(B) + G(), 7*G())
    e = (E(P,P) + 3*E(P,P,B) + E(Q,N) + E(N,N2) + G(B), \
         E(N,N2,0,V) + E(Q,N2) + 2*G() + 2*E(P,P)+E(Q,P) )
    f = (E(P,P) + 3*E(P,P,B) + E(Q,N) + E(N,N2,B) + 3*G(B), \
         E(N,N2,0,V) + E(Q,N2) + 5*G() + 2*E(P,P) + E(Q,P) )
    g = (E(3*Q, N2) + 2*E(N,N2, B) + 6*E(P,P,B) + 2*E(P,P) + G(B) + 2*E(Q,P), \
         2*E(N,N2,0,V) + 2*E(Q,N2) + E(3*Q, N2) + 2*G() + 4*E(P,P) + 2*E(Q,P) )

    presign = 2*b[0] + 7*G() + d[0] + (n-1)*(f[0] + f[1] + 4*b[0] + 2*E(Q,N2) + 2*g[0] + 2*g[1] + e[0] + e[1] + 2*b[1] + d[1])
    return presign
