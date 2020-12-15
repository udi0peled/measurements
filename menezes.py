#!/usr/bin/python3

def f(k, exp_bitlength=2048):
    t = exp_bitlength
    l = t//k
    return (t + l*(2**k-1)//2**k + 2**k-3, t + l*(2**k-1)//2**k + 2**(k-1)) 