#include <assert.h>
#include <stdio.h>
#include "openssl/bn.h"
#include <time.h>

#define udi_check_top(a)
#define udi_pollute(a)
#define BN_FLG_FIXED_TOP 0

int UDI_mod_mul_montgomery(BIGNUM *r, const BIGNUM *a, const BIGNUM *b, BN_MONT_CTX *mont, BN_CTX *ctx);
int udi_mul_mont_fixed_top(BIGNUM *r, const BIGNUM *a, const BIGNUM *b, BN_MONT_CTX *mont, BN_CTX *ctx);
int UDI_from_montgomery(BIGNUM *ret, const BIGNUM *a, BN_MONT_CTX *mont, BN_CTX *ctx);
int udi_from_mont_fixed_top(BIGNUM *ret, const BIGNUM *a, BN_MONT_CTX *mont, BN_CTX *ctx);

int PM_mod_mul_montgomery(BIGNUM *r, BIGNUM *rq, const BIGNUM *a, const BIGNUM *aq, const BIGNUM *b, const BIGNUM *bq, BN_MONT_CTX *mont, const BIGNUM *N, BN_CTX *ctx);
int pm_mul_mont_fixed_top(BIGNUM *r, BIGNUM *rq, const BIGNUM *a, const BIGNUM *aq, const BIGNUM *b, const BIGNUM *bq, BN_MONT_CTX *mont, BN_CTX *ctx);
int PM_from_montgomery(BIGNUM *ret, BIGNUM *retq, const BIGNUM *a, BN_MONT_CTX *mont, BN_CTX *ctx);
int pm_from_mont_fixed_top(BIGNUM *ret, BIGNUM *retq, const BIGNUM *a, BN_MONT_CTX *mont, BN_CTX *ctx);

clock_t start;
clock_t diff;
double  single_ms;

int debug_print = 0;
void debug_print_on() { debug_print = 1; }
void debug_print_off() { debug_print = 0; }

void printBIGNUM(const char * prefix, const BIGNUM *bn, const char * suffix) {
  char *bn_str = BN_bn2dec(bn);
  printf("%s%s%s", prefix, bn_str, suffix);
  free(bn_str);
}

/**
 * 
 * 
 *  OpenSSL replacements
 * 
 */

struct bignum_st {
    BN_ULONG *d;                /* Pointer to an array of 'BN_BITS2' bit
                                 * chunks. */
    int top;                    /* Index of last used d +1. */
    /* The next are internal book keeping for bn_expand. */
    int dmax;                   /* Size of the d array. */
    int neg;                    /* one if the number is negative */
    int flags;
};


/* Used for montgomery multiplication */
struct bn_mont_ctx_st {
    int ri;                     /* number of bits in R */
    BIGNUM RR;                  /* used to convert to montgomery form,
                                   possibly zero-padded */
    BIGNUM N;                   /* The modulus */
    BIGNUM Ni;                  /* R*(1/R mod N) - N*Ni = 1 (Ni is only
                                 * stored for bignum algorithm) */
    BN_ULONG n0[2];             /* least significant word(s) of Ni; (type
                                 * changed with 0.9.9, was "BN_ULONG n0;"
                                 * before) */
    int flags;
};

void udi_correct_top(BIGNUM *a)
{
    BN_ULONG *ftl;
    int tmp_top = a->top;

    if (tmp_top > 0) {
        for (ftl = &(a->d[tmp_top]); tmp_top > 0; tmp_top--) {
            ftl--;
            if (*ftl != 0)
                break;
        }
        a->top = tmp_top;
    }
    if (a->top == 0)
        a->neg = 0;
    a->flags &= ~BN_FLG_FIXED_TOP;
    udi_pollute(a);
}

/** 
 * 
 * 
 *  Udi Montgomery OpenSSL simple computation (no asm or word), 2x slower
 * 
 */

int UDI_MONT_CTX_set(BN_MONT_CTX *mont, const BIGNUM *mod, BN_CTX *ctx)
{
    int i, ret = 0;
    BIGNUM *Ri, *R;

    if (BN_is_zero(mod))
        return 0;

    BN_CTX_start(ctx);
    if ((Ri = BN_CTX_get(ctx)) == NULL)
        goto err;
    R = &(mont->RR);            /* grab RR as a temp */
    if (!BN_copy(&(mont->N), mod))
        goto err;               /* Set N */
    if (BN_get_flags(mod, BN_FLG_CONSTTIME) != 0)
        BN_set_flags(&(mont->N), BN_FLG_CONSTTIME);
    mont->N.neg = 0;

    {                           /* bignum version */
        mont->ri = BN_num_bits(&mont->N);
        BN_zero(R);
        if (!BN_set_bit(R, mont->ri))
            goto err;           /* R = 2^ri */
        /* Ri = R^-1 mod N */
        if ((BN_mod_inverse(Ri, R, &mont->N, ctx)) == NULL)
            goto err;
        if (!BN_lshift(Ri, Ri, mont->ri))
            goto err;           /* R*Ri */
        if (!BN_sub_word(Ri, 1))
            goto err;
        /*
         * Ni = (R*Ri-1) / N
         */
        if (!BN_div(&(mont->Ni), NULL, Ri, &mont->N, ctx))
            goto err;
    }

    /* setup RR for conversions */
    BN_zero(&(mont->RR));
    if (!BN_set_bit(&(mont->RR), mont->ri * 2))
        goto err;
    if (!BN_mod(&(mont->RR), &(mont->RR), &(mont->N), ctx))
        goto err;

    for (i = mont->RR.top, ret = mont->N.top; i < ret; i++)
        mont->RR.d[i] = 0;
    mont->RR.top = ret;
    mont->RR.flags |= BN_FLG_FIXED_TOP;

    ret = 1;
 err:
    BN_CTX_end(ctx);
    return ret;
}

int UDI_mod_mul_montgomery(BIGNUM *r, const BIGNUM *a, const BIGNUM *b, BN_MONT_CTX *mont, BN_CTX *ctx)
{
    int ret = udi_mul_mont_fixed_top(r, a, b, mont, ctx);

    udi_correct_top(r);
    udi_check_top(r);

    return ret;
}

int udi_mul_mont_fixed_top(BIGNUM *r, const BIGNUM *a, const BIGNUM *b, BN_MONT_CTX *mont, BN_CTX *ctx)
{
    BIGNUM *tmp;
    int ret = 0;
    int num = mont->N.top;

    if ((a->top + b->top) > 2 * num)
        return 0;

    BN_CTX_start(ctx);
    tmp = BN_CTX_get(ctx);
    if (tmp == NULL)
        goto err;

    udi_check_top(tmp);
    if (a == b) {
        if (!BN_sqr(tmp, a, ctx))
            goto err;
    } else {
        if (!BN_mul(tmp, a, b, ctx))
            goto err;
    }
    /* reduce from aRR to aR */
    if (!UDI_from_montgomery(r, tmp, mont, ctx))
        goto err;
        
    ret = 1;
 err:
    BN_CTX_end(ctx);
    return ret;
}

int UDI_to_montgomery(BIGNUM *r, const BIGNUM *a, BN_MONT_CTX *mont,
                     BN_CTX *ctx)
{
    return UDI_mod_mul_montgomery(r, a, &(mont->RR), mont, ctx);
}

int UDI_from_montgomery(BIGNUM *ret, const BIGNUM *a, BN_MONT_CTX *mont, BN_CTX *ctx)
{
    int retn;

    retn = udi_from_mont_fixed_top(ret, a, mont, ctx);
    udi_correct_top(ret);
    udi_check_top(ret);

    return retn;
}

int udi_from_mont_fixed_top(BIGNUM *ret, const BIGNUM *a, BN_MONT_CTX *mont, BN_CTX *ctx)
{
    int retn = 0;
    
    BIGNUM *t1, *t2;

    BN_CTX_start(ctx);
    t1 = BN_CTX_get(ctx);
    t2 = BN_CTX_get(ctx);
    if (t2 == NULL)
        goto err;

    if (!BN_copy(t1, a))
        goto err;
    BN_mask_bits(t1, mont->ri);

    if (!BN_mul(t2, t1, &mont->Ni, ctx))
        goto err;
    BN_mask_bits(t2, mont->ri);

    if (!BN_mul(t1, t2, &mont->N, ctx))
        goto err;
    if (!BN_add(t2, a, t1))
        goto err;
    if (!BN_rshift(ret, t2, mont->ri))
        goto err;

    if (BN_ucmp(ret, &(mont->N)) >= 0) {
        if (!BN_usub(ret, ret, &(mont->N)))
            goto err;
    }
    retn = 1;
    udi_check_top(ret);
 err:
    BN_CTX_end(ctx);

    return retn;
}

/** 
 * 
 *  Paillier-Montgomery function
 * 
 */


int PM_mod_mul_montgomery(BIGNUM *r, BIGNUM *rq, const BIGNUM *a, const BIGNUM *aq, const BIGNUM *b, const BIGNUM *bq, BN_MONT_CTX *mont, const BIGNUM *N, BN_CTX *ctx)
{
    int ret = pm_mul_mont_fixed_top(r, rq, a, aq, b, bq, mont, ctx);

    udi_correct_top(r);
    udi_correct_top(rq);
    udi_check_top(r);
    udi_check_top(rq);

    return ret;
}

int pm_mul_mont_fixed_top(BIGNUM *r, BIGNUM *rq, const BIGNUM *a, const BIGNUM *aq, const BIGNUM *b, const BIGNUM *bq, BN_MONT_CTX *mont, BN_CTX *ctx)
{
    BIGNUM *sum;
    BIGNUM *tmp;
    int ret = 0;
    int num = mont->N.top;

    if ((a->top + b->top) > 2 * num)
        return 0;

    BN_CTX_start(ctx);
    tmp = BN_CTX_get(ctx);
    
    if (tmp == NULL)
        goto err;

    udi_check_top(tmp);
    if (a == b) {
        if (!BN_sqr(tmp, a, ctx))
            goto err;
        
    } else {
        if (!BN_mul(tmp, a, b, ctx))
            goto err;
    }

    printBIGNUM("T = ", tmp, "\n");
    
    /* reduce from aRR to aR */
    if (!PM_from_montgomery(r, rq, tmp, mont, ctx))
        goto err;
    printBIGNUM("r = ", r, "\n");
    printBIGNUM("rq = ", rq, "\n");

    if ((a == b) && (aq == bq))
    {
      if (!udi_mul_mont_fixed_top(tmp, a, aq, mont, ctx))
          goto err;
      printBIGNUM("tmp = ", tmp, "\n");
      if (!BN_mul_word(tmp, 2))
          goto err;
      
      if (!BN_add(rq, rq, tmp)) 
          goto err;

    } else {

      if (!udi_mul_mont_fixed_top(tmp, a, bq, mont, ctx))
          goto err;
      if (!BN_add(rq, rq, tmp)) 
          goto err;
      if (!udi_mul_mont_fixed_top(tmp, b, aq, mont, ctx))
          goto err;
      if (!BN_add(rq, rq, tmp)) 
          goto err;
    }
        
    ret = 1;
 err:
    BN_CTX_end(ctx);
    return ret;
}

int PM_from_montgomery(BIGNUM *ret, BIGNUM *retq, const BIGNUM *a, BN_MONT_CTX *mont, BN_CTX *ctx)
{
    int retn;

    retn = pm_from_mont_fixed_top(ret, retq, a, mont, ctx);
    udi_correct_top(ret);
    udi_correct_top(retq);
    udi_check_top(ret);
    udi_check_top(retq);

    return retn;
}

int pm_from_mont_fixed_top(BIGNUM *ret, BIGNUM *retq, const BIGNUM *a, BN_MONT_CTX *mont, BN_CTX *ctx)
{
    int retn = 0;
    
    BIGNUM *t1, *t2, *t2q;

    BN_CTX_start(ctx);
    t1 = BN_CTX_get(ctx);
    t2 = BN_CTX_get(ctx);
    if (t2 == NULL)
        goto err;

    if (!BN_copy(retq, a))
      goto err;
    BN_mask_bits(retq, mont->ri);
    if (!BN_mul(retq, retq, &mont->Ni, ctx))
      goto err;
    if (!BN_rshift(retq, retq, mont->ri))
      goto err;
    BN_mask_bits(retq, mont->ri);

    if (!BN_copy(t1, a))
        goto err;
    BN_mask_bits(t1, mont->ri);

    if (!BN_mul(t2, t1, &mont->Ni, ctx))
        goto err;
    BN_mask_bits(t2, mont->ri);
    
    if (!BN_mul(t1, t2, &mont->N, ctx))
        goto err;
    if (!BN_add(t2, a, t1))
        goto err;
    if (!BN_rshift(ret, t2, mont->ri))
        goto err;

    if (BN_ucmp(ret, &(mont->N)) >= 0) {
        if (!BN_usub(ret, ret, &(mont->N)))
            goto err;
        if (!BN_add_word(retq, 1))
            goto err;
    }
    retn = 1;
    udi_check_top(ret);
    udi_check_top(retq);
 err:
    BN_CTX_end(ctx);

    return retn;
}

/***
 * 
 *  Testing and timing 
 * 
 */

void time_mod_mult(uint64_t reps, BIGNUM *modulus)
{ 
  BN_CTX* bn_ctx = BN_CTX_new();
  BN_CTX_start(bn_ctx);

  BIGNUM *base = BN_CTX_get(bn_ctx);
  BIGNUM *res = BN_CTX_get(bn_ctx);
  BN_rand_range(base, modulus);
  BN_rand_range(res, modulus);

  start = clock();

  for (uint64_t i = 0; i < reps; ++i)
  {
    BN_mod_mul(res, res, base, modulus, bn_ctx);
  }

  diff = clock() - start;
  single_ms = ((double) diff * 1000/ CLOCKS_PER_SEC) / reps;
  printf("mul[%d] = %ld\n", BN_num_bits(modulus), diff * 1000/ CLOCKS_PER_SEC);

  BN_CTX_end(bn_ctx);
  BN_CTX_free(bn_ctx);
}

void test_mod_mult_same(uint64_t reps, BIGNUM *N)
{ 
  // BN_CTX* bn_ctx = BN_CTX_new();
  // BN_CTX* bn_udi_ctx = BN_CTX_new();
  // BN_CTX_start(bn_ctx);
  // BN_CTX_start(bn_udi_ctx);

  // BN_MONT_CTX *mnt_ctx = BN_MONT_CTX_new();
  // BN_MONT_CTX *mnt_udi_ctx = BN_MONT_CTX_new();
  // BN_MONT_CTX_set(mnt_ctx, N, bn_ctx);
  // UDI_MONT_CTX_set(mnt_udi_ctx, modulus, bn_udi_ctx);

  // BIGNUM *base = BN_CTX_get(bn_ctx);
  // BIGNUM *base_udi = BN_CTX_get(bn_udi_ctx);

  // BN_rand_range(base, modulus);

  // start = clock();

  // for (uint64_t i = 0; i < reps; ++i)
  // {
  //   BN_mod_mul_montgomery(base, base, base, mnt_ctx, bn_ctx);
  //   UDI_mod_mul_montgomery(base_udi, base_udi, base_udi, mnt_udi_ctx, bn_udi_ctx);
  // }
  
  // diff = clock() - start;
  // single_ms = ((double) diff * 1000/ CLOCKS_PER_SEC) / reps;
  // printf("Montgomery mod mult (%d-bit modulus)\n%lu repetitions, time: %lu msec, avg: %f msec\n",  BN_num_bits(modulus), reps, diff * 1000/ CLOCKS_PER_SEC, single_ms);

  // BN_CTX_end(bn_ctx);
  // BN_CTX_free(bn_ctx);
  // BN_MONT_CTX_free(mnt_ctx);
}


void time_mod_mult_mont(uint64_t reps, BIGNUM *modulus)
{ 
  BN_CTX* bn_ctx = BN_CTX_new();
  BN_CTX* bn_udi_ctx = BN_CTX_new();
  BN_CTX_start(bn_ctx);
  BN_CTX_start(bn_udi_ctx);

  BN_MONT_CTX *mnt_ctx = BN_MONT_CTX_new();
  BN_MONT_CTX *mnt_udi_ctx = BN_MONT_CTX_new();
  BN_MONT_CTX_set(mnt_ctx, modulus, bn_ctx);
  UDI_MONT_CTX_set(mnt_udi_ctx, modulus, bn_udi_ctx);

  BIGNUM *base = BN_CTX_get(bn_ctx);
  BIGNUM *base_udi = BN_CTX_get(bn_udi_ctx);

  BN_rand_range(base, modulus);

  start = clock();

  for (uint64_t i = 0; i < reps; ++i)
  {
    BN_mod_mul_montgomery(base, base, base, mnt_ctx, bn_ctx);
    UDI_mod_mul_montgomery(base_udi, base_udi, base_udi, mnt_udi_ctx, bn_udi_ctx);
  }
  
  diff = clock() - start;
  single_ms = ((double) diff * 1000/ CLOCKS_PER_SEC) / reps;
  printf("Montgomery mod mult (%d-bit modulus)\n%lu repetitions, time: %lu msec, avg: %f msec\n",  BN_num_bits(modulus), reps, diff * 1000/ CLOCKS_PER_SEC, single_ms);

  BN_CTX_end(bn_ctx);
  BN_CTX_free(bn_ctx);
  BN_MONT_CTX_free(mnt_ctx);
}

void time_mod_mult_udi(uint64_t reps, BIGNUM *modulus)
{ 
  BN_CTX* bn_ctx = BN_CTX_new();
  BN_CTX_start(bn_ctx);

  BN_MONT_CTX *mnt_ctx = BN_MONT_CTX_new();
  UDI_MONT_CTX_set(mnt_ctx, modulus, bn_ctx);

  BIGNUM *base = BN_CTX_get(bn_ctx);

  BN_rand_range(base, modulus);

  start = clock();

  for (uint64_t i = 0; i < reps; ++i)
  {
    UDI_mod_mul_montgomery(base, base, base, mnt_ctx, bn_ctx);
  }
  
  diff = clock() - start;
  single_ms = ((double) diff * 1000/ CLOCKS_PER_SEC) / reps;
  printf("Udi mod mult (%d-bit modulus)\n%lu repetitions, time: %lu msec, avg: %f msec\n",  BN_num_bits(modulus), reps, diff * 1000/ CLOCKS_PER_SEC, single_ms);

  BN_CTX_end(bn_ctx);
  BN_CTX_free(bn_ctx);
  BN_MONT_CTX_free(mnt_ctx);
}


void time_mod_exp_mont(uint64_t reps, BIGNUM *modulus, BIGNUM *exp)
{ 
  BN_CTX* bn_ctx = BN_CTX_new();
  BN_CTX_start(bn_ctx);

  BN_MONT_CTX *mnt_ctx = BN_MONT_CTX_new();
  BN_MONT_CTX_set(mnt_ctx, modulus, bn_ctx);

  BIGNUM *base = BN_CTX_get(bn_ctx);
  BIGNUM *res = BN_CTX_get(bn_ctx);

  BN_rand_range(base, modulus);
  
  start = clock();

  for (uint64_t i = 0; i < reps; ++i)
  {
    BN_mod_exp_mont(res, res, exp, modulus, bn_ctx, mnt_ctx);
  }

  diff = clock() - start;
  single_ms = ((double) diff * 1000/ CLOCKS_PER_SEC) / reps;
  
  printf("Montgomery mod exp (%d-bit modulus)\n%lu repetitions, time: %lu msec, avg: %f msec\n",  BN_num_bits(modulus), reps, diff * 1000/ CLOCKS_PER_SEC, single_ms);

  BN_CTX_end(bn_ctx);
  BN_CTX_free(bn_ctx);
  BN_MONT_CTX_free(mnt_ctx);
}

int main() {
  
  debug_print_off();
  
  BN_CTX *bn_ctx = BN_CTX_new();
  BN_CTX_start(bn_ctx);

  BIGNUM *P   = BN_CTX_get(bn_ctx);
  BIGNUM *Q   = BN_CTX_get(bn_ctx);
  BIGNUM *N   = BN_CTX_get(bn_ctx);
  BIGNUM *N2  = BN_CTX_get(bn_ctx);
  BIGNUM *exp = BN_CTX_get(bn_ctx);
  
  uint64_t prime_bits = 1024;

  BIGNUM *a   = BN_CTX_get(bn_ctx);
  BIGNUM *aq  = BN_CTX_get(bn_ctx);
  BIGNUM *b   = BN_CTX_get(bn_ctx);
  BIGNUM *bq  = BN_CTX_get(bn_ctx);
  BIGNUM *c   = BN_CTX_get(bn_ctx);
  BIGNUM *cq  = BN_CTX_get(bn_ctx);

  srand(time(NULL));
  unsigned long long_N = 5;
  unsigned long long_A = rand() % long_N;
  unsigned long long_Aq = rand() % long_N;
  unsigned long long_B = rand() % long_N;
  unsigned long long_Bq = rand() % long_N;
  unsigned long long_C = 0;
  unsigned long long_Cq = 0;

  BN_set_word(N, long_N);
  BN_set_word(a, long_A);
  BN_set_word(aq, long_Aq);
  BN_set_word(b, long_B);
  BN_set_word(bq, long_Bq);

  printBIGNUM("N  = ", N, "\n");
  printBIGNUM("a  = ", a, "\n");
  printBIGNUM("aq = ", aq, "\n");
  printBIGNUM("b  = ", b, "\n");
  printBIGNUM("bq = ", bq, "\n");

  BN_MONT_CTX *mnt_ctx = BN_MONT_CTX_new();
  UDI_MONT_CTX_set(mnt_ctx, N, bn_ctx);
  printf("mnt->ri = %d\n", mnt_ctx->ri);
  printBIGNUM("mnt->Ni = ", &mnt_ctx->Ni, "\n");
  
  UDI_to_montgomery(a, a, mnt_ctx, bn_ctx);
  UDI_to_montgomery(aq, aq, mnt_ctx, bn_ctx);
  UDI_to_montgomery(b, b, mnt_ctx, bn_ctx);
  UDI_to_montgomery(bq, bq, mnt_ctx, bn_ctx);
  printBIGNUM("a  = ", a, "\n");
  printBIGNUM("aq = ", aq, "\n");
  printBIGNUM("b  = ", b, "\n");
  printBIGNUM("bq = ", bq, "\n");

  PM_mod_mul_montgomery(c, cq, a, aq, a, aq, mnt_ctx, N, bn_ctx );

  UDI_from_montgomery(c, c, mnt_ctx, bn_ctx);
  UDI_from_montgomery(cq, cq, mnt_ctx, bn_ctx);

  printBIGNUM("c  = ", c, "\n");
  printBIGNUM("cq = ", cq, "\n");

  long_C = (long_A * long_A) % long_N;
  long_Cq = ((unsigned long) ((long_A + long_N*long_Aq) * (long_A + long_N*long_Aq)) / long_N) % long_N;
  printf("%d %d\n", BN_is_word(c, long_C), BN_is_word(cq, long_Cq));

  return 0;
  start = clock();

  BN_generate_prime_ex(P, prime_bits, 0, NULL, NULL, NULL);
  BN_generate_prime_ex(Q, prime_bits, 0, NULL, NULL, NULL);
  BN_mul(N, P, Q, bn_ctx);
  BN_mul(N2, N, N, bn_ctx);

  diff = clock() - start;
  single_ms = (double) diff * 1000/ CLOCKS_PER_SEC;

  printf("generating two safe primes P,Q (%ld-bits): %f msec\n", prime_bits, single_ms);
  printBIGNUM("P: ", P, "\n");
  printBIGNUM("Q: ", Q, "\n");
  printBIGNUM("N: ", N, "\n");
  printBIGNUM("N2: ", N2, "\n");

  #define MUL_NUM 100000

  test_mod_mult_same(10, N);
  time_mod_mult_mont(MUL_NUM, N2);
  time_mod_mult_udi(MUL_NUM, N2);

  #define EXP_NUM 100
  
  debug_print_on();
  //BN_set_flags(N2, BN_FLG_CONSTTIME);
  time_mod_exp_mont(EXP_NUM, N2, N);

  BN_mod_exp(N, N, N, N, bn_ctx);
  
  BN_CTX_end(bn_ctx);
  BN_CTX_free(bn_ctx);
}