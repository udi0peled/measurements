#include <string.h>
#include <assert.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>

#include <openssl/bn.h>
#include <openssl/ec.h>
#include <openssl/sha.h>
#include <openssl/err.h>
#include <openssl/objects.h>

typedef BIGNUM* scalar_t;
clock_t start;
clock_t diff;
double  single_ms;

char *res_str;
char temp_str[1500];

void append_str(const char str[]) {
    size_t append_len = strlen(str) + 1;
    size_t res_len = strlen(res_str);
    res_str = realloc(res_str, res_len + append_len);
    memcpy(res_str + res_len, str, append_len);
}

void my_printf(const char *format, ...) {
    va_list args;
    va_start(args, format);
    //vsprintf(temp_str, format, args);
    append_str(temp_str);
    va_end(args);

    va_start(args, format);
    vprintf(format, args);
    va_end(args);
}

void printHexBytes(const char * prefix, const uint8_t *src, unsigned len, const char * suffix, int print_len) {
  if (len == 0) {
    my_printf("%s <0 len char array> %s", prefix, suffix);
    return;
  }

  if (print_len) {
      my_printf("[%u]", len);
  }
  my_printf("%s", prefix);
  unsigned int i;
  for (i = 0; i < len-1; ++i) {
    my_printf("%02x",src[i] & 0xff);
  }
  my_printf("%02x%s",src[i] & 0xff, suffix);
}

void printBIGNUM(const char * prefix, const BIGNUM *bn, const char * suffix) {
  char *bn_str = BN_bn2dec(bn);
  my_printf("%s%s%s", prefix, bn_str, suffix);
  free(bn_str);
}

void sample_in_range(const scalar_t range_mod, scalar_t rnd, int coprime)
{  
    BN_CTX *bn_ctx = BN_CTX_secure_new();
    BN_CTX_start(bn_ctx);

    BN_rand_range(rnd, range_mod);
    if (coprime)
    { 
        BIGNUM *gcd = BN_CTX_get(bn_ctx);
        BN_gcd(gcd, range_mod, rnd, bn_ctx);
        
        while (!BN_is_one(gcd))
        {
        BN_rand_range(rnd, range_mod);
        BN_gcd(gcd, range_mod, rnd, bn_ctx);
        }
    }

    BN_CTX_end(bn_ctx);
    BN_CTX_free(bn_ctx);
}

// Multiplications

void time_binary_field_mul(uint64_t reps, int mod_poly[])
{
    BN_CTX *bn_ctx = BN_CTX_secure_new();
    BN_CTX_start(bn_ctx);

    BIGNUM *f = BN_CTX_get(bn_ctx);
    BIGNUM *a = BN_CTX_get(bn_ctx);
    BIGNUM *b = BN_CTX_get(bn_ctx);

    BN_GF2m_arr2poly(mod_poly, f);
    sample_in_range(f, a, 1);
    sample_in_range(f, b, 1);

    start = clock();

    for (uint64_t i = 0; i < reps; ++i) {
        BN_GF2m_mod_mul_arr(a, a, b, mod_poly, bn_ctx);
    }

    diff = clock() - start;
    single_ms = ((double) diff * 1000/ CLOCKS_PER_SEC)/reps;

    my_printf("computing binary mod mult (%d-bit modulus)\n%lu repetitions, time: %lu msec, avg: %f msec\n", BN_num_bits(f)-1, reps, diff * 1000/ CLOCKS_PER_SEC, single_ms); 

    BN_CTX_end(bn_ctx);
    BN_CTX_free(bn_ctx);
}

void time_binary_field_inv_add(uint64_t reps, int mod_poly[])
{
    BN_CTX *bn_ctx = BN_CTX_secure_new();
    BN_CTX_start(bn_ctx);

    BIGNUM *f = BN_CTX_get(bn_ctx);
    BIGNUM *a = BN_CTX_get(bn_ctx);
    BIGNUM *b = BN_CTX_get(bn_ctx);

    BN_GF2m_arr2poly(mod_poly, f);
    sample_in_range(f, a, 1);
    sample_in_range(f, b, 1);

    start = clock();

    for (uint64_t i = 0; i < reps; ++i) {
        BN_GF2m_mod_inv(a, a, f, bn_ctx);
        BN_GF2m_add(a, a, b);
    }

    diff = clock() - start;
    single_ms = ((double) diff * 1000/ CLOCKS_PER_SEC)/reps;

    my_printf("computing binary mod inv+add (%d-bit modulus)\n%lu repetitions, time: %lu msec, avg: %f msec\n", BN_num_bits(f), reps, diff * 1000/ CLOCKS_PER_SEC, single_ms); 

    BN_CTX_end(bn_ctx);
    BN_CTX_free(bn_ctx);
}


void time_mod_mult(uint64_t reps, scalar_t modulus)
{ 
  BN_CTX* bn_ctx = BN_CTX_secure_new();
  BN_CTX_start(bn_ctx);

  scalar_t base = BN_CTX_get(bn_ctx);
  scalar_t res = BN_CTX_get(bn_ctx);
  sample_in_range(modulus, base, 0);
  sample_in_range(modulus, res, 0);

  start = clock();

  for (uint64_t i = 0; i < reps; ++i)
  {
    BN_mod_mul(res, res, base, modulus, bn_ctx);
    //my_printf("%s\n", BN_bn2dec(base));
  }

  diff = clock() - start;
  single_ms = ((double) diff * 1000/ CLOCKS_PER_SEC) / reps;
  my_printf("computing mod mult (%d-bit modulus)\n%lu repetitions, time: %lu msec, avg: %f msec\n", BN_num_bits(modulus), reps, diff * 1000/ CLOCKS_PER_SEC, single_ms);

  BN_CTX_end(bn_ctx);
  BN_CTX_free(bn_ctx);
}

// Exponentiation

void time_mod_exp(uint64_t reps, scalar_t modulus, scalar_t exp)
{ 
  BN_CTX* bn_ctx = BN_CTX_secure_new();
  BN_CTX_start(bn_ctx);

  scalar_t base = BN_CTX_get(bn_ctx);
  scalar_t res = BN_CTX_get(bn_ctx);

  sample_in_range(modulus, base, 0);
  
  start = clock();

  for (uint64_t i = 0; i < reps; ++i)
  {
    BN_mod_exp(res, res, exp, modulus, bn_ctx);
  }

  diff = clock() - start;
  single_ms = ((double) diff * 1000/ CLOCKS_PER_SEC) / reps;
  
  my_printf("computing exp[%d] mod [%d-bits]\n%lu repetitions, time: %lu msec, avg: %f msec\n",  BN_num_bits(exp), BN_num_bits(modulus), reps, diff * 1000/ CLOCKS_PER_SEC, single_ms);

  BN_CTX_end(bn_ctx);
  BN_CTX_free(bn_ctx);
}

void time_mod_exp_mont(uint64_t reps, scalar_t modulus, scalar_t exp)
{ 
  BN_CTX* bn_ctx = BN_CTX_secure_new();
  BN_CTX_start(bn_ctx);

  BN_MONT_CTX *mnt_ctx = BN_MONT_CTX_new();
  BN_MONT_CTX_set(mnt_ctx, modulus, bn_ctx);

  scalar_t base = BN_CTX_get(bn_ctx);
  scalar_t res = BN_CTX_get(bn_ctx);

  sample_in_range(modulus, base, 0);
  
  start = clock();

  for (uint64_t i = 0; i < reps; ++i)
  {
  
    BN_mod_exp_mont(res, res, exp, modulus, bn_ctx, mnt_ctx);
    //my_printf("%s\n", BN_bn2dec(base));
  }

  diff = clock() - start;
  single_ms = ((double) diff * 1000/ CLOCKS_PER_SEC) / reps;
  
  my_printf("computing Montgomery mod exp (%d-bit modulus)\n%lu repetitions, time: %lu msec, avg: %f msec\n",  BN_num_bits(modulus), reps, diff * 1000/ CLOCKS_PER_SEC, single_ms);

  BN_CTX_end(bn_ctx);
  BN_CTX_free(bn_ctx);
  BN_MONT_CTX_free(mnt_ctx);
}

// Group 
void time_ec_exp(uint64_t reps, scalar_t exp)
{ 
  BN_CTX* bn_ctx = BN_CTX_new();

  const EC_GROUP* ec_ctx = EC_GROUP_new_by_curve_name(NID_secp256k1);
  EC_POINT *base = EC_POINT_dup(EC_GROUP_get0_generator(ec_ctx), ec_ctx);

  scalar_t mod_exp = BN_new();
  unsigned long int delta = 1123;
  BN_mod(mod_exp, exp, EC_GROUP_get0_order(ec_ctx), bn_ctx);
  start = clock();

  for (uint64_t i = 0; i < reps; ++i)
  {
    EC_POINT_mul(ec_ctx, base, NULL, base, mod_exp, bn_ctx);
    BN_add_word(mod_exp, delta);
    //printf("%s\n", EC_POINT_point2hex(ec_ctx, base, POINT_CONVERSION_UNCOMPRESSED, bn_ctx));
  }

  diff = clock() - start;
  single_ms = ((double) diff * 1000/ CLOCKS_PER_SEC) / reps;
  printf("computing Elliptic curve exp (%d-bit exp)\n%lu repetitions, time: %lu msec, avg: %f msec\n",  BN_num_bits(mod_exp), reps, diff * 1000/ CLOCKS_PER_SEC, single_ms);

  EC_POINT_free(base);
  BN_free(mod_exp);
  BN_CTX_free(bn_ctx);
}


//Returns length of string of measurements results
char *ipad_measurements() {
    BN_CTX *bn_ctx = BN_CTX_secure_new();
    res_str = malloc(0);

    int poly8192[] = {8192, 9 , 5 , 2 , 0, -1};
    int poly32[] = {32, 19, 14, 13, 0, -1};
    
    time_binary_field_mul(100000, poly8192);
    time_binary_field_inv_add(1000, poly8192);

    time_binary_field_mul(100000, poly32);
    time_binary_field_inv_add(1000, poly32);

    scalar_t P   = BN_new();
    scalar_t Q   = BN_new();
    scalar_t N   = BN_new();
    scalar_t N2  = BN_new();
    scalar_t exp = BN_new();
    
    BN_rand(exp, 256, 1, 0);
    time_ec_exp(1000, exp);
    
    uint64_t safe_prime_bits = 1024;

    BN_generate_prime_ex(P, safe_prime_bits, 1, NULL, NULL, NULL);
    BN_generate_prime_ex(Q, safe_prime_bits, 1, NULL, NULL, NULL);
    BN_mul(N, P, Q, bn_ctx);
    BN_mul(N2, N, N, bn_ctx);

    diff = clock() - start;
    single_ms = (double) diff * 1000/ CLOCKS_PER_SEC;

    my_printf("generating two safe primes P,Q (%ld-bits): %f msec\n", safe_prime_bits, single_ms);
    printBIGNUM("P: ", P, "\n");
    printBIGNUM("Q: ", Q, "\n");
    printBIGNUM("N: ", N, "\n");
    printBIGNUM("N2: ", N2, "\n");

    #define MULT_NUM 100000
    
    time_mod_mult(MULT_NUM, P);
    time_mod_mult(MULT_NUM, N);
    time_mod_mult(MULT_NUM, N2);

    #define EXP_NUM 1000

    time_mod_exp(EXP_NUM, N2, N);
    
    BN_rand(exp, safe_prime_bits, 1, 0);
    time_mod_exp(EXP_NUM, N2, exp);
    
    BN_rand(exp, safe_prime_bits/2, 1, 0);
    time_mod_exp(EXP_NUM, N2, exp);
    
    BN_rand(exp, safe_prime_bits/4, 1, 0);
    time_mod_exp(EXP_NUM, N2, exp);
    
    BN_rand(exp, 3*safe_prime_bits/4, 1, 0);
    time_mod_exp(EXP_NUM, N2, exp);
    
    BN_rand(exp, 5*safe_prime_bits/4, 1, 0);
    time_mod_exp(EXP_NUM, P, exp);

    BN_rand(exp, 7*safe_prime_bits/4, 1, 0);
    time_mod_exp(EXP_NUM, P, exp);

    BN_free(P);
    BN_free(Q);
    BN_free(N);
    BN_free(N2);
    BN_CTX_free(bn_ctx);
    return res_str;
}

int main()
{
  char *res = ipad_measurements();
  printf("%s", res);
  printf("%lu", strlen(res));
  free(res);
}