#include <string.h>
#include <assert.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>

#include <openssl/bn.h>
#include <openssl/ec.h>
#include <openssl/sha.h>
#include <openssl/err.h>
#include <openssl/objects.h>

typedef EC_POINT *group_el_t;
typedef BIGNUM *scalar_t;

clock_t start;
clock_t diff;

FILE *out_file = NULL;
double base_time = 0;
double single_ms;

void print_to_file(uint64_t exp_bits, uint64_t mod_bits, double time_ms)
{
  if (!out_file) return;

  if (base_time == 0) {
    base_time = time_ms;
    fprintf(out_file, "base_time_ms == %f\n", base_time);
  }

  fprintf(out_file, "if (exp, mod) == (%ld,%ld) : return %f * base_time_ms\n", exp_bits, mod_bits, time_ms/base_time);
}

void set_file_letter(char l)
{
  if (!out_file) return;

  base_time = 0;
  fprintf(out_file, "%c\n", l);
}

void printHexBytes_padded(const char * prefix, const uint8_t *src, unsigned len, const char * suffix) {
  if (len == 0) {
    printf("%s <0 len char array> %s", prefix, suffix);
    return;
  }

  printf("%s", prefix);
  unsigned int i;
  for (i = 0; i < len-1; ++i) {
    printf("%02x",src[i] & 0xff);
  }
  printf("%02x%s",src[i] & 0xff, suffix);
}

void sample_in_range(const scalar_t range_mod, scalar_t rnd, int coprime)
{
  BN_rand_range(rnd, range_mod);

  if (coprime)
  { 
    BIGNUM *gcd = BN_new();
    BN_CTX* bn_ctx = BN_CTX_new();
    BN_gcd(gcd, range_mod, rnd, bn_ctx);
    
    while (!BN_is_one(gcd))
    {
      BN_rand_range(rnd, range_mod);
      BN_gcd(gcd, range_mod, rnd, bn_ctx);
    }
    BN_clear_free(gcd);
    BN_CTX_free(bn_ctx);
  }
}

void time_sampling_scalars(uint64_t reps, const scalar_t range, int coprime)
{
  scalar_t alpha;

  start = clock();

  for (uint64_t i = 0; i < reps; ++i)
  {
    alpha = BN_new();
    sample_in_range(range, alpha, coprime);
    //printf("alpha: %s\n", BN_bn2dec(alphas[i]));
    BN_free(alpha);
  }

  diff = clock() - start;
  single_ms = ((double) diff * 1000/ CLOCKS_PER_SEC) / reps;
  printf("sampling scalars (coprime: %d)\n%lu repetitions, time: %lu msec, avg: %f msec\n", coprime, reps, diff * 1000/ CLOCKS_PER_SEC, single_ms);
}

void time_mod_mult(uint64_t reps, scalar_t modulus)
{ 
  BN_CTX* bn_ctx = BN_CTX_new();
  
  scalar_t base = BN_new();
  sample_in_range(modulus, base, 0);
  
  start = clock();

  for (uint64_t i = 0; i < reps; ++i)
  {
    BN_mod_mul(base, base, base, modulus, bn_ctx);
    //printf("%s\n", BN_bn2dec(base));
  }

  diff = clock() - start;
  single_ms = ((double) diff * 1000/ CLOCKS_PER_SEC) / reps;
  printf("computing mod mult (%d-bit modulus)\n%lu repetitions, time: %lu msec, avg: %f msec\n", BN_num_bits(modulus), reps, diff * 1000/ CLOCKS_PER_SEC, single_ms);

  BN_free(base);
  BN_CTX_free(bn_ctx);
}


void time_mod_mult_mont(uint64_t reps, scalar_t modulus)
{ 
  BN_CTX* bn_ctx = BN_CTX_new();

  BN_MONT_CTX *mnt_ctx = BN_MONT_CTX_new();
  BN_MONT_CTX_set(mnt_ctx, modulus, bn_ctx);

  scalar_t base = BN_new();

  sample_in_range(modulus, base, 0);

  start = clock();

  for (uint64_t i = 0; i < reps; ++i)
  {
    BN_mod_exp_mont(base, base, base, modulus, bn_ctx, mnt_ctx);
    //printf("%s\n", BN_bn2dec(base));
  }
  
  diff = clock() - start;
  single_ms = ((double) diff * 1000/ CLOCKS_PER_SEC) / reps;
  printf("computing mod mult mont (%d-bit modulus)\n%lu repetitions, time: %lu msec, avg: %f msec\n",  BN_num_bits(modulus), reps, diff * 1000/ CLOCKS_PER_SEC, single_ms);

  BN_free(base);
  BN_CTX_free(bn_ctx);
  BN_MONT_CTX_free(mnt_ctx);
}

void time_mod_exp(uint64_t reps, scalar_t modulus, scalar_t exp)
{ 
  BN_CTX* bn_ctx = BN_CTX_new();
  scalar_t base = BN_new();

  sample_in_range(modulus, base, 0);
  
  start = clock();

  for (uint64_t i = 0; i < reps; ++i)
  {
    BN_mod_exp(base, base, exp, modulus, bn_ctx);
  }

  diff = clock() - start;
  single_ms = ((double) diff * 1000/ CLOCKS_PER_SEC) / reps;
  print_to_file(BN_num_bits(exp), BN_num_bits(modulus), single_ms);
  printf("computing exp[%d] mod [%d-bits]\n%lu repetitions, time: %lu msec, avg: %f msec\n",  BN_num_bits(exp), BN_num_bits(modulus), reps, diff * 1000/ CLOCKS_PER_SEC, single_ms);


  BN_free(base);
  BN_CTX_free(bn_ctx);
}

void time_mod_exp_mont(uint64_t reps, scalar_t modulus, scalar_t exp)
{ 
  BN_CTX* bn_ctx = BN_CTX_new();

  BN_MONT_CTX *mnt_ctx = BN_MONT_CTX_new();
  BN_MONT_CTX_set(mnt_ctx, modulus, bn_ctx);

  scalar_t base = BN_new();
  sample_in_range(modulus, base, 0);

  start = clock();

  for (uint64_t i = 0; i < reps; ++i)
  {
    BN_mod_exp_mont(base, base, exp, modulus, bn_ctx, mnt_ctx);
    //printf("%s\n", BN_bn2dec(base));
  }

  diff = clock() - start;
  single_ms = ((double) diff * 1000/ CLOCKS_PER_SEC) / reps;
  print_to_file(BN_num_bits(exp), BN_num_bits(modulus), single_ms);
  printf("computing Montgomery mod exp (%d-bit modulus)\n%lu repetitions, time: %lu msec, avg: %f msec\n",  BN_num_bits(modulus), reps, diff * 1000/ CLOCKS_PER_SEC, single_ms);

  BN_free(base);
  BN_CTX_free(bn_ctx);
  BN_MONT_CTX_free(mnt_ctx);
}

void time_naive_mod_exp(uint64_t reps, scalar_t modulus, scalar_t exp)
{
  scalar_t base = BN_new();
  sample_in_range(modulus, base, 0);

  BN_CTX* bn_ctx = BN_CTX_new();
  scalar_t res = BN_new();
  
  int bit_len_exp = BN_num_bits(exp);

  start = clock();

  for (uint64_t i = 0; i < reps; ++i)
  {
    BN_set_word(res, 1);

    for (int i = bit_len_exp-1; i >= 0; --i)
    {
      BN_mod_sqr(res, res, modulus, bn_ctx);

      if (BN_is_bit_set(exp, i))
        BN_mod_mul(res, res, base, modulus, bn_ctx);
      //printf("%d %s\n", BN_is_bit_set(exp, i), BN_bn2dec(res));
    }
    
    //if (i <= 2) printf("pow(%s, %s, %s)\n%s\n", BN_bn2dec(base), BN_bn2dec(exp), BN_bn2dec(modulus), BN_bn2dec(res));

    BN_copy(base,res);
  }

  diff = clock() - start;
  single_ms = ((double) diff * 1000/ CLOCKS_PER_SEC) / reps;

  printf("computing naive mod exp (%d-bit modulus)\n%ld repetitions, time: %ld msec, avg: %f msec\n",  BN_num_bits(modulus), reps, diff * 1000/ CLOCKS_PER_SEC, single_ms);

  BN_free(res);
  BN_free(base);
  BN_CTX_free(bn_ctx);
}

void time_naive_mod_exp_mont(uint64_t reps, scalar_t modulus, scalar_t exp)
{
  scalar_t base = BN_new();

  sample_in_range(modulus, base, 0);

  BN_CTX* bn_ctx = BN_CTX_new();
  BN_MONT_CTX *mnt_ctx = BN_MONT_CTX_new();
  BN_MONT_CTX_set(mnt_ctx, modulus, bn_ctx);
  
  scalar_t res = BN_new();
  
  int bit_len_exp = BN_num_bits(exp);

  start = clock();

  for (uint64_t i = 0; i < reps; ++i)
  {
    BN_set_word(res, 1);

    for (int i = bit_len_exp-1; i >= 0; --i)
    {
      BN_mod_mul_montgomery(res, res, res, mnt_ctx, bn_ctx);

      if (BN_is_bit_set(exp, i))
        BN_mod_mul_montgomery(res, res, base, mnt_ctx, bn_ctx);
      //printf("%d %s\n", BN_is_bit_set(exp, i), BN_bn2dec(res));
    }
    
    //if (i <= 2) printf("pow(%s, %s, %s)\n%s\n", BN_bn2dec(base), BN_bn2dec(exp), BN_bn2dec(modulus), BN_bn2dec(res));

    BN_copy(base,res);
  }

  diff = clock() - start;
  single_ms = ((double) diff * 1000/ CLOCKS_PER_SEC) / reps;

  printf("computing naive mod exp Montgeomery (%d-bit modulus)\n%ld repetitions, time: %ld msec, avg: %f msec\n",  BN_num_bits(modulus), reps, diff * 1000/ CLOCKS_PER_SEC, single_ms);

  BN_free(res);
  BN_free(base);
  BN_CTX_free(bn_ctx);
}

void time_parralel_mod_exp(uint64_t reps, scalar_t modulus, scalar_t exp)
{
  uint64_t i;
  scalar_t* v_base = calloc(reps, sizeof(scalar_t));
  scalar_t* v_resu = calloc(reps, sizeof(scalar_t));

  for (i = 0; i < reps; ++i)
  {
    v_base[i] = BN_new();
    sample_in_range(modulus, v_base[i], 0);

    v_resu[i] = BN_new();
    BN_set_word(v_resu[i], 1);
  }

  BN_CTX* bn_ctx = BN_CTX_new();
  
  int bit_len_exp = BN_num_bits(exp);

  start = clock();

  for (int bi = bit_len_exp-1; bi >= 0; --bi)
  {
    for (i = 0; i < reps; ++i) BN_mod_sqr(v_resu[i], v_resu[i], modulus, bn_ctx);

    if (BN_is_bit_set(exp, bi))
      for (i = 0; i < reps; ++i) BN_mod_mul(v_resu[i], v_resu[i], v_base[i], modulus, bn_ctx);
    //printf("%d %s\n", BN_is_bit_set(exp, i), BN_bn2dec(res));
  }
  
  //if (i <= 2) printf("pow(%s, %s, %s)\n%s\n", BN_bn2dec(base), BN_bn2dec(exp), BN_bn2dec(modulus), BN_bn2dec(res));

  diff = clock() - start;
  single_ms = ((double) diff * 1000/ CLOCKS_PER_SEC) / reps;

  printf("computing parralel mod exp (%d-bit modulus)\n%ld repetitions, time: %ld msec, avg: %f msec\n",  BN_num_bits(modulus), reps, diff * 1000/ CLOCKS_PER_SEC, single_ms);

  BN_CTX_free(bn_ctx);
  for (i = 0; i < reps; ++i)
  {
    BN_free(v_resu[i]);
    BN_free(v_base[i]);
  }
  free(v_resu);
  free(v_base);
}

void time_hashing(uint64_t reps, const uint8_t* data, uint64_t data_len)
{ 
  unsigned char digest[64];

  start = clock();

  for (uint64_t i = 0; i < reps; ++i)
  {
    SHA512_CTX sha_ctx;
    SHA512_Init(&sha_ctx);
    SHA512_Update(&sha_ctx, data, data_len);
    SHA512_Final(digest, &sha_ctx);
    
    memcpy((uint8_t*) data, digest, data_len < 64 ? data_len : 64);
    //printHexBytes_padded("", digest, sizeof(digest), "\n");
  }

  diff = clock() - start;
  single_ms = ((double) diff * 1000/ CLOCKS_PER_SEC) / reps;

  printf("Sha512 Digest\n%lu repetitions, time: %lu msec, avg: %f msec\n", reps, diff * 1000/ CLOCKS_PER_SEC, single_ms);
}

void time_ec_exp(uint64_t reps, scalar_t exp)
{ 
  BN_CTX* bn_ctx = BN_CTX_new();

  const EC_GROUP* ec_ctx = EC_GROUP_new_by_curve_name(NID_secp256k1);
  group_el_t base = EC_POINT_dup(EC_GROUP_get0_generator(ec_ctx), ec_ctx);

  scalar_t mod_exp = BN_new();
  BN_mod(mod_exp, exp, EC_GROUP_get0_order(ec_ctx), bn_ctx);
  
  start = clock();

  for (uint64_t i = 0; i < reps; ++i)
  {
    EC_POINT_mul(ec_ctx, base, NULL, base, mod_exp, bn_ctx);
    //BN_mul(mod_exp, mod_exp, mod_exp, bn_ctx);
    //printf("%s\n", EC_POINT_point2hex(ec_ctx, base, POINT_CONVERSION_UNCOMPRESSED, bn_ctx));
  }

  diff = clock() - start;
  single_ms = ((double) diff * 1000/ CLOCKS_PER_SEC) / reps;
  printf("computing Elliptic curve exp (%d-bit exp)\n%lu repetitions, time: %lu msec, avg: %f msec\n",  BN_num_bits(mod_exp), reps, diff * 1000/ CLOCKS_PER_SEC, single_ms);

  EC_POINT_free(base);
  BN_free(mod_exp);
  BN_CTX_free(bn_ctx);
}

int main()
{
  out_file = fopen("./timing.txt","w");

  if(out_file == NULL)
  {
    printf("Can't wrtie to file\n");   
  }

  BN_CTX* bn_ctx = BN_CTX_new();

  uint64_t safe_prime_bits = 256;

  scalar_t P = BN_new();
  scalar_t Q = BN_new();

  start = clock();

  BN_generate_prime_ex(P, safe_prime_bits, 1, NULL, NULL, NULL);
  BN_generate_prime_ex(Q, safe_prime_bits, 1, NULL, NULL, NULL);

  diff = clock() - start;
  printf("generating two safe primes (%d-bits): %lu msec\nP: %s\nQ: %s\n", BN_num_bits(P), diff * 1000/ CLOCKS_PER_SEC, BN_bn2dec(P), BN_bn2dec(Q));

  scalar_t N = BN_new();
  scalar_t N2 = BN_new();

  BN_mul(N, P, Q, bn_ctx);
  BN_mul(N2, N, N, bn_ctx);
  printf("N (%d-bits):\t%s\nN2 (%d-bits):\t%s\n", BN_num_bits(N), BN_bn2dec(N), BN_num_bits(N2), BN_bn2dec(N2));

  time_ec_exp(1000, N);

  // time_mod_mult(1000, P);
  // time_mod_mult(1000, N);
  // time_mod_mult(1000, N2);
  
  set_file_letter('E');

  scalar_t exp = BN_new();

  BN_rand(exp, safe_prime_bits, 1, 0);
  time_mod_exp(1000, N2, exp);
  time_mod_exp(1000, P, exp);
  
  BN_rand(exp, 2*safe_prime_bits, 1, 0);  
  time_mod_exp(1000, P,  exp);
  time_mod_exp(1000, N, P);
  time_mod_exp(1000, N2, exp);
  
  BN_rand(exp, 4*safe_prime_bits, 1, 0);
  time_mod_exp(1000, N2, exp);
  
  BN_rand(exp, safe_prime_bits/2, 1, 0);
  time_mod_exp(1000, N2, exp);
  time_mod_exp(1000, N, exp);
  time_mod_exp(1000, P, exp);

  BN_rand(exp, 3*safe_prime_bits/2, 1, 0);
  time_mod_exp(1000, N2, exp);

  set_file_letter('F');

  scalar_t base = BN_new();
  uint64_t max_prime_bits = 2048;
  

  for (uint64_t prime_bits = 32; prime_bits <= max_prime_bits; prime_bits *= 2) 
  {
    BN_generate_prime_ex(P, prime_bits, 0, NULL, NULL, NULL);
    for (uint64_t i = 2; i <= prime_bits-3; i *= 2) {
      BN_rand(exp, i, 1, 0);
      time_mod_exp(1000, P, exp);
    }
  }
  
  // time_mod_exp_mont(1000, P,  exp);
  // time_mod_exp_mont(1000, N,  exp);
  // time_mod_exp_mont(1000, N2, exp);

  // time_naive_mod_exp(1000, P, exp);
  // time_naive_mod_exp(1000, N, exp);
  // time_naive_mod_exp(1000, N2, exp);

  //time_naive_mod_exp_mont(1000, P, exp);

  // time_parralel_mod_exp(1000, P, exp);
  // time_parralel_mod_exp(1000, N, exp);
  // time_parralel_mod_exp(1000, N2, exp);


  // uint8_t data[1000] = {7};
  // printf("%ld\n", sizeof(data));
  // time_hashing(100000, data, sizeof(data));

  BN_free(exp);
  BN_free(base);
  BN_free(P);
  BN_free(Q);
  BN_free(N);
  BN_free(N2);

  BN_CTX_free(bn_ctx);
}

/*
int main() {
  BN_CTX* ctx = BN_CTX_new();

  scalar_t num = BN_new();

  BN_generate_prime_ex(num, 100, 1, NULL, NULL, NULL);
  
  printf("num (%d-bits)\nN: %s\n", BN_num_bits(num), BN_bn2dec(num));

  BN_free(num);
  BN_CTX_free(ctx);
}
*/