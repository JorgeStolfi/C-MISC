/* Finds primes between two 64-bit integers. */
/* Last edited on 2013-07-26 20:55:33 by stolfilocal */

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <assert.h>
#include <limits.h>

/* PROTOTYPES */

uint64_t find_primes(uint64_t a, uint64_t b);
  /* Prints to stdout all primes in {a..b}. Returns the cont of those primes. */ 

uint64_t *alloc_bits(uint64_t n);
  /* Allocates a vector of {uint64_t}s with {n} bits total, at least. Sets all bits to 1. */ 

uint64_t print_primes(uint64_t a, uint64_t v[], uint64_t n);
  /* Prints all integers {a + k} with {k} in {0..n-1} such that bit {k} of {v} is set.
    Returns the number of bits found. */ 

uint64_t isqrt(uint64_t x);
  /* Returns the square root of {x}, rounded down. */ 

#define get_bit(v,k)   ((1ull) & ((v)[(k)/64] >> ((k)%64)))
  /* Returns bit index {k} in {uint64_t} vector {v}. */

#define off_bit(v,k) ((v)[(k)/64]) &= ~((1ull) << ((k)%64));
  /* Clears bit index {k} in {uint64_t} vector {v}. */

#define debug 0
  /* Define as 1 to get a trace of the algorithm. */

/* IMPLEMENTATIONS */

int main(int argc, char** argv)
{
  fprintf(stderr, "parsing the arguments...\n");
  assert(argc == 1+2);
  char *end;
  uint64_t a = strtoull(argv[1], &end, 10); assert((*end) == 0);
  uint64_t b = strtoull(argv[2], &end, 10); assert((*end) == 0);
  uint64_t np = find_primes(a, b);
  fprintf(stderr, "found %lu primes in {%lu .. %lu}\n", np, a, b);
  return 0;
}

uint64_t find_primes(uint64_t a, uint64_t b)
{
  if (a < 2) { a = 2; }
  if (a > b) { return 0; }

  /* Erathostenes sieve for candidates in {a..b}, shifted by {a}: */
  uint64_t n = b - a + 1; /* Number of candidates to check. */
  fprintf(stderr, "will test %lu candidates in {%lu .. %lu}\n", n, a, b);
  uint64_t *hi = alloc_bits(n);
  
  /* Eratosthenes sieve for possible divisors in {2..sqrt(b)}: */
  uint64_t m = isqrt(b);
  fprintf(stderr, "possible prime divisors are in {2 .. %lu}\n", m);
  uint64_t *lo = alloc_bits(m+1);
  
  /* Cross out 0 and 1 from {0..m}, just for consistency: */
  off_bit(lo,0);
  off_bit(lo,1);
  uint64_t p = 2; /* Next candidate for prime. */
  while (1) 
    { /* Scan for next prime: */
      while ((p <= m) && (get_bit(lo,p) == 0)) { p++; }
      if (p > m) 
        { /* All composites in {a..b} have been cleared: */
          fprintf(stderr, "  writing the primes to stdout\n");
          return print_primes(a,hi,n);
        }

      if (debug) { fprintf(stderr, "crossing out all multiples of %lu:", p); }
      
      /* Cross out from {lo} all multiples of {p} in {0..m}, starting with {p^2}: */
      uint64_t k = p*p; /* Smallest multiple not yet crossed out. */
      while (k <= m) 
        { if (debug) { fprintf(stderr, " %lu", k); } 
          off_bit(lo,k);
          k += p;
        }
      
      /* Cross out from {hi} all multiples of {p} in {a..b}, starting with {p^2}: */
      uint64_t q = (a + p - 1)/p;
      if (q < p) { q = p; }
      k = q*p - a; /* Smallest multiple not yet crossed out, minus {a}. */
      while (k < n) 
        { if (debug) { fprintf(stderr, " %lu", k); }  
          off_bit(hi,k); 
          k += p;
        }
      
      if (debug) { fprintf(stderr, "\n"); } 
      /* Done with this divisor: */
      p++;
    }
}

uint64_t print_primes(uint64_t a, uint64_t v[], uint64_t n)
{
  uint64_t k;
  uint64_t np = 0;
  for (k = 0; k < n; k++) 
    { if (get_bit(v,k) != 0)
        { fprintf(stdout, "%22lu\n", a+k); np++; }
    }
  fflush(stdout);
  return np;
}

uint64_t *alloc_bits(uint64_t n)
{ 
  if (n == 0) { return NULL; }
  /* Compute the number {nW} of 64-bit words needed: */
  uint64_t nW = n / 64; 
  if (nW*64 < n) { nW ++;}
  assert(nW*64 >= n); /* Paranoia. */
  /* Compute approx number of bytes {nB} and gigabytes {nGB}: */
  uint64_t nB = nW*sizeof(uint64_t);
  uint64_t nGB = (nB + 500000000ull)/1000000000ull;
  fprintf(stderr, "trying to allocate a vector with %lu bits = %lu words ~= %lu GB\n", n, nW, nGB);
  
  /* Try to allocate: */
  size_t sz = (size_t)nB;
  assert(sz > 0);
  assert(nB == (uint64_t)sz);
  uint64_t *v = malloc(sz); assert(v != NULL);
  uint64_t k;
  for (k = 0; k < nW; k++) { v[k] = ~(0ull); }
  return v;
}

uint64_t isqrt(uint64_t x)
{ 
  if (x <= 0) { return 0; }
  if (x <= 3) { return 1; }
  uint64_t y = 2;
  do 
    { y = (y + x/y)/2;
      y = (y + x/y)/2;
    }
  while((x/y) > y+1);
  return y;
}
