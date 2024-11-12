/* See {logsys.h}. */
/* Last edited on 2013-10-02 03:00:03 by stolfilocal */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <assert.h>

#include <bool.h>
#include <affirm.h>
#include <jsrandom.h>

#include <logsys.h>

#include <logsys_arith.h>
  
/* IMPLEMENTATIONS */

#define ONE64 ((uint64_t)1)
#define ALL64 ((uint64_t)(-1LL))

logsys_t *logsys_build_mul(int nx, int ny, logsys_va_t *x[], logsys_va_t *y[], logsys_va_t *z[])
  {
    /* In the comments below, {num(x[r..s])} means {SUM{ x[i]*2^(i-r) : i IN r..s }}. */
    
    int nz = nx + ny; /* Number of output bits. */
   
    /* Create the system: */
    logsys_t *S = logsys_new_system();
    
    int i, j, k;
    
    /* Create the input variables (bits of {X} and {Y}): */
    for (i = 0; i < nx; i++) { x[i] = logsys_add_variable(S); }
    for (j = 0; j < ny; j++) { y[j] = logsys_add_variable(S); }
    
    /* Initialize the output variables (bits of {Z}) to NULL: */
    for (k = 0; k < nz; k++) { z[k] = NULL; }

    /* Add the circuitry: */
    logsys_va_t *p[ny];  /* Product of one bit of {X} by {Y}. */
    for (i = 0; i < nx; i++)
      { /* Here {z[0..ny+i-1]} are the bits of {num(x[0..i-1])*Y}. */
        /* Create variables {p[0..ny-1]} for the bits of {num(x[i])*Y}: */
        for (j = 0; j < ny; j++)
          { p[j] = logsys_add_variable(S);
            logsys_va_t *va[3] = {p[j], x[i], y[j]};
            (void)logsys_add_equation(S, logsys_op3_AND, 3, va, NULL);
          }
        /* Now create a circuitry that adds {num(p[0..ny-1])*2^i} onto {Z}: */
        logsys_va_t *c = NULL; /* The carry bit. */
        for (j = 0; j < ny; j++)
          { logsys_va_t *u = p[j];
            logsys_va_t *v = z[i+j];
            logsys_va_t *w = c;
            logsys_add_adder_slice(S, u, v, w, &(z[i+j]), &c); }
        /* Store the carry bit (which may be NULL): */
        z[i+ny] = c;
      }
    return S;
  }
  
void logsys_add_adder_slice
  ( logsys_t *S,
    logsys_va_t *u,
    logsys_va_t *v, 
    logsys_va_t *w, 
    logsys_va_t **s, 
    logsys_va_t **c
  )
  {
    /* Make sure that {v == NULL} implies {w == NULL}: */
    if (v == NULL) { v = w; w = NULL; }
    /* Make sure that {u == NULL} implies {v == NULL}: */
    if (u == NULL) { u = v; v = w; w = NULL; }
    if (u == NULL)
      { /* Adder has no inputs, there is nothing to do: */
        (*s) = NULL; (*c) = NULL; 
      }
    else if (v == NULL)
      { /* Adder has only one input {u}: */
        (*s) = u; (*c) = NULL; 
      }
    else
      { /* Adder has at least two inputs, {u} and {v}. */
        /* Form their sum and carry: */
        logsys_va_t *s1 = logsys_add_variable(S); /* Low bit of {u+v}. */
        logsys_va_t *c1 = logsys_add_variable(S); /* Carry of {u+v}. */
        logsys_va_t *ua[3] = {s1, u, v};
        (void)logsys_add_equation(S, logsys_op3_XOR, 3, ua, NULL); 
        logsys_va_t *va[3] = {c1, u, v};
        (void)logsys_add_equation(S, logsys_op3_AND, 3, va, NULL); 
        if (w == NULL)
          { /* Adder has only two inputs, {u} and {v}: */
            (*s) = s1; (*c) = c1;
          }
        else
          { /* Adder has three inputs, {u}, {v}, {w}. */
            /* Form the sum and carry of {u,v,w}: */
            logsys_va_t *s2 = logsys_add_variable(S); /* Low bit of {s1+w}. */
            logsys_va_t *c2 = logsys_add_variable(S); /* Carry of {s1+w} */
            logsys_va_t *c3 = logsys_add_variable(S); /* Carry of {u+v+w} */
            logsys_va_t *ua[3] = {s2, w, s1};
            (void)logsys_add_equation(S, logsys_op3_XOR, 3, ua, NULL); 
            logsys_va_t *va[3] = {c2, w, s1};
            (void)logsys_add_equation(S, logsys_op3_AND, 3, va, NULL); 
            /* We can combine {c1} and {c2} with OR or XOR since they cannot be both 1. */
            logsys_va_t *wa[3] = {c3, c2, c1};
            (void)logsys_add_equation(S, logsys_op3_XOR, 3, wa, NULL); 
            (*s) = s2; (*c) = c3;
          }
      }
  }

#define uint64_dec_fmt "lu"
  /* Decimal format for printing {uint64_t}. */

void logsys_check_mul(logsys_t *S, int nx, int ny, logsys_va_t *x[], logsys_va_t *y[], logsys_va_t *z[])
  {
    int nz = nx + ny; /* Number of bits in product. */
    int nt = 100; /* Number of tests. */
    
    /* Gather all the system's variables: */
    int nv;
    logsys_va_t **v = NULL;
    logsys_get_variables(S, &nv, &v);
    
    /* Logical value ranges: */
    bool_t lo[nv], hi[nv];  
    
    /* Allocate bit vectors: */
    bool_t xb[nx]; /* Values for the bits of {X}. */
    bool_t yb[ny]; /* Values for the bits of {Y}. */
    bool_t zb[nz]; /* Values for the bits of {Z}. */
    
    int t;
    for (t = 0; t < nt; t++)
      {
        fprintf(stderr, "=== Test number %d ===========================================\n", t);
        /* Generate a new test case: */
        fprintf(stderr, "Generating test case with %d bits in X and %d bits in Y\n", nx, ny);
        uint32_t Xmask = (1LU << (uint32_t)nx) - 1LU;
        uint32_t Ymask = (1LU << (uint32_t)ny) - 1LU;
        uint32_t X = Xmask & (uint32_t)uint64_random();  /* A random {nx}-bit integer. */
        uint32_t Y = Ymask & (uint32_t)uint64_random();  /* A random {ny}-bit integer. */
        uint64_t Z = ((uint64_t)X)*((uint64_t)Y);
        fprintf(stderr, "Test values X = %u  Y = %u  Z = %" uint64_dec_fmt "\n", X, Y, Z);
        /* Extract its bits: */
        logsys_get_bits(X, nx, xb);
        logsys_get_bits(Y, ny, yb);
        logsys_get_bits(Z, nz, zb);
        /* Set initial ranges: */
        int i, j, k;
        logsys_va_id_t iv;
        for (iv = 0; iv < nv; iv++) { lo[iv] = FALSE; hi[iv] = TRUE; }
        for (i = 0; i < nx; i++) { iv = logsys_variable_id(x[i]); lo[iv] = hi[iv] = xb[i]; }
        for (j = 0; j < ny; j++) { iv = logsys_variable_id(y[j]); lo[iv] = hi[iv] = yb[j]; }
        /* Apply multiplier circuit and check result: */
        bool_t ok = logsys_find_solution(S, nv, v, NULL, lo, hi, NULL); 
        demand(ok, "** circuit has no valid state"); 
        /* Check result: */
        for (k = 0; k < nz; k++) 
          { iv = logsys_variable_id(z[k]);
            if ((lo[iv] != zb[k]) || (hi[iv] != zb[k]))
              { fprintf(stderr, "** error in output bit z[%02d]", k); 
                logsys_print_variable_id(stderr, " == v", iv, 4, NULL); 
                fprintf(stderr, "  value = [%d..%d] should be %d\n", (int)(lo[iv]), (int)(hi[iv]), (int)(zb[k]));
                demand(FALSE, "aborted");
              }
          }
        fprintf(stderr, "Test succeeded.\n");
      }
    
  }

void logsys_check_solve_mul
  ( logsys_t *S, 
    int ni, 
    logsys_va_t **iv, 
    int nx, 
    int ny, 
    logsys_va_t *x[], 
    logsys_va_t *y[], 
    logsys_va_t *z[]
  )
  {
    affirm(FALSE, "not implemented");
  }
