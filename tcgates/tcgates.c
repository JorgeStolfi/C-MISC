/* See {tcgates.h}. */
/* Last edited on 2008-02-12 12:08:30 by stolfi */

#define tcgates_C_COPYRIGHT "Copyright © 2006 by the State University of Campinas (UNICAMP)"
#define tcgates_C_AUTHORS "Created 2006-may-20 by Jorge Stolfi, IC-UNICAMP"

/* Must define _GNU_SOURCE in order to get {asprintf} */
#define _GNU_SOURCE
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <values.h>
#include <math.h>

#include <bool.h>
#include <jsmath.h>
#include <affirm.h>
#include <argparser.h>

#include <tcgates.h>

#define DEBUG TRUE
  /* Define this as TRUE to print various diagnostics. */

/* INTERNAL PROTOTYPES */

bifun_hash_t gcd32(bifun_hash_t a, bifun_hash_t b);
  /* Greatest common divisor of {a,b}, or 0 if both are zero. */

/* BIT STRINGS */

word_ct_t word_count(int n)
  { word_ct_t nw = (((word_ct_t)1) << n);
    demand(nw > 0, "word_ct_t is too small") ;
    return nw;
  }

void word_print(FILE *wr, char *pre, int m, word_t x, bool_t bin, bool_t dec, char *suf)
  { if ((pre != NULL) && ((*pre) != 0)) { fputs(pre, wr); }
    int xd = x; /* Save {x} for decimal pintout. */
    if (bin)
      { int i;
        for (i = 0; (i < m) || (x != 0); i++)
          { fputc(((x & 1) ? '1' : '0'), wr); x = x >> 1; }
        assert(x == 0);
      }
    if (bin && dec) { fprintf(wr, "="); }
    if (dec) { fprintf(wr, "%d", xd); }
    if ((suf != NULL) && ((*suf) != 0)) { fputs(suf, wr); }
  }

void word_table_print_v(FILE *wr, int ind, int m, int n, word_ct_t ysz, word_t y[])
  { word_ct_t kx;
    for (kx = 0; kx < ysz; kx++) 
      { word_t x = kx;
        fprintf(wr, "%*s", ind, "");
        word_print(wr, "", m, x, TRUE, TRUE, "");
        word_print(wr, " -> ", n, y[x], TRUE, TRUE, "\n");
      }
  }

void word_table_print_h(FILE *wr, char *pre, int m, int n, word_ct_t ysz, word_t y[], char *suf)
  { if ((pre != NULL) && ((*pre) != 0)) { fputs(pre, wr); }
    word_ct_t kx;
    for (kx = 0; kx < ysz; kx++) 
      { word_t x = kx;
        if (x > 0) { fprintf(wr, " "); }
        word_print(wr, "", m, x, FALSE, TRUE, "");
        fprintf(wr, "->");
        word_print(wr, "", n, y[x], TRUE, TRUE, "");
      }
    if ((suf != NULL) && ((*suf) != 0)) { fputs(suf, wr); }
  }

/* TC GATES

  A TC gate is encoded as a binary number where bits {2*k} and {2*k+1}
  describe the element sitting on wire {k}, with the conventions
  00 = pass, 01 = invert, 10 = 0-test, 11 = 1-test.
  
  All TC gates have {n_max} slices, but the encoding is such that
  gates which effectively use only the lowest {m} wires have
  consecutive encodings in {0 .. 4^m-1}. Note that many of those are
  no-ops. */

void get_tcgate_masks(int n, tcgate_t g, word_t *a, word_t *b);
  /* Splits an {n}-wire TC gate {g} into two binary masks {a,b}. The element on
    wire {i} is encoded as the pair {(a[i],b[i])}  */

typedef uint8_t bit_t;

char elem_symbol(bit_t a, bit_t b);
  /* The symbolic code for element {e} with mask bits {a,b}: which is one
    of '×' (invert), '-' (pass), '0' (0-test), and '1' (1-test). */

void get_tcgate_masks(int n, tcgate_t g, word_t *a, word_t *b)
  { tcgate_ct_t gw = tcgate_count(n);
    assert(g < gw);
    /* Masks being built: */
    word_t aa = 0, bb = 0;
    int k = 0;    /* Current wire. */
    word_t z = 1; /* {z = 2^i}. */
    while (k < n)
      { if ((g & 1) != 0) { aa |= z; }
        if ((g & 2) != 0) { bb |= z; }
        g = g >> 2;
        k++; z <<= 1;
      }
    /* Return masks: */
    (*a) = aa; (*b) = bb;
  }

tcgate_ct_t tcgate_count(int n)
  { tcgate_ct_t ng = ((tcgate_ct_t)1) << 2*n;
    demand(ng > 0, "the type tcgate_ct_t is too small");
    return ng;
  }

tcgate_ct_t tcgate_proper_count(int n)
  { /* Assumes that a {tcgate_ct_t} is large enough for {4^n}: */
    tcgate_ct_t ng = 1, og = 1; while(n > 0) { ng *= 4; og *= 3; n--; }
    return ng - og;
  }

word_t tcgate_apply(int n, tcgate_t g, word_t x)
  { bool_t debug = FALSE;
    if (debug) { word_print(stderr, "\nx = ", n, x, TRUE, TRUE, "\n"); }
    word_t a, b;
    get_tcgate_masks(n, g, &a, &b);
    if (debug) { word_print(stderr, "a = ", n, a, TRUE, FALSE, "\n"); }
    if (debug) { word_print(stderr, "b = ", n, b, TRUE, FALSE, "\n"); }
    word_t t = a & (x ^ b);
    if (debug) { word_print(stderr, "t = ", n, t, TRUE, FALSE, "\n"); }
    word_t y = (~a) & b;
    if (debug) { word_print(stderr, "y = ", n, y, TRUE, FALSE, "\n"); }
    if (t == 0) { x = x ^ y; }
    if (debug) { word_print(stderr, "y = ", n, x, TRUE, TRUE, "\n"); }
    return x;
  }

bool_t tcgate_is_no_op(int n, tcgate_t g)
  { word_t a, b;
    get_tcgate_masks(n, g, &a, &b);
    return (((~a) & b) == 0);
  }

char elem_symbol(bit_t a, bit_t b)
  { if (a == 0)
      { return (b == 0 ? '-' : '×'); } 
    else
      { return (b == 0 ? '0' : '1'); } 
  }

void tcgate_print(FILE *wr, char *pre, int n, tcgate_t g, char *suf, bool_t masks)
  { if ((pre != NULL) && ((*pre) != 0)) { fputs(pre, wr); }
    word_t a, b;
    get_tcgate_masks(n, g, &a, &b);
    if (masks)
      { word_print(wr, "", n, a, TRUE, FALSE, "");
        word_print(wr, "|", n, b, TRUE, FALSE, " = ");
      }
    int k;      /* Current wire. */
    for (k = 0; k < n; k++)
      { char c = elem_symbol(a & 1, b & 1);
        fprintf(wr, "%c", c);
        a >>= 1; b >>= 1;
      }
    if ((suf != NULL) && ((*suf) != 0)) { fputs(suf, wr); }
  }

/* BIFUNS
  
  An {bifun_t} that represents a map from {B^m} to {B^n} consists of
  {2^m} fields of {n} bits each, packed into an {int64_t}, starting
  from the low end.
  
  Field number {x} contains {(f(x) - x) % 2^n}. This repersentation
  has the peoeprty that, for any {m,n}, the identity bifun has all
  fields equal to zero. */

bifun_t bifun_set(int m, int n, bifun_t f, word_t x, word_t fx);
  /* Returns the {m,n}-bifun {f} modified so that it maps word {x}
    (which should be in {B^m}) to word {fx} (which should be in {B^n}).
    Beware that the result may not be a valid bifun.  This procedure is
    used only as the basic step of {pack} below. */

bifun_ct_t bifun_count(int m, int n)
  { if (m > n) { return 0; }
    word_ct_t mw = word_count(m); /* Number of {m}-words, namely {2^m}. */
    word_ct_t nw = word_count(n); /* Number of {n}-words, namely {2^n}. */
    bifun_ct_t mnf = 1;
    word_ct_t kx;
    for (kx = 0; kx < mw; kx++) 
      { word_ct_t nch = nw - kx; /* Number of choices for {f(kx)}. */
        bifun_ct_t tt = mnf*nch;
        assert(tt/mnf == nch); /* Check for overflow. */
        mnf = tt;
      }
    return mnf;
  }

word_t bifun_apply(int m, int n, bifun_t f, word_t x)
  { assert(m <= m_max);
    assert(n <= n_max);
    word_ct_t mw = word_count(m);
    assert(x <= mw);
    word_ct_t nw = word_count(n);
    word_t nmask = nw - 1;
    word_t fx = (((f >> x*n) & nmask) + x) & nmask;
    return fx;
  }

bifun_t bifun_set(int m, int n, bifun_t f, word_t x, word_t fx)
  { assert(m <= m_max);
    assert(n <= n_max);
    word_ct_t mw = word_count(m);
    assert(x <= mw);
    word_ct_t nw = word_count(n);
    word_t nmask = nw - 1;
    word_t oslot = (f >> x*n) & nmask; /* Current contents of the {x}-slot of {f}. */
    word_t nslot = ((fx + nw - x ) & nmask); /* new contents of that slot. */
    bifun_t g = f ^ ((oslot ^ nslot) << x*n); 
    return g;
  }

void unpack_bifun(int m, int n, bifun_t f, word_t y[]);
  /* Assumes that {f} is an bifun from {B^m} to {B^n}. Stores in
    {t[x]} the value of {f(x)}, for all words {x} in {0..2^m-1}. */

bifun_t pack_bifun(int m, int n, word_t y[]);
  /* Assumes that {y} is an array of {2^m} words, each with {n} bits.
    Returns the bifun {f} from {B^m} to {B^n} such that {f(x) = y[x]}
    for all {x} in {0..2^m-1}. */

void unpack_bifun(int m, int n, bifun_t f, word_t y[])
  { assert(m <= m_max);
    assert(n <= n_max);
    word_ct_t mw = word_count(m);
    word_ct_t kx;
    for (kx = 0; kx < mw; kx++)
      { word_t x = kx;
        y[x] = bifun_apply(m, n, f, x);
      }
  }

bifun_t pack_bifun(int m, int n, word_t y[])
  { assert(m <= m_max);
    assert(n <= n_max);
    word_ct_t mw = word_count(m);
    word_ct_t kx;
    bifun_t f = 0;
    for (kx = 0; kx < mw; kx++)
      { word_t x = kx;
        f = bifun_set(m, n, f, x, y[x]);
      }
    return f;
  }

bifun_t bifun_ident(int m, int n)
  { return 0; }

bool_t bifun_is_proper(int m, int n, bifun_t f)
  { assert(m <= m_max);
    assert(n <= n_max);
    word_ct_t mw = word_count(m);
    word_ct_t nw = word_count(n);
    /* Get a mask {himask} for the higher {n-m} bits of a word: */
    word_t himask = nw - mw;
    /* Check whether {f(x)} is in {B^m} for all {x} in {B^m}: */ 
    word_ct_t kx;
    for (kx = 0; kx < mw; kx++) 
      { word_t x = kx;
        word_t fx = bifun_apply(m, n, f, x);
        if ((fx & himask) != 0) { return FALSE; } 
      }
    return TRUE;    
  }

bifun_t bifun_tcgate_compose(int m, int n, bifun_t f, tcgate_t g)
  { assert(m <= m_max);
    assert(n <= n_max);
    word_ct_t mw = word_count(m);
    word_ct_t nw = word_count(n);
    word_t nmask = nw - 1;
    word_t yf[mw]; 
    unpack_bifun(m, n, f, yf);
    word_t yc[mw];
    word_ct_t kx;
    for (kx = 0; kx < mw; kx++) 
      { word_t x = kx;
        yc[x] = tcgate_apply(n, g, yf[x] & nmask) & nmask;
      }
    return pack_bifun(m, n, yc);
  }

bifun_hash_t bifun_hash(bifun_t f, bifun_hash_t M)
  { assert(sizeof(bifun_t)*8 <= 64);
    return ((uint64_t)f) % ((uint64_t)M);
  }

bifun_hash_t bifun_hash_step(bifun_t f, bifun_hash_t M)
  { /* The formula is based on the Linux {drand48} generator: */
    uint64_t a = 0x5DEECE66Dll;
    uint64_t c = 0xBll;
    uint64_t ff = (f + (f >> 31) + (f >> 63)) & 0xFFFFFFFFll; 
    bifun_hash_t h = (a*(uint64_t)ff + c) % (uint64_t)M;
    while (gcd32(M,h) != 1) { h = (h + 1) % M; }
    return h;
  }

void bifun_packed_print(FILE *wr, char *pre, int m, int n, bifun_t f, char *suf)
  { 
    if ((pre != NULL) && ((*pre) != 0)) { fputs(pre, wr); }
    assert(m <= m_max);
    assert(n <= n_max);
    word_ct_t mw = word_count(m);
    word_ct_t nw = word_count(n);
    word_t yf[mw]; unpack_bifun(m, n, f, yf);
    int d = digits(nw-1);
    char sep = '#';
    word_ct_t kx;
    for (kx = 0; kx < mw; kx++) 
      { word_t x = kx;
        fputc(sep, wr); sep = '.';
        fprintf(wr, "%0*d", d, yf[x]);
      }
    if ((suf != NULL) && ((*suf) != 0)) { fputs(suf, wr); }
  }

void bifun_print_h(FILE *wr, char *pre, int m, int n, bifun_t f, char *suf)
  { 
    assert(m <= m_max);
    assert(n <= n_max);
    word_ct_t mw = word_count(m);
    bifun_packed_print(wr, pre, m, n, f, NULL);
    word_t yf[mw]; unpack_bifun(m, n, f, yf);
    word_table_print_h(wr, " = ", m, n, mw, yf, suf);
  }

void bifun_print_v(FILE *wr, char *head, int ind, int m, int n, bifun_t f, char *foot)
  { assert(m <= m_max);
    assert(n <= n_max);
    word_ct_t mw = word_count(m);
    if (head != NULL) { fprintf(wr, "%s\n", head); }
    bifun_packed_print(wr, NULL, m, n, f, NULL);
    word_t yf[mw]; 
    unpack_bifun(m, n, f, yf);
    word_table_print_v(wr, ind + 2, m, n, mw, yf);
    if (foot != NULL) { fprintf(wr, "%s\n", foot); }
  }

bifun_hash_t gcd32(bifun_hash_t a, bifun_hash_t b)
  { while (b > 0) { bifun_hash_t r = a % b; a = b; b = r; }
    return a;
  }
