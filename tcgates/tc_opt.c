/* See {tc_opt.h}. */
/* Last edited on 2008-02-12 13:01:21 by stolfi */

#define tc_opt_C_COPYRIGHT "Copyright © 2006 by the State University of Campinas (UNICAMP)"
#define tc_opt_C_AUTHORS "Created 2006-may-20 by Jorge Stolfi, IC-UNICAMP"

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
#include <affirm.h>
#include <argparser.h>

#include <tcgates.h>
#include <tc_opt.h>

#define nprobes_max 50

typedef uint32_t bifun_index_t;
  /* The index of a bifun in a {tc_opt_table_t}. */

struct tc_opt_table_t
  {
    int m;                           /* Number of input wires. */
    int n;                           /* Number of internal state wires. */
    bifun_ct_t max_bifuns;           /* Max number of bifuns to generate. */
    int max_gates;                   /* Max number of gates in any circuit. */
    bifun_hash_t hash_size;          /* Hash table size. */

    bifun_t bifun_null;              /* Null {bifun_t} in table. */
    bifun_index_t bifun_index_null;  /* Null {bifun_index_t} in table. */
    tcgate_t tcgate_null;            /* Null {tcgate_t} in table. */
    bifun_t *fh;                     /* Bifuns, indexed by {bifun_index_t}. */
    tcgate_t *ge;                    /* The optimal circuit for bifun {fh[h]} ends with gate {ge[h]}. */
    bifun_index_t *hr;               /* The optimal circuit for bifun {fh[h]} begins with the circuit of {fh[hr[h]]}. */
    
    bifun_ct_t nq;        /* Number of bifuns in table. */
    bifun_index_t *hq;    /* Indices of all bifuns in table. */
  }; 
  /* A table of non-identity bifuns from {B^m} to {B^n},
    with their optimal circuits.
  
    The parallel arrays {fh[0..hash_size-1]}, {ge[0..hash_size-1]},
    and {hr[0..hash_size-1]} contain data on all non-identity bifuns
    with a certain input size {m} and internal state size {n}. They
    are indexed by the bifun index {h}, as computed by
    {tc_opt_table_locate}.
    
    Any entry with {fh[h] == bifun_null} is vacant; otherwise {fh[h]}
    is a non-identity bifun whose optimal circuit has been computed.
    An optimal TCC for the bifun {fh[h]} then consists of an optimal
    TCC for {fh[hr[h]]} (which is empty if {hr[h] ==
    bifun_index_null}), followed by the TC gate {ge[h]}.

    The list of bifun indices of all distinct bifuns in the table,
    including the identity bifun, without gaps, is {hq[0..nq-1]}.
    This list is sorted in order of increasing complexity.
  */

/* INTERNAL PROTOTYPES */

tc_opt_table_t *tc_opt_table_new(int m, int n, bifun_ct_t max_bifuns, int max_gates);
  /* Allocates a new {tc_opt_table_t} record for at least {max_bifuns}
    bifuns from {B^m} to {B^n}, initially containing no bifuns. */

void tc_opt_table_fill(tc_opt_table_t *tb, int np, tcgate_t gp[]);
  /* Fills the table {tb}, assumed to be empty, with all bifuns that can be obtained by 
    combining the TC gates {gp[0..np-1]}, in all possible ways.  The bifuns are
    generated in order of increasing complexity, until all possible bifuns
    have been generated, or {tb->max_bifuns} have been generated. */

bifun_index_t tc_opt_table_locate(tc_opt_table_t *tb, bifun_t f);
  /* Given a bifun {f}, locates its index {h} in the table {tb}. If there is
    no such entry, returns the {h} of a vacant entry {fh} where {f}
    should be stored. */

void tc_opt_table_print_tcc_from_index_v(FILE *wr, tc_opt_table_t *tb, bifun_index_t h);
  /* Writes to {wr} the optimal TCC for the bifun with index {h} in {tb}, in horizontal format. */
    
void tc_opt_table_print_tcc_from_index_h(FILE *wr, tc_opt_table_t *tb, bifun_index_t h);
  /* Writes to {wr} the optimal TCC for the bifun with index {h} in {tb}, in vertical format. */
    
int tc_opt_table_tcc_size_from_index(tc_opt_table_t *tb, bifun_index_t h);
  /* Returns the size of the optimal TCC for the bifun with index {h} in {tb}. */

int tc_opt_table_tcc_size_max(tc_opt_table_t *tb);
  /* Maximum value that can be returned by {tc_opt_table_tcc_size_from_index}. */

bifun_hash_t tc_opt_table_hash_size(bifun_ct_t max_bifuns);
  /* Adequate size of a hash table for storing {max_bifuns} bifuns. */

/* IMPLEMENTATIONS */

tc_opt_table_t *tc_opt_table_new(int m, int n, bifun_ct_t max_bifuns, int max_gates)
  {
    /* Create the table head record: */
    tc_opt_table_t *tb = notnull(malloc(sizeof(tc_opt_table_t)), "no mem for {tb}");
    
    /* Fill its fields: */
    tb->m = m;
    tb->n = n;
    tb->max_bifuns = max_bifuns;   /* Max number of bifuns to generate. */
    tb->max_gates = max_gates;     /* Max number of gates allowed in any bifun. */
    tb->hash_size = tc_opt_table_hash_size(max_bifuns);   /* Hash table size. */

    tb->bifun_null = bifun_ident(m, n);    /* Null bifun in hash table. */
    tb->bifun_index_null = tb->hash_size;        /* Null hash index of bifun in hash table. */
    tb->tcgate_null = 0;                   /* Null TC gate in hash table. */
    tb->fh = notnull(malloc(tb->hash_size*sizeof(bifun_t)), "no mem for {fh}");       /* Bifuns. */
    tb->ge = notnull(malloc(tb->hash_size*sizeof(tcgate_t)), "no mem for {ge}");      /* End gates. */
    tb->hr = notnull(malloc(tb->hash_size*sizeof(bifun_index_t)), "no mem for {hr}");  /* TCC prefix. */

    tb->hq = notnull(malloc(tb->max_bifuns*sizeof(bifun_index_t)), "no mem for {hq}");  /* Hash indices of built bifuns. */
    
    /* Make all hash table entries {fh[h]} vacant, and the list {hq} empty: */
    bifun_index_t h;
    for (h = 0; h < tb->hash_size; h++)
      { tb->fh[h] = tb->bifun_null;      /* Mark slot as vacant. */
        tb->ge[h] = tb->tcgate_null;     /* For tidiness. */
        tb->hr[h] = tb->bifun_index_null; /* For tidiness. */
      }
    tb->nq = 0;                                                             /* Number of bifuns built so far. */

    return tb;
  }

bifun_hash_t tc_opt_table_hash_size(bifun_ct_t max_bifuns)
  {
    bifun_hash_t hash_size = 2*max_bifuns + 1;
    assert((hash_size - 1)/2 == max_bifuns); /* Overflow check. */
    return hash_size;
  }

bifun_index_t tc_opt_table_locate(tc_opt_table_t *tb, bifun_t f)
  { /* Compute initial index {h} by hashing {f}: */
    bifun_index_t h = bifun_hash(f, tb->hash_size);
    if ((tb->fh[h] == f) || (tb->fh[h] == 0)) { return h; }
    /* Search along {f}'s path: */
    int nprobes = 1;
    bifun_index_t step = bifun_hash_step(f, tb->hash_size); /* Actually computed later. */
    while (nprobes < nprobes_max)
      { h = (h + step) % tb->hash_size;
        if ((tb->fh[h] == f) || (tb->fh[h] == 0)) { return h; }
        nprobes++;
      }
    assert(FALSE); /* if stopped by {nprobes >= nprobes_max}, fail. */
    return 0;
  }

void tc_opt_table_print_tcc_from_index_h(FILE *wr, tc_opt_table_t *tb, bifun_index_t h)
  { if (h != tb->bifun_index_null)
      { tc_opt_table_print_tcc_from_index_h(wr, tb, tb->hr[h]);
        tcgate_print(wr, " ", tb->n, tb->ge[h], NULL, FALSE);
      }
  }

void tc_opt_table_print_tcc_from_index_v(FILE *wr, tc_opt_table_t *tb, bifun_index_t h)
  { if (h != tb->bifun_index_null)
      { tc_opt_table_print_tcc_from_index_v(wr, tb, tb->hr[h]);
        tcgate_print(wr, "  ", tb->n, tb->ge[h], "\n", FALSE);
      }
  }

int tc_opt_table_tcc_size_max(tc_opt_table_t *tb)
  { /* Assumes that {tb->hq[0..tb->nq-1]} is sorted by increasing circuit size: */
    return tc_opt_table_tcc_size_from_index(tb, tb->hq[tb->nq-1]);
  }

int tc_opt_table_tcc_size_from_index(tc_opt_table_t *tb, bifun_index_t h)
  { if (h == tb->bifun_index_null)
      { return 0; }
    else
      { return 1 + tc_opt_table_tcc_size_from_index(tb, tb->hr[h]); }
  }

void tc_opt_table_print_table(FILE *wr, tc_opt_table_t *tb)
  {
    int szmax = tc_opt_table_tcc_size_max(tb); /* Max TC complexity of any bifun. */
    bifun_ct_t nsz[szmax+1]; /* {nsz[sz]} is the number of {(m,n)}-bifuns with complexity {sz}. */
    bifun_ct_t psz[szmax+1]; /* {psz[sz]} is the number of proper {(m,n)}-bifuns with complexity {sz}. */
    int sz;
    for (sz = 0; sz <= szmax; sz++) { nsz[sz] = 0; psz[sz] = 0; }
    bifun_ct_t kq;
    for (kq = 0; kq < tb->nq; kq++)
      { 
        bifun_index_t h = tb->hq[kq];
        bifun_t f = tb->fh[h];
        sz = tc_opt_table_tcc_size_from_index(tb, h);
        assert(sz <= szmax); 
        bool_t proper = bifun_is_proper(tb->m, tb->n, f); 
        if (proper) { psz[sz]++; }
        nsz[sz]++;
        assert(nsz[sz] != 0); /* Overflow check. */
        fprintf(wr, "%s fq[%010d]", (proper ? "*" : "-"), kq);
        bifun_packed_print(wr, "  ", tb->m, tb->n, f, "  ");
        tc_opt_table_print_tcc_from_index_h(wr, tb, h);
        fputc('\n', wr);
      }

    fprintf(wr, "# \n");
    fprintf(wr, "# optimal (%d,%d)-TC circuits, general and proper, by size:\n", tb->m, tb->n);
    fprintf(wr, "# \n");
    for (sz = 0; sz <= szmax; sz++) 
      { fprintf(wr, "# %2d %10d %10d\n", sz, nsz[sz], psz[sz]); }
    fflush(wr);
  }

void tc_opt_table_fill(tc_opt_table_t *tb, int np, tcgate_t gp[])
  {
    /* Start with the identity function and the empty circuit: */
    tb->hq[0] = tb->bifun_index_null;
    tb->nq = 1;
    int nsize = 0; /* Size of last circuit in queue. */
    int csize = -1; /* Size of last circuit that was processed. */
    /* Enumerate all optimal TCCs in order  */
    bifun_ct_t kq = 0;
    while (kq < tb->nq)
      { bool_t debug = ((kq < 50) || (kq >= tb->max_bifuns - 50));
        bifun_index_t hk = tb->hq[kq];
        int ksize = tc_opt_table_tcc_size_from_index(tb, hk); 
        if (ksize > csize)
          { fprintf(stderr, "begin extending circuits of size %d to size %d (kq = %d).\n", ksize, ksize+1, kq);
            csize = ksize; 
            if (ksize >= tb->max_gates)
              { fprintf(stderr, "interruped -- reached circuit size limit.\n");
                return;
              }
          }
        /* Try to compose {tb->fq[kq]} with all TC gates: */
        tcgate_ct_t jp;
        for (jp = 0; jp < np; jp++)
          { 
            bifun_t fk = (hk == tb->bifun_index_null ? bifun_ident(tb->m, tb->n) : tb->fh[hk]);
            tcgate_t gj = gp[jp];
            if (debug) { bifun_print_h(stderr, "  fk = ", tb->m, tb->n, fk, "\n"); }
            if (debug) { tcgate_print(stderr, "  gj =  ", tb->n, gj, "\n", FALSE); }
            bifun_t fc = bifun_tcgate_compose(tb->m, tb->n, fk, gj);
            if (debug) { bifun_print_h(stderr, "  fc = ", tb->m, tb->n, fc, "\n"); }
            if (fc != tb->bifun_null)
              { bifun_index_t hc = tc_opt_table_locate(tb, fc);
                if (tb->fh[hc] == tb->bifun_null)
                  { /* Bifun {fc} is new; save it and its implementation: */
                    if (debug) { fprintf(stderr, "  new"); }
                    /* Check whether we can add it: */
                    if (tb->nq >= tb->max_bifuns) 
                      { fprintf(stderr, "interruped -- too may bifuns.\n");
                        return;
                      }
                    /* Bifun {fc} is new; save it and its implementation: */
                    if (debug) { fprintf(stderr, " - saved in fh[%010u]", hc); }
                    tb->fh[hc] = fc; tb->ge[hc] = gj; tb->hr[hc] = hk;
                    if (debug) { fprintf(stderr, " - enqueued in hq[%010d]", tb->nq); }
                    tb->hq[tb->nq] = hc; (tb->nq)++;
                    assert((nsize == csize) || (nsize == csize + 1));
                    nsize = csize + 1;
                  }
                else
                  { if (debug) { fprintf(stderr, "  bifun is old"); } }
              }
            else
              { if (debug) { fprintf(stderr, "  bifun is identity"); } }
            if (debug) { fprintf(stderr, "\n"); }
          }
        kq++;
      }
    fprintf(stderr, "success -- enumerated all bifuns.\n");
  }

tc_opt_table_t *tc_opt_table_build(int m, int n, bifun_ct_t max_bifuns, int max_gates)
  {
    fprintf(stderr, "listing optimal TC circuits with %d inputs and %d wires.\n", m, n);
    
    demand(m <= m_max, "too many input bits");
    demand(n <= n_max, "too many state bits");
  
    /* Checking whether the data types are large enough: */
    demand(sizeof(word_t)*8 >= n, "the type word_t is too small");
    demand(((word_ct_t) -1) < 0, "the type word_ct_t is unsigned");
    
    word_ct_t mw = word_count(m); /* Includes a check for overflow. */
    word_ct_t nw = word_count(n); /* Includes a check for overflow. */
    assert(nw > 0); /* To pacify the compiler. */

    demand(sizeof(tcgate_t)*8 >= 2*n, "the type tcgate_t is too small");
    demand(((tcgate_ct_t) -1) < 0, "the type tcgate_ct_t is unsigned");

    tcgate_ct_t ng = tcgate_count(n); /* Number of {n}-gates (including no-ops). */
    tcgate_ct_t np = tcgate_proper_count(n); /* Number of {n}-gates (excluding no-ops). */
    
    fprintf(stderr, "there are %d TC gates, of which %d are proper.\n", ng, np);

    demand(sizeof(bifun_t)*8 >= mw*n, "");
    demand(((bifun_ct_t) -1) < 0, "the type bifun_ct_t is unsigned");
    
    /* Make list of non-trivial {n}-gates: */
    tcgate_t gp[np];
    { tcgate_ct_t jp = 0, kg;
      for (kg = 0; kg < ng; kg++)
        { tcgate_t g = kg;
          if (! tcgate_is_no_op(n, g)) { assert(jp < np); gp[jp]= g; jp++; }
        }
      assert(jp == np);
    }

    /* Compute the number of distinct bifuns: */
    bifun_ct_t tot_bifuns = bifun_count(m, n);    /* Number of bifuns from {B^m} to {B^n}. */
    fprintf(stderr, "there are %d distinct (%d,%d)-bifuns.\n", tot_bifuns, m, n);

    /* Choose the table size: */
    if (max_bifuns > tot_bifuns) 
      { max_bifuns = tot_bifuns; }
    else if (max_bifuns < tot_bifuns) 
      { fprintf(stderr, "generating only the first %d of them.\n", max_bifuns); }

    if (max_gates < INT_MAX) 
      { fprintf(stderr, "generating only circuits with at most %d gates.\n", max_gates); }
    
    /* Create the optimal TCC table, initially empty: */
    tc_opt_table_t *tb = tc_opt_table_new(m, n, max_bifuns, max_gates);
    
    /* Fill it: */
    tc_opt_table_fill(tb, np, gp);
    fprintf(stderr, "generated %d functions\n", tb->nq);
    if (max_bifuns >= tot_bifuns) { assert(tb->nq == tot_bifuns); }
    
    return tb;      
  }
