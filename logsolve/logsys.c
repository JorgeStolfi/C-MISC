/* See {logsys.h}. */
/* Last edited on 2013-10-02 02:59:11 by stolfilocal */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <values.h>

#include <bool.h>
#include <affirm.h>
#include <jsmath.h>
#include <vec.h>

#include <logsys.h>
#include <logsys_def.h>
#include <logsys_sort_H0_H1.h>
#include <logsys_sort_H2.h>
#include <logsys_sort_H3.h>
#include <logsys_sort_H4.h>
  
/* INTERNAL PROTOTYPES */

logsys_op_t logsys_op_mask_for_arg(int n, int j, bool_t val);
  /* Returns the truth table for {n} variables which is
    1 for all cases where variable {j} has value {val}.
    If {j >= n} or {j < 0} ignores {j,val}. */

void logsys_set_arg(logsys_eq_t *eq, int j, logsys_va_t *va);
  /* Sets argument {j} of {eq} as variable {va}, fixing pointers as
    needed. That argument of {eq} must be currently null. (This may
    happen only while the equation is being built, modified or
    dismantled.) */

logsys_va_t *logsys_unset_arg(logsys_eq_t *eq, int j);
  /* Sets argument {j} of {eq} to null, fixing pointers as needed. (This
    operation is valid only while the equation is being built, modified
    or dismantled.) */

void logsys_print_op0(FILE *wr, char *pre, logsys_op0_t op, char *suf);
void logsys_print_op1(FILE *wr, char *pre, logsys_op1_t op, char *suf);
void logsys_print_op2(FILE *wr, char *pre, logsys_op2_t op, char *suf);
void logsys_print_op3(FILE *wr, char *pre, logsys_op3_t op, char *suf);
void logsys_print_op4(FILE *wr, char *pre, logsys_op4_t op, char *suf);
void logsys_print_op5(FILE *wr, char *pre, logsys_op5_t op, char *suf);
void logsys_print_op6(FILE *wr, char *pre, logsys_op6_t op, char *suf);
  /* Specialized opcode printing operations. */

void logsys_intersect_eq_args(logsys_eq_t *a, logsys_eq_t *b, int *snP, logsys_va_t *vc[], int ka[], int kb[]);
  /* Gathers in {vc[0..sn-1]} the intersection of the argument variables of equations {a} and {b}.
    The count {sn} is returned in {*snP}. Also stores in {ka[i]} and {kb[i]},
    for each {i} in {0..sn-1}, the position of {vc[i]} in those argument lists.
    The arrays {vc,ka,kb} had better be allocated with {a->n + b->n} or more 
    elements by the caller. */

void logsys_unite_eq_args(logsys_eq_t *a, logsys_eq_t *b, int *unP, logsys_va_t *vc[], int ka[], int kb[]);
  /* Gathers in {vc[0..un-1]} the union of the argument variables of equations {a} and {b}.
    The count {un} is returned in {*unP}. Also stores in {ka[i]} and {kb[i]},
    for each {i} in {0..un-1}, the position of {vc[i]} in those argument lists.
    The index {ka[i]} is set to -1 if the variable {vc[i]} is not used in {a}.
    Ditto for {kb} and {b}. The arrays {vc,ka,kb} had better be allocated
    with {a->n + b->n} or more elements by the caller. */

bool_t logsys_try_to_combine_with_others(logsys_eq_t *eq, logsys_eq_t **fqP, int *nstP, logsys_eq_vec_t *steq);
  /* Tries to combine {**fqP} with every equation {gq} that has one or more
    variables in common with {fq} and smaler or equal arity. If
    any equations are modified by this process, tries to split nd
    simplify them. 
    
    Returns TRUE as soon as it succeeds in changing some equation.
    May delete the equation {**fqP} itself; in that case sets {*fqP} to NULL.
    If {steq} isnot NULL, any equations that are modified but not deleted are stored into
    {steq} starting at position {steq.e[(*nstP)]}, and {*nstP} is incremented
    accordingly. */
  
logsys_eq_vec_t logsys_collect_equations(logsys_t *S, bool_t tight);
  /* Returns a vector {veq} with all the equations in the system {S}.
    
    If {tight} is true, the elements {veq.e[0..veq.ne-1]} are
    all the equations, in the order they appear in the circular 
    list {S->eq,.nxeq},without repetitions or NULLs.  
    
    If {tight} is false, each equation {eq} of {S} is stored 
    into {veq.e[eq->id]}; unused elements of {veq} are set to NULL. */
    
void logsys_sort_equations(logsys_t *S);
  /* Rearranges the list of equations of {S} (defined by {S->eq,eq->nxeq,eq->preq}) 
    so that they are in increasing order of arity, with {S->eq} pointing
    to the firstequation with minimum arity. */

void logsys_stack_equation(int *nstP, logsys_eq_vec_t *steq, logsys_eq_t *eq);
  /* Assumes that {steq.e[0..(*nstP)-1]} is a stack of equations. If
    [eq} is not in the srack already, pushes {eq} onto that stack,
    expanding {steq} if needed. A no-op if both {nstP} and {steq} are
    NULL. */

void logsys_remove_equation_from_stack(int *nstP, logsys_eq_vec_t *steq, logsys_eq_t *eq);
  /* Removes all occurrences of {eq} from the stack {steq.e[0..(*nstP)-1]}, 
     rearranging the other entries to fill the holes and updating {*nstP} accordingly. */

/* IMPLEMENTATIONS */

int logsys_ones(uint64_t X)
  {
    int n = 0;
    while (X != ZERO64) { n += (X & 1); X = X >> 1; }
    return n;
  }

bool_t logsys_op_is_valid(logsys_op_t op, int n)
  {
    demand(n <= logsys_n_MAX, "invalid n");
    int NX = 1 << n; /* Number of assignments to {n} variables. */
    uint64_t nmsk = ALL64 + (n == 6 ? ZERO64 : (ONE64 << (uint64_t)NX));
    return ((op & nmsk) == op);
  }

logsys_op_t logsys_op_mask_for_arg(int n, int j, bool_t val)
  {
    demand(n <= logsys_n_MAX, "invalid n");
    if (n == 0) { return ONE64; }
    int k = n-1;
    bool_t alok = (k == j ? val : FALSE );
    bool_t ahik = (k == j ? val : TRUE);
    uint64_t KX = ONE64 << (uint64_t)k; /* Shift in mask for toggling arg {k}. */
    uint64_t msk0 = logsys_op_mask_for_arg(n-1, j, val);
    return (alok == FALSE ? msk0 : ZERO64) | (ahik == TRUE ? msk0 << KX : ZERO64);
  }

bool_t logsys_op_depends(logsys_op_t op, int n, int j)
  {
    demand(n <= logsys_n_MAX, "invalid n");
    demand((j >= 0) && (j < n), "invalid arg index");
    uint64_t JX = ONE64 << (uint64_t)j;
    uint64_t msk = logsys_op_mask_for_arg(n, j, FALSE);
    return ((op ^ (op >> JX)) & msk) != ZERO64;
  }

bool_t logsys_op_is_functional(logsys_op_t op, int n, int j)
  {
    demand(n <= logsys_n_MAX, "invalid n");
    demand((j >= 0) && (j < n), "invalid arg index");
    uint64_t msk = logsys_op_mask_for_arg(n, j, FALSE);
    uint64_t JX = ONE64 << (uint64_t)j;
    return ((op ^ (op >> JX)) & msk) == msk;
  }

bool_t logsys_op_is_solvable(logsys_op_t op, int n, int j)
  {
    demand(n <= logsys_n_MAX, "invalid n");
    demand((j >= 0) && (j < n), "invalid arg index");
    uint64_t msk = logsys_op_mask_for_arg(n, j, FALSE);
    uint64_t JX = ONE64 << (uint64_t)j;
    return ((op | (op >> JX)) & msk) == msk;
  }

double logsys_entropy(double m0, double m1)
  { 
    demand((m0 >= 0) && (m1 >= 0), "invalid counts");
    double m = m0 + m1;
    if (m == 0)
      { return 1.0; }
    else
      { return log2(m) 
          - (m0 == 0 ? 0.0 : m0*log2(m0)/m) 
          - (m1 == 0 ? 0.0 : m1*log2(m1)/m);
      }
  }

double logsys_op_entropy_reduction(logsys_op_t op, int n, int j)
  {
    demand(op != ZERO64, "entropy reduction not meaningful for FAL operator");
    demand((n >= 0) && (n <= logsys_n_MAX), "invalid n");
    demand((j >= 0) && (j < n), "invalid j");
    double red = 0;  /* Total entropy reduction in variables other than {j}. */
    uint32_t N = (1 << n);
    uint32_t J = (1 << j);
    int k;
    for (k = 0; k < n; k++)
      { /* Arg {k} had some entropy to lose: */
        uint32_t K = (1 << k);
        /* Count in {mk{a}j{b}} the valid assignments of {op} with {x[k]=a} and {x[j]=b} */ 
        int mk0j0 = 0, mk0j1 = 0, mk1j0 = 0, mk1j1 = 0;
        uint32_t X;
        for (X = 0; X < N; X++)
          { if ((op & (ONE64 << X)) != 0)
              { /* {X} is a valid arg assignment for {op}. */
                if ((X & K) == 0) 
                  { if ((X & J) == 0) { mk0j0++; } else { mk0j1++; } }
                else
                  { if ((X & J) == 0) { mk1j0++; } else { mk1j1++; } }
              }
          }
        int mk0 = mk0j0 + mk0j1; /* Valid assignments for {op} with{x[k]=0}. */
        int mk1 = mk1j0 + mk1j1; /* Valid assignments for {op} with{x[k]=1}. */
        int mj0 = mk0j0 + mk1j0; /* Valid assignments for {op} with{x[j]=0}. */
        int mj1 = mk0j1 + mk1j1; /* Valid assignments for {op} with{x[j]=1}. */
        int m = mk0 + mk1; /* Valid assignments for {op}. */
        assert(m > 0);
        /* Entropy of {x[k]} before setting {x[j]}: */
        double epre = logsys_entropy(mk0, mk1);
        /* Entropy of {x[k]} after setting {x[j]}: */
        double epos =
          (mj0 == 0 ? 0.0 : mj0*logsys_entropy(mk0j0, mk1j0)/m) + 
          (mj1 == 0 ? 0.0 : mj1*logsys_entropy(mk0j0, mk1j0)/m);
        red += (epre - epos);
      }
    return red;
  }

logsys_t *logsys_new_system(void)
  {
    logsys_t *S = (logsys_t *)notnull(malloc(sizeof(logsys_t)), "out of mem");
    S->eq = NULL;
    S->va = NULL;
    S->next_va_id = 0;
    S->next_eq_id = 0;
    return S;
  }

logsys_va_t *logsys_add_variable(logsys_t *S)
  { 
    /* Allocate node: */
    logsys_va_t *va = (logsys_va_t *)notnull(malloc(sizeof(logsys_va_t)), "out of mem");
    /* Assign the next available variable number: */
    va->id = S->next_va_id; S->next_va_id++;
    /* No equation currently uses this variable: */
    va->use = NULL;
    /* Insert variable into the list of variables in system: */
    va->sys = S;
    if (S->va == NULL)
      { /* System has no variables: */
        va->nxva = va;
        va->prva = va;
      }
    else
      { /* Insert {va} into the doubly-linked list: */
        logsys_va_t *prva = S->va;
        logsys_va_t *nxva = prva->nxva;
        va->nxva = nxva;
        va->prva = prva;
        prva->nxva = va;
        nxva->prva = va;
      }
    S->va = va;
    return va;
  }

logsys_eq_t *logsys_add_equation(logsys_t *S, logsys_op_t op, int n, logsys_va_t *va[], logsys_eq_t *preq)
  { 
    demand((n > 0) && (n <= logsys_n_MAX), "invalid variable count");
    demand(logsys_op_is_valid(op, n), "invalid opcode for this arity");
    /* Allocate node: */
    logsys_eq_t *eq = (logsys_eq_t *)notnull(malloc(sizeof(logsys_eq_t)), "out of mem");
    /* Assign the next available equation number: */
    eq->id = S->next_eq_id; S->next_eq_id++;
    /* Set the opcode: */
    eq->op = op;
    eq->n = n;
    /* Set all operands {eq->va[0..logsys_n_MAX]} to NULL initially: */
    int j;
    for (j = 0; j < n; j++) 
      { eq->arg[j] = (logsys_arg_t){ .va = NULL, .nxeq = NULL, .preq = NULL }; }
    /* Append equation to list of equations in system: */
    eq->sys = S;
    if (S->eq == NULL)
      { /* System has no equations: */
        demand(preq == NULL, "invalid {preq}");
        eq->nxeq = eq;
        eq->preq = eq;
      }
    else
      { /* Insert {eq} into the doubly-linked list: */
        if (preq != NULL)
          { demand(preq->sys == S, "{preq} from wrong system"); }
        else
          { preq = S->eq; }
        logsys_eq_t *nxeq = preq->nxeq;
        eq->nxeq = nxeq;
        eq->preq = preq;
        preq->nxeq = eq;
        nxeq->preq = eq;
      }
    S->eq = eq;

    /* Set the arguments and insert {eq} into the corresp usage list: */
    for (j = 0; j < n; j++) { logsys_set_arg(eq, j, va[j]); }
    return eq;
  }

void logsys_set_arg(logsys_eq_t *eq, int j, logsys_va_t *va)
  { 
    logsys_t *S = eq->sys;
    int n = eq->n;
    demand((j >= 0) && (j < n), "invalid operand index");
    demand(va != NULL, "variable must be non-null");
    demand(va->sys == S, "variable and equation from different systems");
    demand(eq->arg[j].va == NULL, "argument is not NULL");
    assert(eq->arg[j].preq == NULL);
    assert(eq->arg[j].nxeq == NULL);
    
    /* Arguments must be distinct variables in order for usage links to work: */
    int i;
    for (i = 0; i < n; i++) { demand(eq->arg[j].va != va, "arguments must be distinct vars"); }
        
    /* Set argument {eq->va[j]} to {va}: */
    eq->arg[j].va = va;
    /* Append equation {eq} to list of equations that use {va}: */
    if (va->use == NULL)
      { /* This is the first equation to use {va}. */
        /* Make {eq} into a singleton list: */
        eq->arg[j].nxeq = eq;
        eq->arg[j].preq = eq;
      }
    else
      { /* Find two consecutive eqs {preq,nxeq} in list of users of {va} and their arg indices: */
        logsys_eq_t *preq = va->use;
        int prj = logsys_find_arg(preq, va);
        assert((prj >= 0) && (prj < preq->n) && (preq->arg[prj].va == va));
        
        logsys_eq_t *nxeq = preq->arg[prj].nxeq;
        int nxj = logsys_find_arg(nxeq, va);
        assert((nxj >= 0) && (nxj < nxeq->n) && (nxeq->arg[nxj].va == va));
        
        /* Insert {eq} between {preq} and {nxeq}: */
        eq->arg[j].nxeq = nxeq;
        eq->arg[j].preq = preq;
        preq->arg[prj].nxeq = eq;
        nxeq->arg[nxj].preq = eq;
      } 
    va->use = eq;
  }

logsys_va_t *logsys_unset_arg(logsys_eq_t *eq, int j)
  { 
    logsys_t *S = eq->sys;
    int n = eq->n;
    demand((j >= 0) && (j < n), "invalid operand index");
    
    logsys_va_t *ua = eq->arg[j].va; /* Variable that was operand {j} of {eq}. */
    assert(ua != NULL);
    assert(ua->sys == S);
    assert(ua->use != NULL);
    
    /* Get adjacent equations on list of equations that use {ua}: */
    logsys_eq_t *preq = eq->arg[j].preq;
    logsys_eq_t *nxeq = eq->arg[j].nxeq;
    assert(nxeq != NULL);
    assert(preq != NULL);
    
    /* Remove equation {eq} from list of equations that use {ua}: */
    if (nxeq == eq)
      { /* This was the only equation using variable {ua}: */
        assert(preq == eq);
        assert(nxeq == eq);
        assert(ua->use == eq);
        ua->use = NULL; 
      }
    else
      { /* This was NOT the only equation using variable {ua}. */
        /* Find arg indices of {ua} in {preq,nxeq}: */
        int prj = logsys_find_arg(preq, ua);
        assert((prj >= 0) && (prj < preq->n) && (preq->arg[prj].va == ua));
        assert(preq->arg[prj].nxeq == eq);
        
        int nxj = logsys_find_arg(nxeq, ua);
        assert((nxj >= 0) && (nxj < nxeq->n) && (nxeq->arg[nxj].va == ua));
        assert(nxeq->arg[nxj].preq == eq);
        
        /* Fix pointers in list of {j}-users of {ua} so as to bypass {eq}: */
        nxeq->arg[nxj].preq = preq;
        preq->arg[prj].nxeq = nxeq;
        /* Set the user pointer of {ua} to the equation immediately before {eq}: */
        ua->use = preq;
      } 
    /* Set argument {j} of {eq} to NULL: */
    eq->arg[j].va = NULL;
    eq->arg[j].preq = NULL;
    eq->arg[j].nxeq = NULL;
    
    return ua;
  }

logsys_va_t *logsys_replace_arg(logsys_eq_t *eq, int j, logsys_va_t *va)
  {
    logsys_va_t *ua = logsys_unset_arg(eq, j);
    logsys_set_arg(eq, j, va);
    return ua;
  }

void logsys_delete_equation(logsys_t *S, logsys_eq_t *eq)    
  { 
    demand(eq->sys == S, "removing equation from wrong system");
    
    /* Set all operands to NULL to unlink {eq} from their usage lists: */
    int j;
    for (j = 0; j < eq->n; j++) { (void)logsys_unset_arg(eq, j); }
     
    /* Remove equation from list of equations in system: */
    assert(eq->nxeq != NULL);
    assert(eq->preq != NULL);
    logsys_eq_t *preq = eq->preq;
    logsys_eq_t *nxeq = eq->nxeq;
    
    if (preq == eq)
      { /* This was the only equation in the system: */
        assert(nxeq == eq);
        assert(S->eq == eq);
        S->eq = NULL;
      }
    else
      { /* This was NOT the only equation in the system. */
        /* Fix pointers in list equations so as to bypass {eq}: */
        preq->nxeq = nxeq;
        nxeq->preq = preq;
        /* Set the last equation of the system to be the equation immediately before {eq}: */
        S->eq = preq;
      }
    eq->preq = NULL;
    eq->nxeq = NULL;
    eq->sys = NULL;
    /* Reclaim the space: */
    free(eq);
  }

int logsys_num_args(logsys_eq_t *eq)
  { 
    return eq->n;
  }

int logsys_find_arg(logsys_eq_t *eq, logsys_va_t *va)
  {
    int k = eq->n - 1;
    while((k > 0) && (eq->arg[k].va != va)) { k--; }
    return k;
  }
  
void logsys_print_system(FILE *wr, char *head, logsys_t *S, char *foot)
  {
    if ((head != NULL) && ((*head) != 0)) { fprintf(wr, "%s\n", head); }
    if (S->eq != NULL)
      { /* Print equations: */
        logsys_eq_t *eq = S->eq->nxeq;  /* First equation, nominally. */
        do {
          logsys_print_equation(wr, "  ", eq, "\n");
          eq = eq->nxeq;
        } while (eq != S->eq->nxeq);
      }
    if ((foot != NULL) && ((*foot) != 0)) { fprintf(wr, "%s\n", foot); }
  }

#define uint64_hex_fmt "lx"
  /* Hex format for printing {uint64_t}. */

#define uint64_dec_fmt "lu"
  /* Decimal format for printing {uint64_t}. */

void logsys_print_equation_id(FILE *wr, char *pre, logsys_eq_id_t id, int d, char *suf)
  {
    if ((pre != NULL) && ((*pre) != 0)) { fputs(pre, wr); }
    fprintf(wr, "%0*" uint64_dec_fmt "", d, id);
    if ((suf != NULL) && ((*suf) != 0)) { fputs(suf, wr); }
  }

void logsys_print_variable_id(FILE *wr, char *pre, logsys_va_id_t id, int d, char *suf)
  {
    if ((pre != NULL) && ((*pre) != 0)) { fputs(pre, wr); }
    fprintf(wr, "%0*" uint64_dec_fmt "", d, id);
    if ((suf != NULL) && ((*suf) != 0)) { fputs(suf, wr); }
  }

void logsys_print_equation(FILE *wr, char *pre, logsys_eq_t *eq, char *suf)
  {
    if ((pre != NULL) && ((*pre) != 0)) { fputs(pre, wr); }
    logsys_print_equation_id(wr, "e", eq->id, 4, ": ");
    logsys_print_op(wr, NULL, eq->n, eq->op, FALSE, NULL);
    fputc('(', wr); 
    int j;
    char *sep = "";
    for (j = eq->n-1; j >= 0; j--)
      { logsys_va_t *va = eq->arg[j].va; 
        fputs(sep, wr);
        if (va == NULL)
          { fputs("-----", wr); }
        else
          { fputc('x', wr);
            fprintf(wr, "%d", j);
            fputc(':', wr);
            logsys_print_variable_id(wr, "v", va->id, 4, NULL);
          }
        sep = ",";
      }
    fputc(')', wr);
    if ((suf != NULL) && ((*suf) != 0)) { fputs(suf, wr); }
  }

void logsys_print_op0(FILE *wr, char *pre, logsys_op0_t op, char *suf)
  {
    if ((pre != NULL) && ((*pre) != 0)) { fputs(pre, wr); }
    uint32_t uop = op;
    fprintf(wr, "o%01X", uop);
    if ((suf != NULL) && ((*suf) != 0)) { fputs(suf, wr); }
  }

void logsys_print_op1(FILE *wr, char *pre, logsys_op1_t op, char *suf)
  {
    if ((pre != NULL) && ((*pre) != 0)) { fputs(pre, wr); }
    uint32_t uop = op;
    switch (op)
      {
        case logsys_op1_ZR0: fputs("ZR0", wr); break; 
        case logsys_op1_UN0: fputs("UN0", wr); break; 
        default: fprintf(wr, "o%01X", uop);
      }
    if ((suf != NULL) && ((*suf) != 0)) { fputs(suf, wr); }
  }

void logsys_print_op2(FILE *wr, char *pre, logsys_op2_t op, char *suf)
  {
    if ((pre != NULL) && ((*pre) != 0)) { fputs(pre, wr); }
    uint32_t uop = op;
    switch (op)
      {
        case logsys_op2_ZR0: fputs("ZR0", wr); break; 
        case logsys_op2_UN0: fputs("UN0", wr); break;
        
        case logsys_op2_ZR1: fputs("ZR1", wr); break; 
        case logsys_op2_UN1: fputs("UN1", wr); break; 
        
        case logsys_op2_W01: fputs("W01", wr); break; 
        case logsys_op2_N01: fputs("N01", wr); break; 
        
        case logsys_op2_G01: fputs("G01", wr); break; 
        case logsys_op2_L01: fputs("L01", wr); break; 
        default: fprintf(wr, "o%01X", uop);
      }
    if ((suf != NULL) && ((*suf) != 0)) { fputs(suf, wr); }
  }

void logsys_print_op3(FILE *wr, char *pre, logsys_op3_t op, char *suf)
  {
    if ((pre != NULL) && ((*pre) != 0)) { fputs(pre, wr); }
    uint32_t uop = op;
    switch (op)
      {
        case logsys_op3_ZR0: fputs("ZR0", wr); break; 
        case logsys_op3_UN0: fputs("UN0", wr); break; 
        
        case logsys_op3_ZR1: fputs("ZR1", wr); break; 
        case logsys_op3_UN1: fputs("UN1", wr); break; 
        
        case logsys_op3_ZR2: fputs("ZR2", wr); break; 
        case logsys_op3_UN2: fputs("UN2", wr); break; 
        
        case logsys_op3_AND: fputs("AND", wr); break; 
        case logsys_op3_ORR: fputs("ORR", wr); break; 
        case logsys_op3_XOR: fputs("XOR", wr); break; 
        case logsys_op3_EQV: fputs("EQV", wr); break; 
        
        case logsys_op3_W01: fputs("W01", wr); break; 
        case logsys_op3_N01: fputs("N01", wr); break; 
        
        case logsys_op3_W02: fputs("W02", wr); break; 
        case logsys_op3_N02: fputs("N02", wr); break; 
        
        case logsys_op3_W12: fputs("W12", wr); break; 
        case logsys_op3_N12: fputs("N12", wr); break; 
        
        case logsys_op3_NIL: fputs("NIL", wr); break; 
        case logsys_op3_ANY: fputs("ANY", wr); break; 
        case logsys_op3_FAN: fputs("FAN", wr); break; 
        default: fprintf(wr, "o%02X", uop);
      }
      
    if ((suf != NULL) && ((*suf) != 0)) { fputs(suf, wr); }
  }

void logsys_print_op4(FILE *wr, char *pre, logsys_op4_t op, char *suf)
  {
    if ((pre != NULL) && ((*pre) != 0)) { fputs(pre, wr); }
    uint32_t uop = op;
    fprintf(wr, "o%04X", uop);
    if ((suf != NULL) && ((*suf) != 0)) { fputs(suf, wr); }
  }

void logsys_print_op5(FILE *wr, char *pre, logsys_op5_t op, char *suf)
  {
    if ((pre != NULL) && ((*pre) != 0)) { fputs(pre, wr); }
    uint32_t uop = op;
    fprintf(wr, "o%08X", uop);
    if ((suf != NULL) && ((*suf) != 0)) { fputs(suf, wr); }
  }

void logsys_print_op6(FILE *wr, char *pre, logsys_op6_t op, char *suf)
  {
    if ((pre != NULL) && ((*pre) != 0)) { fputs(pre, wr); }
    uint64_t uop = op;
    fprintf(wr, "o%016" uint64_hex_fmt "", uop);
    if ((suf != NULL) && ((*suf) != 0)) { fputs(suf, wr); }
  }

void logsys_print_ops(FILE *wr, char *pre, int n, int m, uint64_t op[], char *suf)
  {
    if ((pre != NULL) && ((*pre) != 0)) { fputs(pre, wr); }
    int i, j;
    /* Print the symbolic or hexadecimal opcodes: */
    for (i = 0; i < m; i++) { logsys_print_op(wr, "  ", n, op[i], FALSE, NULL); }
    fprintf(wr, "\n");
    uint64_t NX = (ONE64 << (uint64_t)n); /* Number of argument assignments. */
    uint64_t X;
    for (X = 0; X < NX; X++)
      { bool_t x[n];
        logsys_get_bits(X, n, x);
        fprintf(wr, "  ");
        for (j = n-1; j >= 0; j--)
          { fprintf(wr, "%c", (x[j] ? '1' : '0')); }
        for (i = 0; i < m; i++)
          { fprintf(wr, " %c", ((op[i] & (ONE64 << X)) != ZERO64 ? '1' : '0')); }
        fprintf(wr, "\n");
      }
    fflush(wr); 
    if ((suf != NULL) && ((*suf) != 0)) { fputs(suf, wr); }
  }

void logsys_print_op(FILE *wr, char *pre, int n, uint64_t op, bool_t verbose, char *suf)
  {
    demand(n <= logsys_n_MAX, "invalid n");
    if (verbose)
      { logsys_print_ops(wr, pre, n, 1, &(op), suf); }
    else
      { if ((pre != NULL) && ((*pre) != 0)) { fputs(pre, wr); }
        int NX = 1 << n;
        uint64_t nmsk = ALL64 + (n == 6 ? ZERO64 : (ONE64 << (uint64_t)NX));   
        if (op == ZERO64)
          { /* Equation that always fails: */
            fputs("FAL", wr);
          }
        else if (op == nmsk)
          { /* Equation that always succeeds: */
            fputs("TRU", wr);
          }
        else
          { switch(n)
              { case 0: logsys_print_op0(wr, NULL, op, NULL); break;
                case 1: logsys_print_op1(wr, NULL, op, NULL); break;
                case 2: logsys_print_op2(wr, NULL, op, NULL); break;
                case 3: logsys_print_op3(wr, NULL, op, NULL); break;
                case 4: logsys_print_op4(wr, NULL, op, NULL); break;
                case 5: logsys_print_op5(wr, NULL, op, NULL); break;
                case 6: logsys_print_op6(wr, NULL, op, NULL); break;
                default: demand(FALSE, "bad arg count"); 
              }
          }
        if ((suf != NULL) && ((*suf) != 0)) { fputs(suf, wr); }
      }
  }

void logsys_print_assignments(FILE *wr, char *pre, int n, bool_t x[], char *suf)
  {
    if ((pre != NULL) && ((*pre) != 0)) { fputs(pre, wr); }
    int i;
    uint64_t V = 0; 
    for (i = n-1; i >= 0; i--)
      { fputc("01"[x[i]], wr); 
        if (x[i]) { V |= (ONE64 << (uint64_t)i); }
      }
    if (n <= 64)
      { uint64_t Vmax = (ONE64 << (uint64_t)n) + ALL64;
        int ndigs = digits(Vmax); 
        fprintf(wr, " = %0*" uint64_dec_fmt "", ndigs, V);
        int nhex = (n + 3)/4;
        fprintf(wr, " = %0*" uint64_hex_fmt "", nhex, V);
      }
    if ((suf != NULL) && ((*suf) != 0)) { fputs(suf, wr); }
  }

void logsys_print_var_assignments(FILE *wr, char *pre, int n, logsys_va_t *va[], bool_t x[], char *suf)
  {
    int m = 0;   /* Number of non-NULL variables. */
    bool_t y[n]; /* The bits to print are {y[0..m-1]}. */
    int i;
    for (i = 0; i < n; i++) 
      { if (va[i] != NULL)
          { logsys_va_id_t id = va[i]->id;
            y[m] = x[id];
            m++; 
          }
      }
    logsys_print_assignments(wr, pre, m, y, suf);
  }

void logsys_print_ranges(FILE *wr, char *pre, int n, bool_t lo[], bool_t hi[], char *suf)
  {
    if ((pre != NULL) && ((*pre) != 0)) { fputs(pre, wr); }
    uint64_t V = 0; 
    bool_t def = TRUE; /* TRUE iff {lo..hi} is singleton for all bits. */
    int i;
    for (i = n-1; i >= 0; i--) 
      { char c;
        V = V << ONE64;
        if (lo[i] < hi[i]) 
          { c = '*'; def = FALSE; }
        else if (lo[i] > hi[i])
          { c = '!'; def = FALSE; }
        else if (lo[i])
          { c = '1'; V |= ONE64; }
        else
          { c = '0'; }
        fputc(c, wr);
      }
    if ((n <= 64) && def)
      { uint64_t Vmax = (ONE64 << (uint64_t)n) + ALL64;
        int ndigs = digits(Vmax); 
        fprintf(wr, " = %0*" uint64_dec_fmt "", ndigs, V);
        int nhex = (n + 3)/4;
        fprintf(wr, " = %0*" uint64_hex_fmt "", nhex, V);
      }
    if ((suf != NULL) && ((*suf) != 0)) { fputs(suf, wr); }
  }

void logsys_print_var_ranges(FILE *wr, char *pre, int n, logsys_va_t *va[], bool_t lo[], bool_t hi[], char *suf)
  {
    int m = 0;   /* Number of non-NULL variables. */
    bool_t loy[n], hiy[n]; /* The ranges to print are {loy[0..m-1] .. hiy[0..m-1]}}. */
    int i;
    for (i = 0; i < n; i++) 
      { if (va[i] != NULL)
          { logsys_va_id_t id = va[i]->id;
            loy[m] = lo[id];
            hiy[m] = hi[id];
            m++; 
          }
      }
    logsys_print_ranges(wr, pre, m, loy, hiy, suf);
  }
  
void logsys_print_arg_ranges(FILE *wr, char *pre, logsys_eq_t *eq, bool_t lo[], bool_t hi[], char *suf)
  {
    int n = eq->n;
    logsys_va_t *va[n];
    int j;
    for (j = 0; j < n; j++) { va[j] = eq->arg[j].va; }
    logsys_print_var_ranges(wr, pre, n, va, lo, hi, suf);
  }

void logsys_print_arg_set
  ( FILE *wr, 
    char *pre, 
    char *vtag_old, 
    char *vtag_new, 
    int n, 
    int k[],
    bool_t x[], 
    char *sep, 
    char *suf
  )
  { 
    if ((pre != NULL) && ((*pre) != 0)) { fputs(pre, wr); }
    char *cur_sep = "";
    int i;
    for (i = n-1; i >= 0; i--) 
      { char *vtag;
        int vnum;
        if (k == NULL)
          { vtag = vtag_old; vnum = i; }
        else if (k[i] >= 0)
          { vtag = vtag_old; vnum = k[i]; }
        else
          { vtag = vtag_new; vnum = i; }
        fprintf(wr, "%s%s%d", cur_sep, vtag, vnum);
        if (x != NULL) { fprintf(wr, "=%c", "01"[x[i]]); }
        cur_sep = (sep == NULL ? "" : sep); 
      }
    if ((suf != NULL) && ((*suf) != 0)) { fputs(suf, wr); }
  }

#define logsys_MAX_CALLS (ONE64 << 35LLU)
  /* Max calls to {do_find_solution}. */

void logsys_sort_variables_for_solver
  ( int heur, 
    logsys_t *S, 
    int nv, 
    logsys_va_t *va[], 
    int vix[], 
    bool_t guess[]
  )
  {
    /* Default guesses if the heuristic does not care toset them: */
    if (guess != NULL)
      { /* Suggest a first {guess} for each variable: */
        int iv;
        for (iv = 0; iv < nv; iv++) 
          { guess[iv] = (((iv ^ (iv >> 1) ^ (iv >> 2)) & 1) == 1); } 
      }

    switch(heur)
      { 
        case 0: 
        case 1: 
          logsys_sort_variables_for_solver_H0_H1(heur, S, nv,va, vix, guess); break;
        case 2: 
          logsys_sort_variables_for_solver_H2(heur, S, nv,va, vix, guess); break;
        case 3: 
          logsys_sort_variables_for_solver_H3(heur, S, nv,va, vix, guess); break;
        case 4: 
          logsys_sort_variables_for_solver_H4(heur, S, nv,va, vix, guess); break;
        default: assert(FALSE); /* No such heuristic. */
      }
  }
    
void logsys_sort_variables_by_score(int nv, logsys_va_t *va[], double score[], int sign, int vix[])
  {
    int i;
    for (i = 0; i < nv; i++) { vix[i] = i; }
    /* Sort {vix[0..nv]} by score decreasing: */
    auto int score_cmp(const void *a, const void *b);
    qsort(vix, nv, sizeof(int), &score_cmp);
    /* Validate sort: */
    for (i = 1; i < nv; i++) 
      { if (va[vix[i]] != NULL)
          { /* Make sure that all NULLs are at the end: */
            assert(va[vix[i-1]] != NULL); 
            /* Check score ordering of non-NULL entries: */
            assert(sign*score[vix[i-1]] <= sign*score[vix[i]]); 
          }
      }
    return;
    
    
    /* INTERNAL IMPLS */
    
    int score_cmp(const void *aP, const void *bP)
      { int ia = *((int *)aP);
        int ib = *((int *)bP);
        bool_t nulla = (va[ia] == NULL);
        bool_t nullb = (va[ib] == NULL);
        
        if ((!nulla) && nullb)
          { return -1; }
        else if (nulla && (!nullb))
          { return +1; }
        else if (nulla && nullb)
          { return 0; }
        else if (sign*score[ia] < sign*score[ib])
          { return -1; }
        else if (sign*score[ia] > sign*score[ib])
          { return +1; }
        else
          { return 0; }
      }
  }

bool_t logsys_find_solution(logsys_t *S, int nv, logsys_va_t *va[], int vix[], bool_t lo[], bool_t hi[], bool_t guess[])
  {
    bool_t verbose = TRUE;
    
     /* Let's say that a variable {ua} is /fixed/ if its range 
       {lo[ua->id]..hi[ua->id]}, whether  set by the client or by this procedure,
       is a single value. */
    
    /* The list of currently fixed variables: */
    int nf = 0; /* Count of fixed variables. */
    logsys_va_t *f[nv]; /* The fixed variables of {S} are {f[0..nf-1]}. */
    
    auto bool_t gather_fixed_variables(void);
      /* Gathers all variables in {va[0..nv-1]} that are effectively
        fixed (that is, that have {lo[i] == hi[i]}, and puts them in
        the queue {f}. Assumes that {nf} is zero at start. If any
        non-null variable has empty range, returns FALSE, else returns
        TRUE. */
    
    auto bool_t scan_all_equations(void);
      /* Scans all equations in the system, checking their admissible
        solutions given the ranges {lo[0..nv-1]} and {hi[0..nv-1]}.
        
        If the procedure detects that the system has no solution 
        with the given range, it returns FALSE, leaving {lo}
        and {hi} unchanged.  Otherwise it returns TRUE,
        after appending to the queue {f[0..nf-1]} any variables
        which have been fixed by the scan. */
    
    auto void unfix_variables(int nf0);
      /* Pops all variables {f[nf0..nf-1]} from the queue {f},
        restores their ranges to indeterminate ({FALSE..TRUE}),
        leaving {nf} equal to {nf0}. */
        
    int last_var_found = -1; /* ID of last variable found by {find_unfixed_variable} */ 
    
    auto logsys_va_t *find_unfixed_variable(void);
      /* Returns a variable {ua} from the list {va[0..nv]} which is
        currently not fixed. If all variables are fixed, returns NULL.
        If {vix} was specified, uses it to search for {ua}; otherwise
        starts the search after the last variable found. */
        
    auto bool_t apply_equation(logsys_eq_t *eq);
      /* Evaluates the effect of {eq} on the ranges of its argument
        variables, given their current ranges.
        
        If the equation cannot be satified at all, the procedure
        returns FALSE without modifying {lo} or {hi}; otherwise it
        returns TRUE. In this case, if the equation implies range for
        any variable that is smaller than the current one (i.e.
        fixed), the procedure updates {lo} and {hi} accordingly, and
        appends that variable to the list {f[0..nf]}. */ 
    
    auto bool_t propagate_values(int nc);
      /* Narrows the range {lo[i]..hi[i]} of each variable {va[i]} as
        implied by each equation of {S}, considered independently of
        any other equations. 
        
        Assumes that all fixed variables are in {f[0..nf-1]}, but
        equations that involve only variables {f[0..nc-1]} (the /cold/
        variables) are already as tight as possible. So
        the procedure will only check equations that involve variables
        {f[nc..nf-1]}, or any variables that become fixed during 
        the call (the /warm/ variables).
        
        If the procedure discovers that no assignment {x[0..nv-1]} with
        {lo[i]<=x[i]<=hi[i]} can satisfy the system, it returns FALSE
        leaving the vectors {lo[0..nv-1]} and {hi[0..nv-1]} unchanged.
        Otherwise it returns TRUE after narrowing the range
        {lo[i]..hi[i]} of each variable as much as possible.
        Note that a TRUE result does not mean that the system
        has a solution, but only that value propagation cannot
        narrow the ranges any further.
        
        Let {nf0} be the value of {nf} upon entry, and {nq1} be its
        value upon exit. If the procedure returns FAlSE, then {nq1 ==
        nq0}. If it returns TRUE, then {nq1 >= nq0}, and
        {f[nq0..nq1-1]} will be all the variables which had its range
        narrowed by this procedure. In either case, {nc} will be
        irrelevant upon exit. */
    
    auto bool_t do_find_solution(int nc, int lev);
      /* Tries to find a solution within the current
        ranges {lo[i]..hi[i]}.
        
        Assumes that all fixed variables are in {f[0..nf-1]}, but
        equations that involve only variables {f[0..nc-1]} (the /cold/
        variables) are already as tight as possible. So
        the procedure will only check equations that involve variables
        {f[nc..nf-1]} or any variables that become fixed 
        during the call (the /warm/ variables).
        
        If the procedure finds an assignment {x[0..nv-1]} with
        {lo[i]<=x[i]<=hi[i]} that satisfies the system, it returns
        TRUE after setting {lo[i] == hi[i] == x[i]}. Otherwise, the
        system has no solution; the procedure then returns FALSE,
        leaving the vectors {lo[0..nv-1]} and {hi[0..nv-1]} as they
        were upon entry.
        
        Let {nf0} be the value of {nf} upon entry, and {nq1} be its
        value upon exit. If the procedure returns FALSE, then {nq1 ==
        nq0}. If it returns TRUE, then {nq1 >= nq0}, and
        {f[nq0..nq1-1]} will be all the variables which had its range
        narrowed by this procedure. In either case, {nc} will be
        irrelevant upon exit.  */
    
    /* Some recursion statistics for {do_find_solution}: */
    uint64_t ncalls = 0;
    int maxlev = -1;
    
    bool_t gather_fixed_variables(void)
      { assert(nf == 0);
        logsys_va_id_t iv;
        for (iv = 0; iv < nv; iv++)
          { logsys_va_t *ua = va[iv]; 
            if (ua != NULL) 
              { assert(ua->id == iv);
                if (lo[iv] > hi[iv]) 
                  { /* Given range is empty: */
                    if (verbose) { fprintf(stderr, "(%s: v%04" uint64_dec_fmt " has empty range)", __FUNCTION__, iv); }
                    return FALSE;
                  }
                if (lo[iv] == hi[iv]) { f[nf] = ua; nf++; }
              }
          }
        return TRUE;
      }

    bool_t scan_all_equations(void)
      { /* Save the queue state for undoing fixes in case of failure: */
        int nf0 = nf;
        /* Scan all equations of the system: */
        if (S->eq != NULL)
          { logsys_eq_t *eq = S->eq->nxeq;
            do {
              bool_t ok = apply_equation(eq);
              if (! ok) 
                { /* System has no solution.  Restore {lo,hi} and return falure: */
                  if (verbose) { fprintf(stderr, "(%s: e%04" uint64_dec_fmt " is not satisfiable)", __FUNCTION__, eq->id); }
                  unfix_variables(nf0);
                  return FALSE;
                }
              eq = eq->nxeq;
            } while (eq != S->eq->nxeq);
          }
        return TRUE;
      }

    void unfix_variables(int nf0)
      { while (nf > nf0)
          { nf--; 
            logsys_va_t *wa = f[nf];
            assert(wa != NULL);
            logsys_va_id_t iw = wa->id;
            lo[iw] = FALSE; hi[iw] = TRUE; 
          }
      }

    logsys_va_t *find_unfixed_variable(void)
      { /* !!! This should be made more efficient !!! */
        /* Starting index for search when {vix} is NULL: */
        int start = (vix != NULL ? 0 : last_var_found + 1);
        int r;
        for (r = 0; r < nv; r++)
          { /* Get the next candidate ID {iu} in the search order: */
            logsys_va_id_t iu = (vix == NULL ? (r + start) % nv : vix[r]);
            assert((iu >= 0) && (iu < nv));
            logsys_va_t *ua = va[iu]; 
            if (ua != NULL) 
              { /* If this variable is non-null and unfixed, return it: */
                assert(ua->id == iu);
                if (lo[iu] < hi[iu]) 
                  { last_var_found = iu; 
                    return ua;
                  }
              }
            else
              { /* If {vix} is not NULL, we must find an unfixed var before the first NULL: */
                assert(vix == NULL);
              }
          }
        /* No unfixed variables remain: */
        return NULL;
      }
      
    bool_t apply_equation(logsys_eq_t *eq)
      {
        bool_t debug = FALSE /* (eq->id == 10) */;
        if (debug) { fprintf(stderr, "(e%04" uint64_dec_fmt "", eq->id); logsys_print_arg_ranges(stderr, NULL, eq, lo, hi, NULL); }
        /* Extract the current ranges of the arguments: */
        int n = eq->n;
        bool_t lop[n], hip[n];
        int j;
        for (j = 0; j < n; j++) 
          { logsys_va_t *ua = eq->arg[j].va;
            assert(ua != NULL);
            logsys_va_id_t iu = ua->id;
            lop[j] = lo[iu];
            hip[j] = hi[iu];
          }
        /* Contract the ranges according to {eq}: */
        bool_t sat = logsys_narrow_ranges_by_equation(eq, lop, hip);
        if (! sat) { /* Equation cannot be satisfied: */ return FALSE; }
        /* Update the global ranges of the argument variables: */
        for (j = 0; j < n; j++) 
          { logsys_va_t *ua = eq->arg[j].va;
            assert(ua != NULL);
            logsys_va_id_t iu = ua->id; 
            assert(lop[j] >= lo[iu]);
            assert(hip[j] <= hi[iu]);
            assert(lop[j] <= hip[j]); /* Since {logsys_narrow_ranges_by_equation} succeeded. */
            if ((lop[j] > lo[iu]) || (hip[j] < hi[iu]))
              { /* The range of {ua} was reduced by {eq}. */
                /* Fix the variable {ua}: */
                assert(lop[j] == hip[j]);
                lo[iu] = lop[j]; hi[iu] = hip[j];
                /* Insert {ua} into the queue: */
                /* Note that this step can be executed at most once for each variable. */
                assert(nf < nv);
                f[nf] = ua; nf++;
              }
          }
        if(debug) { logsys_print_arg_ranges(stderr, " --> ", eq, lo, hi, ")\n"); }
        return TRUE;
      }
    
    bool_t propagate_values(int nc)
      {
        /* Save current queue pointer for backtracking: */
        int nf0 = nf;
        
        /* Loop on queued variables: */
        while (nc < nf)
          { /* Get a warm variable {ua} from {f}, make it cold: */;
            logsys_va_t *ua = f[nc]; nc++; 
            assert(ua != NULL);
            bool_t debug = FALSE /* ((ua->id == 2) || (ua->id == 8) || (ua->id == 20)) */; 
            if (debug) { logsys_print_variable_id(stderr, "[v", ua->id, 4, "]"); }
            /* Enumerate the equations that use {ua}: */
            logsys_eq_t *eq = ua->use; 
            if (eq != NULL)
              { do
                  { if (debug) { fputc('.', stderr); } 
                    bool_t okeq = apply_equation(eq); 
                    if (! okeq)
                      { /* The equation cannot be satisfied with the current ranges. */
                        /* Unfix all variables that we had fixed, and return with failure: */
                        if (verbose) { fprintf(stderr, "(%s: e%04" uint64_dec_fmt " is not satisfiable)", __FUNCTION__, eq->id); }
                        unfix_variables(nf0);
                        return FALSE;
                      }
                    /* Go to the next equation: */
                    int uj = logsys_find_arg(eq, ua);
                    assert((uj >= 0) && (uj < eq->n) && (eq->arg[uj].va == ua));
                    eq = eq->arg[uj].nxeq;
                  }
                while (eq != ua->use);
              }
          }
        /* We propagated as far as possible without detecting a contradiction. */
        if (verbose) { fprintf(stderr, "(%s: fixed %d vars)", __FUNCTION__, nf - nf0); }
        return TRUE;
      }
   
    bool_t do_find_solution(int nc, int lev)
      {
        /* Save current queue pointer for backtracking: */
        int nf0 = nf;
        
        /* Update statistics: */
        demand(ncalls < logsys_MAX_CALLS, "too many recursive calls, aborted");
        ncalls++;
        if (lev > maxlev) { maxlev = lev; }
        
        if (verbose) { fprintf(stderr, "%*s%s: propagating ... ", lev, "", __FUNCTION__); }
        bool_t ok = propagate_values(nc);
        if (verbose) { fprintf(stderr, "\n"); }
        if (! ok)
          { /* Propagation led to a contradiction. */
            if (verbose) { fprintf(stderr, "%*s%s: failure\n", lev, "", __FUNCTION__); }
            assert(nf == nf0); 
            return FALSE;
          }
        else
          { /* Check whether there are any unfixed variables. */
            logsys_va_t *ua = find_unfixed_variable();
            if (ua == NULL)
              { /* Success! Return TRUE, leaving all variables fixed. */ 
                if (verbose) { fprintf(stderr, "%*s%s: no more unfixed variables\n", lev, "", __FUNCTION__); }
                return TRUE;
              }
            else
              { /* Must try both values and recurse: */
                logsys_va_id_t iu = ua->id;
                if (verbose) { fprintf(stderr, "%*s%s: variable v%04" uint64_dec_fmt " still undef\n", lev, "", __FUNCTION__, ua->id); }
                /* Choose a random boolean value {val} to start. */
                bool_t val;
                if (guess != NULL)
                  { /* Use the suggested first guess: */
                    val = guess[iu];
                  }
                else
                  { /* Pick a pseudorandom boolean: */
                    val = ((((iu >> 2) ^(iu >> 1) ^ iu ^ ncalls) & 1) == 0); 
                  }
                /* Stack {ua} in {f} as a warm variable: */
                f[nf] = ua; nf++;
                int try;
                for (try = 0; try < 2; try++)
                  { 
                    /* Try setting {ua} to {val}: */
                    if (verbose) { fprintf(stderr, "%*s%s: trying v%04" uint64_dec_fmt " = %d\n", lev, "", __FUNCTION__, ua->id, val); }
                    lo[iu] = hi[iu] = val; 
                    if (do_find_solution(nf-1, lev+1)) 
                      { /* Success: */
                        if (verbose) { fprintf(stderr, "%*s%s: success\n", lev, "", __FUNCTION__); }
                        return TRUE; 
                      }
                    /* Setting {ua} to {val} failed. */
                    /* Variable {ua} must still be at the end of the queue: */
                    assert(nf > 0);
                    assert(f[nf-1] == ua);
                    /* Try setting {ua} to the opposite value: */
                    val = (! val);
                  }
                if (verbose) { fprintf(stderr, "%*s%s: failure\n", lev, "", __FUNCTION__); }
                /* Unfix all variables that were fixed in the call, and return with failure: */
                unfix_variables(nf0);
                return FALSE;
              }
          }
      }
      
    bool_t report = TRUE; /* TRUE to print the top-level diagnostics. */
    
    /* Find all fixed variables, put them in {f}: */
    if (report) { fprintf(stderr, "%s: gathering variables fixed by given ranges ... ", __FUNCTION__); } 
    bool_t ok1 = gather_fixed_variables();
    if (ok1)
      { if (report) { fprintf(stderr, " found %d fixed vars\n", nf); } }
    else
      { /* The given range is empty: */
        if (report) { fprintf(stderr, " given range is empty\n"); }
        return FALSE;
      }
    
    /* Scan all equations (some of them may fix variables spontaneously): */
    if (report) { fprintf(stderr, "%s: initial equation scan ... ", __FUNCTION__); } 
    int nf1 = nf;
    bool_t ok2 = scan_all_equations();
    if (ok2)
      { if (report) { fprintf(stderr, " fixed %d vars\n", nf - nf1); } }
    else
      { /* The given range is incompatible with some equation: */
        if (report) { fprintf(stderr, " scan failed\n"); }
        return FALSE;
      }
    
    /* Main call (root of recursion tree): */
    int nf2 = nf;
    bool_t success = do_find_solution(0, 0);
    if (report) { fprintf(stderr, "%s: recursion calls = %" uint64_dec_fmt "  max depth = %d\n", __FUNCTION__, ncalls, maxlev); }
    if (success)
      { if (report) { fprintf(stderr, "%s: recursion succeded and fixed %d vars\n", __FUNCTION__, nf - nf2); } 
        /* Do an extra check of all equations, just to be sure: */
        if (report) { fprintf(stderr, "%s: final eq scan ... ", __FUNCTION__); } 
        int nf3 = nf;
        bool_t ok3 = scan_all_equations();
        if (ok3) 
          { if (report) { fprintf(stderr, " fixed %d vars\n", nf - nf3); } }
        else
          { fprintf(stderr, " ** failed!\n"); 
            demand(FALSE, "aborted");
          }
      }
    else
      { if (report) { fprintf(stderr, "%s: recursion failed to find solution\n", __FUNCTION__); } }
      
    return success;
  }

int logsys_get_variable_table_size(logsys_t *S)
  { logsys_va_id_t idmax = 0;
    if (S->va != NULL)
      { logsys_va_t *ua = S->va;
        do {
          logsys_va_id_t iu = ua->id;
          if (iu >= idmax) { idmax = iu; }
          ua = ua->nxva;
        } while (ua != S->va); 
      }
    assert(idmax < MAXINT);
    return 1 + (int)idmax;
  }
  
int logsys_get_equation_table_size(logsys_t *S)
  { logsys_eq_id_t idmax = 0;
    if (S->eq != NULL)
      { logsys_eq_t *eq = S->eq;
        do {
          logsys_eq_id_t ie = eq->id;
          if (ie >= idmax) { idmax = ie; }
          eq = eq->nxeq;
        } while (eq != S->eq); 
      }
    assert(idmax < MAXINT);
    return 1 + (int)idmax;
  }

void logsys_get_variables(logsys_t *S, int *nvP, logsys_va_t ***vaP)
  { int nv = logsys_get_variable_table_size(S);
    logsys_va_t **va;
    if (nv == 0)
      { va = NULL; }
    else
      { va = (logsys_va_t **)notnull(malloc(nv*sizeof(logsys_va_t *)), "out of mem");
        logsys_va_id_t iu; 
        for (iu = 0; iu < nv; iu++) { va[iu] = NULL; }
        logsys_va_t *ua = S->va;
        do {
          iu = ua->id;
          assert((iu >= 0) && (iu < nv));
          va[iu] = ua;
          ua = ua->nxva;
        } while (ua != S->va);
      }
    (*vaP) = va;
    (*nvP) = nv;
  }
  
void logsys_get_equations(logsys_t *S, int *neP, logsys_eq_t ***eqP)
  { int ne = logsys_get_equation_table_size(S);
    logsys_eq_t **eq;
    if (ne == 0)
      { eq = NULL; }
    else
      { eq = (logsys_eq_t **)notnull(malloc(ne*sizeof(logsys_eq_t *)), "out of mem");
        logsys_eq_id_t ie; 
        for (ie = 0; ie < ne; ie++) { eq[ie] = NULL; }
        logsys_eq_t *fq = S->eq;
        do {
          ie = fq->id;
          assert((ie >= 0) && (ie < ne));
          eq[ie] = fq;
          fq = fq->nxeq;
        } while (fq != S->eq);
      }
    (*eqP) = eq;
    (*neP) = ne;
  }

bool_t logsys_equation_is_satisfied(logsys_eq_t *eq, int nv, logsys_va_t *va[], bool_t x[])
  {
    int n = eq->n;
    
    /* Extract the values {xp[0..n-1]} of the arguments from {x}, pack as integer {X}: */
    int j;
    uint64_t X = 0;
    for (j = n-1; j >= 0; j--) 
      { logsys_va_t *ua = eq->arg[j].va;
        assert(ua != NULL);
        logsys_va_id_t iu = ua->id;
        assert(va[iu] == ua);
        X = X << ONE64;
        if (x[iu]) { X = X | ONE64; }
      }
    /* Check whether {X} satisfies {eq->op}: */
    return ((eq->op & (ONE64 << X)) != ZERO64); 
  }

bool_t logsys_is_solution(logsys_t *S, int nv, logsys_va_t *va[], bool_t x[])
  {
    if (S->eq == NULL) { return TRUE; }
    logsys_eq_t *eq = S->eq;
    do
      {
        if (! logsys_equation_is_satisfied(eq, nv, va, x)) { return FALSE; }
        eq = eq->nxeq;
      }
    while (eq != S->eq);
    return TRUE;
  }

void logsys_enum_solutions(logsys_t *S, int nv, logsys_va_t *va[], bool_t lo[], bool_t hi[], logsys_solution_proc_t *proc)
  {
    bool_t x[nv];
    int iv;
    for (iv = 0; iv < nv; iv++) 
      { bool_t lov = (lo != NULL ? lo[iv] : FALSE);
        bool_t hiv = (hi != NULL ? hi[iv] : TRUE);
        if (lov > hiv) { fprintf(stderr, "!! empty range for variable"); return; }
        x[iv] = lov;
      }
    do
      { bool_t is_sol = logsys_is_solution(S, nv, va, x);
        if (is_sol) 
          { bool_t ok = proc(nv, va, x);
            if (! ok) { return; }
          }
        /* Get next assignment: */
        iv = nv - 1;
        while (iv >= 0)
          { if (va[iv] != NULL)
              { bool_t lov = (lo != NULL ? lo[iv] : FALSE);
                bool_t hiv = (hi != NULL ? hi[iv] : TRUE);
                if (x[iv] > lov)
                  { x[iv] = lov; }
                else if (x[iv] < hiv)
                  { x[iv] = hiv; break; }
                else
                  { /* Leave {x[iv]} unchanged. */ }
              }
            iv--;
          }
      }
    while (iv >= 0);
  }

bool_t logsys_narrow_ranges_by_equation(logsys_eq_t *eq, bool_t lop[], bool_t hip[])
  {
    int n = eq->n;

    /* Tight ranges for the set {A}, initially empty: */
    int j, k;
    bool_t lox[n], hix[n];
    for (j = 0; j < n; j++) { lox[j] = TRUE;  hix[j] = FALSE; }
            
    /* Enumerates all assignments allowed by {op}: */
    assert(n < 31);
    uint64_t NX = ONE64 << (uint64_t)n; /* Number of assignments for {n} variables: */
    assert (NX <= logsys_NX_MAX);
    uint64_t X;
    for (X = 0; X < NX; X++)
      { /* See if the bits of {X} satisfy the opcode {op}: */
        if ((eq->op & (ONE64 << X)) != ZERO64)
          { /* X satisfies the opcode, unpack {X} as a Boolean vector {x[0..n-1]}: */
            bool_t x[n];
            logsys_get_bits(X, n, x);
            /* Check whether it is compatible with the range and aliasing: */
            bool_t ok = TRUE; /* A priori, but... */
            for (j = 0; (j < n) && ok; j++)
              { logsys_va_t *va = eq->arg[j].va;
                assert(va != NULL);
                for (k = 0; k < j; k++) { assert(eq->arg[k].va != va); }
                /* Check whether {x[j]} is in the given range: */
                if ((x[j] < lop[j]) || (x[j] > hip[j])) { ok = FALSE; }
              }
            if (ok)
              { /* The assignment {x[0..n-1]} is in {A}. */
                /* Expand the new range to include it: */
                for (j = 0; j < n; j++) 
                  { if (x[j] < lox[j]) { lox[j] = x[j]; }
                    if (x[j] > hix[j]) { hix[j] = x[j]; }
                  }
              }
          }
      }
    /* Now {lox[0..n-1], hix[0..n-1]} is the smallest range that encloses {A}. */
    /* Update the given ranges {lop,hip}: */
    bool_t satisfiable = TRUE; /* A priori. */
    for(j = 0; j < n; j++)
      { if (lox[j] > hix[j]) { satisfiable = FALSE; }
        lop[j] = lox[j]; hip[j] = hix[j];
      }
    return satisfiable;
  }

void logsys_get_bits(uint64_t X, int n, bool_t x[])
  {
    demand(n <= 64, "not enough bits in X");
    int i;
    for (i = 0; i < n; i++) 
      { x[i] = ((X & ONE64) != ZERO64); X >>= ONE64; }
    demand(X == ZERO64, "** arg has more than {n} bits");
  }


logsys_va_id_t logsys_variable_id(logsys_va_t *va)
  { 
    return va->id;
  }

logsys_eq_id_t logsys_equation_id(logsys_eq_t *eq)
  { 
    return eq->id;
  }

void logsys_debug_scramble(bool_t enter, bool_t scram, logsys_op_t op_in, int n_old, int n_new, int k[], logsys_op_t op_out);

void logsys_debug_scramble(bool_t enter, bool_t scram, logsys_op_t op_in, int n_old, int n_new, int k[], logsys_op_t op_out)
  { 
    int indent = 1 + 2*logsys_n_MAX - (n_old + n_new);
    char *fun = (scram ? "scramble" : "unscramble");
    char *tag = (enter ? "+" : "-");
    int nhex_in = 1 << imax((scram ? n_old : n_new) - 2, 0);
    int nhex_out = 1 << imax((scram ? n_new : n_old) - 2, 0);
    fprintf(stderr, "%*s%s %s(%0*" uint64_hex_fmt ",%d,%d,[", indent,"", tag, fun, nhex_in,op_in, n_old, n_new);
    int j;
    char *sep = "";
    for (j = 0; j < n_new; j++) { fprintf(stderr, "%s%d", sep, k[j]); sep = ","; }
    fprintf(stderr, "])");
    if (! enter) { fprintf(stderr, " = %0*" uint64_hex_fmt "", nhex_out,op_out); }
    fprintf(stderr, "\n");
  }

logsys_op_t logsys_scramble_args_in_op(logsys_op_t op_old, int n_old, int n_new, int k[])
  { 
    bool_t debug = FALSE;
    if (debug) { logsys_debug_scramble(TRUE, TRUE, op_old, n_old, n_new, k, 0); }
    logsys_op_t op_new;
    if (op_old == 0)
      { /* Old equation was FAIL, return FAIL: */
        op_new =  (logsys_op_t)0;
      }
    else if (n_new == 0)
      { /* There are no new arguments, and equation was not FAIL: */
        op_new =  (logsys_op_t)ONE64;
      }
    else if (k[n_new-1] < 0)
      { /* The last new argument is a new variable: */
        logsys_op_t op_tmp = logsys_scramble_args_in_op(op_old, n_old, n_new-1, k);
        int HY = (1 << (n_new-1)); /* Number of assignments for {n_new-1} args. */
        assert ((uint64_t)HY <= logsys_NX_MAX);
        op_new =  op_tmp | (op_tmp << (uint64_t)HY);
      }
    else 
      { /* The last new argument is some old variable: */
        int i_old = k[n_new-1];
        assert(i_old < n_old);
        int IX = (1 << i_old); /* Number of assignments for {i_old} args. */
        /* Masks that select in {op_old} the outcomes where {k[n-1]} is set to 0: */
        logsys_op_t omsk_0 = (ONE64 << (uint64_t)IX) + ALL64;
        int SX = 2*IX; while ((uint64_t)SX < logsys_NX_MAX) { omsk_0 = omsk_0 | (omsk_0 << (uint64_t)SX); SX = SX << 1; }
        /* Get old opcodes {op_old_0,op_old_1} assuming old arg {k[n-1]} is set to {0,1}: */
        logsys_op_t op_old_0 = op_old & omsk_0;
        logsys_op_t op_old_1 = (op_old >> (uint64_t)IX) & omsk_0;
        /* Get new opcodes {op_new_0,op_new_1} for {k[0..n-2]} assuming old arg {k[n-1]} is set to {0,1}: */
        logsys_op_t op_new_0 = logsys_scramble_args_in_op(op_old_0, n_old, n_new-1, k);
        logsys_op_t op_new_1 = logsys_scramble_args_in_op(op_old_1, n_old, n_new-1, k);
        /* Join the two partial opcodes: */
        int HY = (1 << (n_new-1)); /* Number of assignments for {n_new-1} args. */
        assert ((uint64_t)HY <= logsys_NX_MAX);
        op_new =  op_new_0 | (op_new_1 << (uint64_t)HY);
      }
    if (debug) { logsys_debug_scramble(FALSE, TRUE, op_old, n_old, n_new, k, op_new); }
    return op_new;
  }

logsys_op_t logsys_unscramble_args_in_op(logsys_op_t op_new, int n_old, int n_new, int k[])
  { 
    bool_t debug = FALSE;
    if (debug) { logsys_debug_scramble(TRUE, FALSE, op_new, n_old, n_new, k, 0); }
    logsys_op_t op_old;
    if (op_new == 0)
      { /* New equation is FAIL, return FAIL: */
        op_old =  (logsys_op_t)0;
      }
    else if (n_new == 0)
      { /* There are no new arguments, and the new equation is not FAIL, so it is trivial: */
        int NX = (1 << n_old); /* Number of assignments for {n_old} args. */
        assert ((uint64_t)NX <= logsys_NX_MAX);
        op_old = ALL64 + (n_old == 6 ? ZERO64 : (ONE64 << (uint64_t)NX));
      }
    else if (k[n_new-1] < 0)
      { /* The last new argument is a new variable: */
        assert(! logsys_op_depends(op_new, n_new, n_new-1));
        int HY = (1 << (n_new-1)); /* Number of assignments for {n_new-1} args. */
        assert ((uint64_t)HY <= logsys_NX_MAX);
        uint64_t msk = ((ONE64 << (uint64_t)HY) + ALL64); /* Mask for {op_new} cases where new arg {n_new-1} is 0. */
        op_old =  logsys_unscramble_args_in_op(op_new & msk, n_old, n_new-1, k);
      }
    else 
      { /* The last new argument is some old variable: */
        int HY = (1 << (n_new-1)); /* Number of assignments for {n_new-1} args. */
        assert (HY < logsys_NX_MAX);
        uint64_t hmsk = ((ONE64 << (uint64_t)HY) + ALL64); /* Mask for op on {n_new-1} args. */
        /* Get new opcodes {op_new_0,op_new_1} when new arg {n-1} is set to {0,1}: */
        logsys_op_t op_new_0 = op_new & hmsk;
        logsys_op_t op_new_1 = (op_new >> (uint64_t)HY) & hmsk;
        /* Get a mask {omsk_0,omks_1} that selects in {op_old} the outcomes where old arg {k[n-1]} is set to {0,1}: */
        int i_old = k[n_new-1];
        assert(i_old < n_old);
        int IX = (1 << i_old); /* Number of assignments for {i_old} args. */
        logsys_op_t omsk_0 = (ONE64 << (uint64_t)IX) + ALL64; 
        int SX = 2*IX; while ((uint64_t)SX < logsys_NX_MAX) { omsk_0 = omsk_0 | (omsk_0 << (uint64_t)SX); SX = SX << 1; }
        logsys_op_t omsk_1 = ALL64 ^ omsk_0;
        /* Get old opcodes {op_old_0,op_old_1} for the cases where old arg {k[n-1]} is set to {0,1}: */
        logsys_op_t op_old_0 = logsys_unscramble_args_in_op(op_new_0, n_old, n_new-1, k);;
        logsys_op_t op_old_1 = logsys_unscramble_args_in_op(op_new_1, n_old, n_new-1, k);
        /* Merge them: */
        op_old = (omsk_0 & op_old_0) | (omsk_1 & op_old_1);
        if (debug) { fprintf(stderr, "  i_old = %d IX = %d\n", i_old, IX); }
        if (debug) { fprintf(stderr, "  omsk_0   = %016" uint64_hex_fmt " omsk_1   = %016" uint64_hex_fmt "\n", omsk_0, omsk_1); }
        if (debug) { fprintf(stderr, "  op_new_0 = %016" uint64_hex_fmt " op_new_1 = %016" uint64_hex_fmt "\n", op_new_0, op_new_1);}
        if (debug) { fprintf(stderr, "  op_old_0 = %016" uint64_hex_fmt " op_old_1 = %016" uint64_hex_fmt "\n", op_old_0, op_old_1);}
      }
    if (debug) { logsys_debug_scramble(FALSE, FALSE, op_new, n_old, n_new, k, op_old); }
    return op_old;
  }

bool_t logsys_try_to_simplify_equation(logsys_eq_t *eq)
  { 
    int n_old = eq->n;
    
    /* Determine the {m} arguments that do matter, save their indices in {k[0..n-1]} */
    int k[n_old];
    int n_new = 0;
    int j_old;
    for (j_old = 0; j_old < n_old; j_old++)
      { if (logsys_op_depends(eq->op, n_old, j_old)) 
          { k[n_new] = j_old; n_new++; }
      }
    if(n_new == n_old) { return FALSE; }
    
    /* Remove old args and put the relevant ones only: */
    int j_new = 0;
    for (j_old = 0; j_old < n_old; j_old++) 
      { logsys_va_t *va = logsys_unset_arg(eq, j_old);
        if ((j_new < n_new) && (k[j_new] == j_old)) { logsys_set_arg(eq, j_new, va); j_new++; }
      }
    assert(j_new == n_new);

    /* Build a new op for the relevant arguments only: */
    eq->op = logsys_scramble_args_in_op(eq->op, n_old, n_new, k);
    eq->n = n_new;
    return TRUE;
  }

void logsys_stack_equation(int *nstP, logsys_eq_vec_t *steq, logsys_eq_t *eq)
  { demand((nstP != NULL) == (steq != NULL), "either both null or both non-null");
    if (steq != NULL) 
      { int nst = (*nstP);
        assert((nst >= 0) & (nst <= steq->ne));
        /* Check if already stacked: */
        /* !!! Should use a boolean table to avoid scanning the stack. !!! */
        int k = nst-1;
        while ((k >= 0) && (eq != steq->e[k])) { k--; }
        if (k < 0)
          { /* Not there yet, stack it: */
            logsys_eq_vec_expand(steq, nst);
            steq->e[nst] = eq; nst++;
            (*nstP) = nst;
          }
      }
  }

void logsys_remove_equation_from_stack(int *nstP, logsys_eq_vec_t *steq, logsys_eq_t *eq)
  { demand((nstP != NULL) == (steq != NULL), "either both null or both non-null");
    if (steq != NULL) 
      { int nst = (*nstP);
        assert((nst >= 0) & (nst <= steq->ne));
        /* Remove all occurrences (should be only one): */
        /* !!! Should use a boolean table to avoid scanning the stack. !!! */
        int k = nst-1;
        int ndel = 0;
        while (k >= 0) 
          { if (eq == steq->e[k]) { steq->e[k] = steq->e[nst-1]; nst--; ndel++; }
            k--;
          }
        assert(ndel <= 1);
        (*nstP) = nst;
      }
  }

bool_t logsys_try_to_split_equation(logsys_eq_t *eq, int *nstP, logsys_eq_vec_t *steq)
  {
    bool_t debug = TRUE;
    bool_t verbose = TRUE;
    if (debug) 
      { fprintf(stderr, "        considering splitting of eq ");
        logsys_print_equation(stderr, NULL, eq, "\n");
      }
    logsys_t *S = eq->sys;
    /* Use {logsys_try_to_simplify_equation} to get rid of trivial {TRU} factors: */
    bool_t changed = logsys_try_to_simplify_equation(eq);
    /* Now we try to split the equation {eq} as much as we can: */
    while (eq->n > 0)
      { /* Try to split the opcode: */
        int n = eq->n;
        assert(n > 0);
        logsys_op_t opa, opb;
        int ka[n], kb[n];
        int na, nb;
        logsys_split_op(eq->op, eq->n, &opa, &na, ka, &opb, &nb, kb);
        if (debug) 
          { fprintf(stderr, "        factorization returned na = %d nb = %d\n", na, nb); }
        assert((na > 0) && (na <= n));
        assert((nb >= 0) && (na <= n));
        assert(na + nb == n);
        /* If the first factor is the whole eq, give up: */
        if (na == n) { break; }
        /* We found a non-trivial split of {eq->op}, split the equation: */
        changed = TRUE;
        if (debug) 
          { fprintf(stderr, "        opa = ");
            logsys_print_op(stderr, NULL, na, opa, FALSE, NULL);
            logsys_print_arg_set(stderr, "(", "x", "u", na, ka, NULL, ",", ")\n");
            fprintf(stderr, "        opb = ");
            logsys_print_op(stderr, NULL, nb, opb, FALSE, NULL);
            logsys_print_arg_set(stderr, "(", "x", "v", nb, kb, NULL, ",", ")\n");
            if (verbose)
              { logsys_op_t zopa = logsys_unscramble_args_in_op(opa, n, na, ka);
                logsys_op_t zopb = logsys_unscramble_args_in_op(opb, n, nb, kb);
                logsys_op_t zop = zopa & zopb;
                fprintf(stderr, "        operations with original argument list:\n");
                logsys_op_t prtop[4] = { eq->op, zopa, zopb, zop };
                logsys_print_ops(stderr, "        original op, factors, product:\n", n, 4, prtop, NULL);
              }
          }
        /* Get the list of variables of the two factors: */
        logsys_va_t *va[na];
        logsys_va_t *vb[nb];
        int ia, ib, j;
        for (ia = 0; ia < na; ia++) { j = ka[ia]; assert((j >= 0) && (j < n)); va[ia] = eq->arg[j].va; }
        for (ib = 0; ib < nb; ib++) { j = kb[ib]; assert((j >= 0) && (j < n)); vb[ib] = eq->arg[j].va; }
        /*  Create a separate equation {fq} for the first factor: */
        logsys_eq_t *fq = logsys_add_equation(S, opa, na, va, eq->preq);
        /* Stack it: */
        logsys_stack_equation(nstP, steq, fq);
        /* Modify the equation {eq} to be just the second factor: */
        for (j = 0; j < n; j++) { (void)logsys_unset_arg(eq, j); }
        for (ib = 0; ib < nb; ib++) { logsys_set_arg(eq, ib, vb[ib]); }
        eq->op = opb;
        eq->n = nb;
      }
    /* We have stacked all factors except for {eq} itself: */
    if (changed) { logsys_stack_equation(nstP, steq, eq); }
    return changed;
  }
  
void logsys_try_to_combine_equations(logsys_eq_t **aP, bool_t *amodP, logsys_eq_t **bP, bool_t *bmodP)
  {
    bool_t debug = TRUE;
    bool_t verbose = FALSE;
    
    demand(aP != bP, "aliased parameters");
    
    logsys_eq_t *a = (*aP);
    logsys_eq_t *b = (*bP);
    
    demand(a != b, "equations must be distinct");
    demand(a->sys == b->sys, "not in the same system");
    logsys_t *S = a->sys;
    
    demand(amodP != bmodP, "aliased parameters");
    (*amodP) = FALSE; /* For now. */
    (*bmodP) = FALSE; /* For now. */
    bool_t adel = FALSE; /* Set if node must be deleted. */
    bool_t bdel = FALSE; /* Set if node must be deleted. */
    
    /* Some subset of variables appearing in the equations: */
    int cnmax = a->n + b->n; /* Max variables in subset */
    logsys_va_t *vc[cnmax]; /* The subset is {vc[0..cn-1]}. */
    int ka[cnmax]; /* Variable {vc[j]} is argument {ka[j]} of {a}. */
    int kb[cnmax]; /* Variable {vc[j]} is argument {kb[j]} of {b}. */
    
    /* We assume that the equations are already simplified individually. */
    
    /* Extract the common variables {vc[0..sn-1]}: */
    int sn = 0; /* Count of shared variables. */
    logsys_intersect_eq_args(a, b, &sn, vc, ka, kb);
    assert(sn <= imin(a->n, b->n));
    if (sn == 0) { /* No common variables: */ return; }
    if (debug) 
      { fprintf(stderr, "        pair has %d shared variables - considering merge\n", sn);
        logsys_print_equation(stderr, "          {a} = ", a, "\n");
        logsys_print_equation(stderr, "          {b} = ", b, "\n");
      }

    /* There are common variables, worth looking further: */
    /* Obtain the valid combinations of the common variables: */
    logsys_op_t ropa = logsys_scramble_args_in_op(a->op, a->n, sn, ka);
    logsys_op_t ropb = logsys_scramble_args_in_op(b->op, b->n, sn, kb);
    logsys_op_t rop = ropa & ropb;
    if (debug) 
       { fprintf(stderr, "        opcodes rearranged for shared variables:");
         logsys_print_op(stderr, "  ropa = ", sn, ropa, FALSE, "");
         logsys_print_op(stderr, "  ropb = ", sn, ropb, FALSE, "");
         logsys_print_op(stderr, "   rop = ", sn,  rop, FALSE, "\n");
         if (verbose)
           { logsys_op_t aop[3] = { ropa, ropb, rop };
             logsys_print_ops(stderr, NULL, sn, 3, aop, NULL);
           }
       }
    if (ropa != rop)
      { /* Some cases allowed by {a} are excluded by {b}: */
        logsys_op_t xopa = logsys_unscramble_args_in_op(rop, a->n, sn, ka);
        assert((a->op & xopa) != a->op);
        a->op = a->op & xopa;
        (*amodP) = TRUE;
        int n_saved = a->n;
        (void)logsys_try_to_simplify_equation(a);
        if (debug) 
          { fprintf(stderr, "        equation {a} modified by induction");
            if (a->n < n_saved) { fprintf(stderr, " (args reduced from %d to %d)", n_saved, a->n); }
            logsys_print_op(stderr, "  xopa = ", sn, xopa, FALSE, "");
            logsys_print_op(stderr, "  a.op = ", a->n, a->op, FALSE, "\n");
          }
        if ((a->n == 0) && (a->op == 1)) { adel = TRUE; }
      }
    if (ropb != rop)
      { /* Some cases allowed by {b} are excluded by {a}: */
        logsys_op_t xopb = logsys_unscramble_args_in_op(rop, b->n, sn, kb);
        assert((b->op & xopb) != b->op);
        b->op = b->op & xopb;
        (*bmodP) = TRUE;
        int n_saved = b->n;
        (void)logsys_try_to_simplify_equation(b);
        if (debug) 
          { fprintf(stderr, "        equation {b} modified by induction");
            if (b->n < n_saved) { fprintf(stderr, " (args reduced from %d to %d)", n_saved, b->n); }
            logsys_print_op(stderr, "  xopb = ", sn, xopb, FALSE, "");
            logsys_print_op(stderr, "  b.op = ", b->n, b->op, FALSE, "\n");
          }
        if ((b->n == 0) && (b->op == 1)) { bdel = TRUE; }
      }
    
    /* Consider merging the nodes: */
    /* Note that {a->n}, {b->n} may have decreased. */
    /* Get the union {vc[0..un-1]} of all variables: */
    int un;
    logsys_unite_eq_args(a, b, &un, vc, ka, kb);
    assert(un <= cnmax);
    
    /* Update count of shared args: */
    sn = (a->n + b->n) - un; 
    if (debug) { fprintf(stderr, "        total %d variables in args (%d U %d), %d shared\n", un, a->n, b->n, sn); }
    /* Not worth merging an unary equation: */
    /* It has already had its effect on the other one, and it will immediately factor out. */
    if ((sn > 0) && (un <= logsys_n_MAX) && (a->n > 1) && (b->n > 1))
      { /* Merge is possible and perhaps profitable: */
        if (debug) { fprintf(stderr, "        attempting merge\n"); }
        logsys_op_t mopa = logsys_scramble_args_in_op(a->op, a->n, un, ka);
        logsys_op_t mopb = logsys_scramble_args_in_op(b->op, b->n, un, kb);
        logsys_op_t mop = mopa & mopb;
        if (debug) 
           { fprintf(stderr, "        opcodes rearranged for united variables");
             logsys_print_op(stderr, "  mopa = ", un, mopa, FALSE, "");
             logsys_print_op(stderr, "  mopb = ", un, mopb, FALSE, "");
             logsys_print_op(stderr, "   mop = ", un,  mop, FALSE, "\n");
             if (verbose)
               { logsys_op_t aop[3] = { mopa, mopb, mop };
                 logsys_print_ops(stderr, NULL, un, 3, aop, NULL);
               }
           }
    
        /* Mark {a} for deletion: */
        if (debug) { fprintf(stderr, "        merging equation {a} into {b}\n"); }
        (*amodP) = TRUE;
        adel = TRUE;

        /* Modify {b} to have the union of variables and the merged op: */
        int ib;
        for (ib = 0; ib < b->n; ib++) { (void)logsys_unset_arg(b, ib); }
        b->n = un;
        for (ib = 0; ib < un; ib++) { logsys_set_arg(b, ib, vc[ib]); }
        b->op = mop;
        int n_saved = b->n;
        (void)logsys_try_to_simplify_equation(b);
        if (debug) 
          { fprintf(stderr, "        equation {b} modified by merge");
            if (b->n < n_saved) { fprintf(stderr, " (args reduced from %d to %d)", n_saved, b->n); }
            logsys_print_op(stderr, "  b.op = ", b->n, b->op, FALSE, "\n");
          }
        (*bmodP) = TRUE;
        if ((b->n == 0) && (b->op == 1)) { bdel = TRUE; }
      }

    /* Delete what needs to be deleted: */
    if (adel)
      { /* Delete {a}: */
        if (debug) { fprintf(stderr, "        deleting equation {a}\n"); }
        logsys_delete_equation(S, a);
        (*aP) = NULL;
      }
    if (bdel)
      { /* Delete {b}: */
        /* Note that if {a} was deleted it may be the only node now: */
        if (debug) { fprintf(stderr, "        deleting equation {b}\n"); }
        logsys_delete_equation(S, b);
        (*bP) = NULL;
      }
  }

void logsys_unite_eq_args(logsys_eq_t *a, logsys_eq_t *b, int *unP, logsys_va_t *vc[], int ka[], int kb[])
  {
    bool_t debug = FALSE;
    int un = 0; /* Count of variables in union of args. */
    int iab; /* Scans the concatenation of arg lists. */
    for (iab = 0; iab < (a->n + b->n); iab++) 
      { if (iab < a->n)
          { int ia = iab;
            vc[un] = a->arg[ia].va; 
            ka[un] = ia; kb[un] = -1;
            if (debug) 
              { fprintf(stderr, "      joint arg variable vc[%d] from a.arg[%d] = ", un, ia);
                logsys_print_variable_id(stderr, "v", logsys_variable_id(vc[un]), 4, "\n");
              }
            un++;
          }
        else
          { int ib = iab - a->n;
            logsys_va_t *vb = b->arg[ib].va;
            int ic = 0; while ((ic < a->n) && (vc[ic] != vb)) { ic++; }
            if (ic < a->n)
              { kb[ic] = ib;
                if (debug) 
                  { fprintf(stderr, "          joint arg variable vc[%d] also from b.arg[%d] = ", ic, ib);
                    logsys_print_variable_id(stderr, "v", logsys_variable_id(vc[ic]), 4, "\n");
                  }
              }
            else
              { vc[un] = b->arg[ib].va; 
                ka[un] = -1; kb[un] = ib;
                if (debug) 
                  { fprintf(stderr, "          joint arg variable vc[%d] from b.arg[%d] = ", un, ib);
                    logsys_print_variable_id(stderr, "v", logsys_variable_id(vc[un]), 4, "\n");
                  }
                un++;
              }
          }
      }
    (*unP) = un;
  }

void logsys_intersect_eq_args(logsys_eq_t *a, logsys_eq_t *b, int *snP, logsys_va_t *vc[], int ka[], int kb[])
  {
    bool_t debug = FALSE;
    int sn = 0;
    int ia;
    for (ia = 0; ia < a->n; ia++)
      { int ib = 0; 
        while ((ib < b->n) && (a->arg[ia].va != b->arg[ib].va)) { ib++; }
        if (ib < b->n)
          { vc[sn] = a->arg[ia].va;
            ka[sn] = ia; 
            kb[sn] = ib;
            if (debug) 
              { fprintf(stderr, "          shared arg variable vc[%d] a.arg[%d] and b.arg[%d] = ", sn, ia, ib);
                logsys_print_variable_id(stderr, "v", logsys_variable_id(vc[sn]), 4, "\n");
              }
            sn++;
          }
      }
    (*snP) = sn;
  }

vec_typeimpl(logsys_eq_vec_t,logsys_eq_vec,logsys_eq_t*);

logsys_eq_vec_t logsys_collect_equations(logsys_t *S, bool_t tight)
  { 
    if (S->eq == NULL)
      { return logsys_eq_vec_new(0); }
    else 
      { int neq, ieq;
        logsys_eq_t *eq;
        /* Get the number of elements {neq} required in the vector: */
        if (tight)
          { /* Count the equations: */
            neq = 0;
            eq = S->eq; 
            do { neq++; eq = eq->nxeq; } while (eq != S->eq);
          }
        else
          { /* Get the largest equation ID plus one: */
            neq = S->next_eq_id;
          }
        /* Allocate the vector: */
        logsys_eq_vec_t veq = logsys_eq_vec_new(neq);
        /* Collect the equations: */
        if (tight)
          { /* Store in consecutive positions: */
            ieq = 0;
            eq = S->eq; 
            do { veq.e[ieq] = eq; ieq++; eq = eq->nxeq; } while (eq != S->eq);
            assert(ieq == neq);
          }
        else
          { /* Store each equation in the element defined by its ID: */
            for (ieq = 0; ieq < veq.ne; ieq++) { veq.e[ieq] = NULL; }
            eq = S->eq; 
            do { 
              ieq = eq->id;
              assert((ieq >= 0) && (ieq < neq)); /* ID must be valid. */
              assert(veq.e[ieq] == NULL); /* Equations must have unique IDs. */
              veq.e[ieq] = eq;
              eq = eq->nxeq;
            } while (eq != S->eq);
          }
        return veq;
      } 
  }

void logsys_condense_and_split(logsys_t *S)
  {
    bool_t debug = TRUE;
    bool_t verbose = TRUE;

    if (debug) { fprintf(stderr, "  starting {logsys_condense_and_split}\n"); }

    /* Collect all equations in the system into a stack {steq[0..neq-1]}: */
    logsys_eq_vec_t steq = logsys_collect_equations(S, TRUE);
    int nst = steq.ne;
    
    /* Loop over the equations trying to improve the system, until the stack is empty: */
    while (nst > 0)
      { if (verbose || debug) { logsys_print_stats(stderr, S, nst, &steq); }
        /* Get the next pivot equation {fq}: */
        nst--; logsys_eq_t *fq = steq.e[nst];
        assert(fq != NULL);
        if (debug) { logsys_print_equation(stderr, "    pivot {fq} = ", fq, "\n"); }
        /* Try to combine {fq} with all equations that share variables with it: */
        bool_t changed = logsys_try_to_combine_with_others(S->eq, &fq, &nst, &steq);
        if (debug & changed) { fprintf(stderr, "    system changed  nst = %d\n", nst); }
        if (debug & (fq == NULL)) { fprintf(stderr, "    pivot {fq} was deleted\n"); }
      }
    if (debug) { fprintf(stderr, "  finished {logsys_condense_and_split}\n"); }
  }
  
void logsys_print_stats(FILE *wr, logsys_t *S, int nst, logsys_eq_vec_t *steq)
  {
    auto void do_print_stats(bool_t all, char *tag);
      /* Prints stats of {S} if {all} is true, of {steq} if {all} is false. */
    
    void do_print_stats(bool_t all, char *tag)
      { /* Initialize the he counts by arity: */
        int nmax = logsys_n_MAX;
        int count[nmax+1]; /* Count of equations by arity. */
        int totcount = 0; /* Total eq count. */
        int j; 
        for (j = 0; j <= nmax; j++) { count[j] = 0; }

        /* Scan the equations: */
        if (all)
          { /* Scan the system's list: */
            assert(S != NULL);
            if (S->eq != NULL)
              { logsys_eq_t *eq = S->eq;
                do
                  { int n = eq->n;  /* Arity of {eq}. */
                    assert((n >= 0) && (n <= nmax));
                    count[n]++; totcount++; 
                    eq = eq->nxeq;
                  }
                while (eq != S->eq);
              }
          }
        else
          { assert(steq != NULL);
            int ie;
            for (ie = 0; ie < nst; ie++)
              { logsys_eq_t *eq = steq->e[ie];
                if (eq != NULL)
                  { int n = eq->n; /* Arity of {eq}. */
                    assert((n >= 0) && (n <= nmax));
                    count[n]++; totcount++;
                  }
              }
          }
            
        /* Print the counts: */
        fprintf(wr, "  total %d equations in %s - by arity: ", totcount, tag);
        for (j = 0; j <= nmax; j++) 
          { fprintf(stderr, " %d:%d", j,count[j]); }
        fprintf(stderr, "\n");
      }    
    if (S != NULL) { do_print_stats(TRUE, "system"); }
    if (steq != NULL) { do_print_stats(FALSE, "stack"); }
  }
  
void logsys_sort_equations(logsys_t *S)
  {
    bool_t debug = FALSE;
    bool_t print_counts = TRUE;
    
    /* Initialize the tables of equations by arity: */
    int nmax = logsys_n_MAX;
    logsys_eq_t *sq = NULL; /* First equation in sorted list. */
    logsys_eq_t *afirst[nmax+1]; /* Eqn. {afirst[r]} is the first eq in list with arity {r}. */
    int count[nmax+1]; /* Count of equations by arity. */
    int j; 
    for (j = 0; j <= nmax; j++) { afirst[j] = NULL; count[j] = 0; }

    /* Remove them from the system's list and add to the sorted list {sq}: */
    while (S->eq != NULL)
      { /* Remove {S->eq} from the system's list,save in {pq}: */
        logsys_eq_t *pq = S->eq;
        { logsys_eq_t *nxeq = pq->nxeq;
          logsys_eq_t *preq = pq->preq;
          if (nxeq == pq)
            { assert(preq == pq);
              S->eq = NULL;
            }
          else
            { assert(preq != pq);
              assert(preq->nxeq == pq);
              assert(nxeq->preq == pq);
              preq->nxeq = nxeq;
              nxeq->preq = preq;
              S->eq = preq;
            }
          pq->nxeq = pq;
          pq->preq = pq;
        }
        if (debug) { logsys_print_equation(stderr, "    got equation {pq} = ", pq, "\n"); }
        /* Get the arity {n} of {pq}: */
        int n = pq->n;
        assert((n >= 0) && (n <= nmax));
        /* Find the first equation {mq} in the sorted list with arity {m>=n}, or {sq} if none: */
        int m = n;
        while ((m <= nmax) && (afirst[m] == NULL)) { m++; }
        if (debug) { fprintf(stderr, "    bin {m} = %d\n", m); }
        logsys_eq_t *mq = (m <= nmax ? afirst[m] : sq);
        if (mq == NULL)
          { /* Sorted list must be empty: */
            assert(sq == NULL);
          }
        else
          { /* Insert {pq} in sorted list before {mq}: */
            assert(sq != NULL);
            logsys_eq_t *aq = mq->preq;
            assert(aq->nxeq == mq);
            aq->nxeq = pq; pq->preq = aq;
            mq->preq = pq; pq->nxeq = mq;
            if (debug)
              { if (debug) 
                { logsys_print_equation(stderr, "    inserting after {aq} = ", aq, NULL);
                  logsys_print_equation(stderr, " and {mq} = ", mq, "\n");
                }
              }
          }
        /* Save in arity table: */
        afirst[n] = pq; count[n]++;
        /* Update {sq} if needed: */
        if ((sq == NULL) || (sq->n >= n)) { sq = pq; } 
      }
    /* Now put back the sorted list: */
    S->eq = sq;            
    if (debug || print_counts)
      { fprintf(stderr, "  sorted equations by arity");
        for (j = 0; j <= nmax; j++) 
          { fprintf(stderr, " %d:%d", j,count[j]); }
        fprintf(stderr, "\n");
      }
    if (debug) { logsys_print_equation(stderr, "  root equation = ", S->eq, "\n"); }
  }

bool_t logsys_try_to_combine_with_others(logsys_eq_t *eq, logsys_eq_t **fqP, int *nstP, logsys_eq_vec_t *steq)
  {
    bool_t debug = TRUE;
    demand(fqP != NULL, "must give an address");
    logsys_eq_t *fq = (*fqP);
    demand(fq != NULL, "pivot must be non-null");
    demand(eq->sys == fq->sys, "not in the same system");
    
    auto void check_changes(logsys_eq_t *aq, logsys_eq_t *aq_saved, bool_t aq_mod, char *tag);
      /* Checks for changes in {aq}. If changed, tries to split it.
        Updates the stack as appropriate for the changes. */
      
    void check_changes(logsys_eq_t *aq, logsys_eq_t *aq_saved, bool_t aq_mod, char *tag)
      { if (aq_mod)
          { /* Equation {aq} was modified or deleted: */
            if (aq != NULL)
              { /* Apparently {aq} was modified but not deleted. */
                assert(aq_saved == aq);
                /* Try to split it into independent factors: */
                bool_t aq_split = logsys_try_to_split_equation(aq, nstP, steq);
                if (debug)
                  { if (aq_split) 
                      { fprintf(stderr, "      %s was changed and split\n", tag); }
                    else
                      { fprintf(stderr, "      %s was changed but not split\n", tag); }
                  }
                /* Stack unless the split did it: */
                if (! aq_split) { logsys_stack_equation(nstP, steq, aq); }
              }
            else
              { fprintf(stderr, "      %s was deleted\n", tag);
                logsys_remove_equation_from_stack(nstP, steq, aq_saved);
              }
          }
        else
          { assert(aq_saved == aq); }
      }
    
    int nf = fq->n;
    /* Scan the variables of {fq}: */
    int jf;
    for (jf = 0; jf < nf; jf++)
      { logsys_va_t *va = fq->arg[jf].va;
        /* Scan every other equation {gq} that also uses variable {va}: */
        logsys_eq_t *gq = fq->arg[jf].nxeq;
        while (gq != fq)
          { if (debug) { logsys_print_equation(stderr, "      spoke {gq} = ", gq, "\n"); }
            /* Save the pointers for consistency: */
            bool_t gq_mod, fq_mod;
            logsys_eq_t *fq_saved = fq;
            logsys_eq_t *gq_saved = gq;
            /* Try to obtain some simplification from {fq} and {gq}: */
            logsys_try_to_combine_equations(&gq, &gq_mod, &fq, &fq_mod); 
            /* If modification succeeded, try splitting the modified eqs: */
            check_changes(fq, fq_saved, fq_mod, "pivot {fq}");
            check_changes(gq, gq_saved, gq_mod, "spoke {gq}");
            if (fq_mod || gq_mod) 
              { (*fqP) = fq; return TRUE; }
            else
              { /* No change; advance {gq} to next equation of variable {va}: */
                assert(fq_saved == fq); 
                assert(gq_saved == gq);
                if (debug) { fprintf(stderr, "      failed to combine pivot and spoke\n"); }
                int jg = logsys_find_arg(gq, va);
                assert((jg >= 0) && (jg < gq->n));
                gq = gq->arg[jg].nxeq;
              }
          }
      }
    return FALSE;
  }

void logsys_split_op
  ( logsys_op_t op, 
    int n, 
    logsys_op_t *opaP, 
    int *naP, 
    int ka[], 
    logsys_op_t *opbP, 
    int *nbP, 
    int kb[]
  )
  {
    bool_t debug = TRUE;
    demand((n >= 0) && (n <= logsys_n_MAX), "bad num args");
    /* Ready for failure: */
    (*naP) = n; (*opaP) = op; 
    (*nbP) = 0; (*opbP) = logsys_op0_TRU; 
    /* Place argument 0 in the {ka} set:*/
    int na = 0;
    ka[na] = 0; na++;
    if (debug) { fprintf(stderr, "        starting subset with x%d\n", ka[0]); }
    /* Temporarily place arguments {1..n} in the {kb} set: */
    int j;
    int nb = 0;
    for (j = 1; j < n; j++) { kb[nb] = j; nb++; }
    /* Now make sure that elements in {ka} are all independent from those in {kb}: */
    if (n <= 1) { return; }
    bool_t closed = FALSE;
    while (! closed)
      { /* Check whether arguments {ka[0..na-1]} are independent of {kb[0..nb-1]}: */
        int ib = logsys_check_arg_dependency(op, n, na, ka, nb, kb);
        if (ib >= 0)
          { assert(ib < nb);
            if (debug) { fprintf(stderr, "        found dependency with x%d\n", kb[ib]); }
            /* Argument {kb[ib]} depends on {ka[0..na-1]}. */
            /* Remove that element from {kb}: */
            int j = kb[ib];
            kb[ib] = kb[nb-1];
            nb--;
            /* Add it to {ka}: */
            ka[na] = j;
            na++;
          }
        else
          { closed = TRUE; }
      }
    /* Reorder {ka,kb} for nicety: */
    logsys_sort_ints(na, ka);
    logsys_sort_ints(nb, kb);
    /* Split the opcode and return: */
    (*naP) = na; (*opaP) = logsys_scramble_args_in_op(op, n, na, ka); 
    (*nbP) = nb; (*opbP) = logsys_scramble_args_in_op(op, n, nb, kb);
  }

void logsys_sort_ints(int n, int k[])
  { int i, j;
    for (i = 1; i < n; i++)
      { int t = k[i];
        j = i; while ((j > 0) && (k[j-1] > t)) { k[j] = k[j-1]; j--; }
        k[j] = t;
      }
  }

int logsys_check_arg_dependency(logsys_op_t op, int n, int na, int ka[], int nb, int kb[])
  {
    demand((n >= 0) && (n <= logsys_n_MAX), "bad num args");
    demand((na >= 0) && (na <= n), "bad na");
    demand((nb >= 0) && (nb <= n), "bad nb");
    demand(na + nb == n, "not a partition");
    /* If either set is empty, they are independent: */
    if ((na == 0) || (nb == 0)) { return -1; }
    /* Build mask {mska} for all assignments of {ka[0..na-1]} with {kb[0..nb-1]} set to 0: */
    logsys_op_t mska = logsys_build_args_mask(n, na, ka);
    logsys_op_t mskb = logsys_build_args_mask(n, nb, kb);
    /* A {ka}-assignment is an assignment of the {n} arguments where all
      args {kb[0..nb-1]} are zero. A {kb}-assignment is defined
      symmetrically. Thus every assignment of values to the {n} args is
      the OR (or sum) of a {ka}-assignment {X} and a {kb}-assignment
      {Y}; and these two assignments are unique.
      
      A {ka}-assignment {X} is /valid/ if there exist some {kb}-assignment {Y}
      such that {X+Y} is accepted by {op}. The valid {kb}-assignments are defined
      symmetrically.  Note that a valid {ka}- or {kb}-assignment by itself is not
      necessarily accepted by {op}.
      
      The sets {ka,kb} are independent iff {X+Y} is valid for every valid {ka}-assignment {X}
      and every valid {kb}-assignment {Y}. */
    logsys_op_t vala = ZERO64; /* Valid {ka}-assignments, or zero if not found yet. */
    uint32_t Yb = 0; /* If {vala != 0}, this is a valid {kb}-assignment. */
    /* Enumerate all {kb}-assignments {Xb}: */
    uint32_t Nrb = (1 << nb);
    assert(Nrb <= 64);
    uint32_t Xrb;
    for (Xrb = 0; Xrb < Nrb; Xrb++)
      { /* Expand {Xrb} to a {kb}-assignment {Xb}: */
        uint32_t Xb = 0;
        int ib;
        for (ib = 0; ib < nb; ib++) { if ((Xrb & (1 << ib)) != 0) { Xb |= (1 << kb[ib]); } }
        assert(Xb < 64);
        if (mskb && (ONE64 << Xb))
          { /* {Xb} is some {kb}-assignment. */
            /* Get the set {seta} all {ka}-assignments compatible with {Xb}: */
            logsys_op_t seta = (op >> (uint64_t)Xb) & mska;
            /* Check independency condition: */
            if (seta != ZERO64) 
              { /* {Xb} is a valid {kb}-assignment. */
                if (vala == ZERO64)
                  { /* Save the valid {ka}-assignments: */ 
                    vala = seta; Yb = Xb;
                  }
                else if (seta != vala)
                  { /* The argument sets {ka} and {kb} are not independent.
                      Let {Za} be any {ka}-assignment that is accepted by exactly one of {seta,vala}.
                      {Xb} and {Yb} are two valid {kb}-assignments such that exactly 
                      one {Za+Xb} and {Za+Yb} is valid. Any bit where 
                      {Xb} and {Yb} differ is a variable that depends
                      on the {ka} set. */
                    /* Find {ib} so that where {Xb} and {Yb} differ on arg {kb[ib]}: */
                    uint32_t Db = Xb ^ Yb;
                    ib = 0;
                    while ((ib < nb) && (((1 << kb[ib]) & Db) == 0)) { ib++; }
                    assert(ib < nb);
                    return ib;
                  }
              }
          }
      }
    return -1;
  }

logsys_op_t logsys_build_args_mask(int n, int na, int ka[])
  { 
    demand((n >= 0) && (n <= logsys_n_MAX), "bad num args");
    demand((na >= 0) && (na < n), "bad na");
    logsys_op_t op = ONE64;  /* All args set to 0. */
    int ia;
    for (ia = 0; ia < na; ia++)
      { int j = ka[ia];
        demand((j >= 0) && (j < n), "bad arg index");
        /* Expand the set {op} with the assignments where arg{j} is 1: */
        int JX = (1 << j);
        assert(JX < 64);
        op = op | (op << JX);
      }
    return op;
  }

void logsys_normalize(logsys_t *S, int *nu, logsys_va_t ***ua)
  {
    affirm(FALSE, "not implemented");
  }

