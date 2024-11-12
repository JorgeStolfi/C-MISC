/* Last edited on 2018-03-04 23:00:46 by stolfilocal */
/* Checks {logsys.h} system manipulation functions. */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <assert.h>

#include <bool.h>
#include <sign.h>
#include <jsrandom.h>
#include <affirm.h>

#include <logsys.h>

#define ONE64 ((uint64_t)1)
#define ZERO64 ((uint64_t)0)
#define ALL64 ((uint64_t)(-1LL))

void logmrgck_uint64_tests(void);
  /* Tests bit hacks. */

void logmrgck_scramble_tests(void);
  /* Tests {logsys_scramble_args_in_op} and {logsys_unscramble_args_in_op} on various arguments. */
  
void logmrgck_split_tests(void);
  /* Tests {logsys_split_op} on various arguments. */
  
void logmrgck_eq_simplify_tests(void);
  /* Tests equation simplification on some systems. */
  
void logmrgck_combine_and_split_tests(bool_t iterative);
  /* Tests equation merging on some systems. */
  
void logmrgck_scramble_single_test
  ( char *optag, 
    logsys_op_t op_old, 
    int n_old, 
    int n_new, 
    char *ktag, 
    int k[], 
    bool_t verbose,
    logsys_op_t op_new_exp
  );
  /* Tests {logsys__scramble_args_in_op} and
    {logsys_unscramble_args_in_op} on a single set of arguments. Namely,
    scrambles the operation {op_old} with {n_old} arguments to to some
    {op_new} with {n_new} arguments according to the scrambling table
    {k[0..n_new-1]}. Checks that that {op_new} is equal to {op_new_exp}
    and is the result of projecting {op_old} (with OR) on the preserved
    variables and replicating on the new variables. Then unscrambles
    {op_new} and checks whether the result is compatible with
    {op_old}. */

void logmrgck_split_single_test
  ( int n, 
    logsys_op_t opa, 
    int na, 
    int ka[], 
    logsys_op_t opb,  
    int nb,  
    int kb[],
    bool_t verbose
  );
  /* Tests {logsys_split_op} with an equation on {n} arguments that is
    the conjunction of {opa} applied to arguments {ka[0..na-1]} and
    {opb} applied to arguments {kb[0..n-b-1]}. Assumes that {ka,kb} are
    a partition of {0..n-1}. The procedure is aware that
    {logsys_split_op} may find a different factorization, even when
    {na==0} or {nb==0}, but must find SOME non-trivial factorization
    when {na,nb} are both positive. In any case, checks whether
    the factorization found is equivalent to {opa/\opb}. */

void logmrgck_check_scrambled_ops
  ( logsys_op_t op_old,
    int n_old, 
    char *vtag_old,
    logsys_op_t op_new, 
    int n_new, 
    char *vtag_new, 
    int k[], 
    bool_t strict
  );
  /* Checks whether the operator {op_old} on {n_old} arguments called
    "{vtag_old}{0..n_old-1}" is consistent with {op_new} on {n_new}
    arguments called "{vtag_new}{0..n_new-1}" that are related to those
    of {n_old} by the table {k[0..n_new-1]}. If {strict} is true, checks
    that {op_old} is independent of any old argument that is not listed
    in {k[0..n_new-1]}. */

void logmrgck_check_system_equiv(logsys_t *S0, logsys_t *S1);
  /* Checks whether {S0} and {S1} have the same set of solutions.
    Systems {S0,S1} must have the same number of variables with the same sets 
    of IDs. The procedure enumerates all possible assignments
    of boolean values to those variables, such that
    variables with same {id} are assigned the same value,
    and complains if one system accepts a solution
    that the other one rejects. */

void logmrgck_create_empty_systems(int nv, int ne, logsys_t *S[], logsys_va_t **va[], logsys_eq_t **eq[]);
  /* The first level of the arrays {S,va,eq} must be allocated by caller
    with size 2. For each {i} in {0..1}, the procedure stores into
    {S[i]} a new empty system, into {va[i]} a newly allocated vector
    with {nv} entries, and into {eq[i]} a newly allocated vector with
    {ne} entries. It then adds {nv} variables to system {S[i]} and
    stores their addresses into {va[i][0..nv-1]}. Does not add any
    equation to the systems.
    
    Note that the vector {va[i][0..nv-1]} is NOT suitable for
    {logsys_print_var_assignments} and the like, since the ID of a
    variable {va[i][k]} is usually NOT {k}. */

int main(int argc, char **argv)
  {
    /* logmrgck_uint64_tests();  */   
    logmrgck_scramble_tests();
    logmrgck_split_tests();
    logmrgck_eq_simplify_tests();
    logmrgck_combine_and_split_tests(FALSE);
    logmrgck_combine_and_split_tests(TRUE);
    return 0;
  }
  
#define uint64_hex_fmt "lX"
  /* Hex format suitable for {uint64_t} */

void logmrgck_uint64_tests(void)
  {
    int I64 = (1 << 6);
    uint64_t op6 = (logsys_op_t)0XF731F731F731F731LLU;
    
    uint64_t test;

    int i;
    for (i = -2; i < 65; i++)
      { 
        test = (ONE64 << i);
        fprintf(stderr, "test = %016" uint64_hex_fmt " = (ONE64 << i), i = %d\n", test, i);

        test = (ONE64 << (uint64_t)i);
        fprintf(stderr, "test = %016" uint64_hex_fmt " = (ONE64 << (uint64_t)i), i = %d\n", test, i);

        test = (ONE64 << (uint64_t)i);
        fprintf(stderr, "test = %016" uint64_hex_fmt " = (ONE64 << (uint64_t)i), i = %d\n", test, i);

        test = (ONE64 << (uint64_t)i) + ALL64;
        fprintf(stderr, "test = %016" uint64_hex_fmt " = (ONE64 << (uint64_t)i) + ALL64, i = %d\n", test, i);

        test = (op6 << i);
        fprintf(stderr, "test = %016" uint64_hex_fmt " = (op6 << i), i = %d\n", test, i);

        test = (op6 >> i);
        fprintf(stderr, "test = %016" uint64_hex_fmt " = (op6 >> i), i = %d\n", test, i);

        test = (op6 >> (uint64_t)i);
        fprintf(stderr, "test = %016" uint64_hex_fmt " = (op6 >> (uint64_t)i), i = %d\n", test, i);
        
        fprintf(stderr, "\n");
      }

    test = (ONE64 << 33) - 1;
    fprintf(stderr, "test = %016" uint64_hex_fmt " = (ONE64 << 33) - 1\n", test);

    test = (ONE64 << (uint64_t)33) - 1;
    fprintf(stderr, "test = %016" uint64_hex_fmt " = (ONE64 << (uint64_t)33) - 1\n", test);

    test = (ONE64 << 33) - ONE64;
    fprintf(stderr, "test = %016" uint64_hex_fmt " = (ONE64 << 33) - ONE64\n", test);

    test = (ONE64 << (uint64_t)33) - ONE64;
    fprintf(stderr, "test = %016" uint64_hex_fmt " = (ONE64 << (uint64_t)33) - ONE64\n", test);

    test = (logsys_op_t)0X0029400040000029LLU;
    fprintf(stderr, "test = %016" uint64_hex_fmt " = (logsys_op_t)0X0029400040000029LLU\n", test);

    test = (op6 >> I64);
    fprintf(stderr, "test = %016" uint64_hex_fmt " = (op6 >> ttn)\n", test);

    test = (op6 >> (uint64_t)I64);
    fprintf(stderr, "test = %016" uint64_hex_fmt " = (op6 >> (uint64_t)ttn)\n", test);
  }
  
void logmrgck_scramble_tests(void)
  { 
    fprintf(stderr, "=== testing op argument scrambling and unscrambling ================\n");
    bool_t verbose = FALSE;
    
    /* Pick an operation on 2 variables that depends on both: */
    logsys_op_t op2 = ((logsys_op2_t)0x2u);  /* Requires {v0==1, v1==0} */
    assert(logsys_op_is_valid(op2, 2));
    assert(logsys_op_depends(op2, 2, 0));
    assert(logsys_op_depends(op2, 2, 1));
    assert(! logsys_op_is_functional(op2, 2, 0));
    assert(! logsys_op_is_functional(op2, 2, 1));
    assert(! logsys_op_is_solvable(op2, 2, 0));
    assert(! logsys_op_is_solvable(op2, 2, 1));
    
    /* Trivial scramble: */
    int k2a[2] = { 0, 1 };
    logmrgck_scramble_single_test("op2", op2, 2, 2, "k2a", k2a, verbose, ((logsys_op2_t)0x2u));
                                                     
    /* Arg swap: */                                  
    int k2b[2] = { 1, 0 };                             
    logmrgck_scramble_single_test("op2", op2, 2, 2, "k2b", k2b, verbose, ((logsys_op2_t)0x4u));
                                                     
    /* Arg 0 replacement: */                         
    int k2c[2] = { -1, 1 };                            
    logmrgck_scramble_single_test("op2", op2, 2, 2, "k2c", k2c, verbose, ((logsys_op2_t)0x3u));
                                                     
    /* Arg 1 replacement: */                         
    int k2d[2] = { 0, -1 };                            
    logmrgck_scramble_single_test("op2", op2, 2, 2, "k2d", k2d, verbose, ((logsys_op2_t)0xAu));
    
    /* Swap and arg 0 replacement: */
    int k2e[2] = { -1, 0 };
    logmrgck_scramble_single_test("op2", op2, 2, 2, "k2e", k2e, verbose, ((logsys_op2_t)0xCu));
    
    /* Swap and arg 1 replacement: */
    int k2f[2] = { 1, -1 };
    logmrgck_scramble_single_test("op2", op2, 2, 2, "k2f", k2f, verbose, ((logsys_op2_t)0x5u));
                                                     
    /* Insertion in the middle: */                   
    int k3g[3] = { 0, -1, 1 };                         
    logmrgck_scramble_single_test("op2", op2, 2, 3, "k3g", k3g, verbose, ((logsys_op3_t)0x0Au));

    /* Pick an operation on 3 variables that depends on all of them: */
    logsys_op_t op3 = logsys_op3_AND;
    assert(logsys_op_is_valid(op3, 3));
    assert(logsys_op_depends(op3, 3, 0));
    assert(logsys_op_depends(op3, 3, 1));
    assert(logsys_op_depends(op3, 3, 2));
    assert(logsys_op_is_functional(op3, 3, 0));
    assert(! logsys_op_is_functional(op3, 3, 1));
    assert(! logsys_op_is_functional(op3, 3, 2));
    assert(logsys_op_is_solvable(op3, 3, 0));
    assert(! logsys_op_is_solvable(op3, 3, 1));
    assert(! logsys_op_is_solvable(op3, 3, 2));

    /* Trivial scramble: */
    int k3a[3] = { 0, 1, 2 };
    logmrgck_scramble_single_test("op3", op3, 3, 3, "k3a", k3a, verbose, logsys_op3_AND);
                                                     
    /* Cyclic perm: */
    int k3b[3] = { 2, 0, 1 };
    logmrgck_scramble_single_test("op3", op3, 3, 3, "k3b", k3b, verbose, ((logsys_op3_t)0x93u));

    /* Replace arg 1: */
    int k3c[3] = { 0, -1, 2 };
    logmrgck_scramble_single_test("op3", op3, 3, 3, "k3c", k3c, verbose, ((logsys_op3_t)0xF5u));

    /* Insert arg 1: */
    int k4d[4] = { 0, -1, 1, 2 };
    logmrgck_scramble_single_test("op3", op3, 3, 4, "k4d", k4d, verbose, ((logsys_op4_t)0xA555u));

    /* Delete arg 1: */
    int k2h[2] = { 0, 2 };
    logmrgck_scramble_single_test("op3", op3, 3, 2, "k2h", k2h, verbose, ((logsys_op2_t)0xDu));

    /* Pick an operation on 6 variables: */
    logsys_op_t op6 = (logsys_op_t)0X0029400040000029LLU;
    assert(logsys_op_is_valid(op6, 6));

    /* Trivial scramble: */
    int k6a[6] = { 0, 1, 2, 3, 4, 5 };
    logmrgck_scramble_single_test("op6", op6, 6, 6, "k6a", k6a, verbose, op6);

    /* Extract {x0,x5}:: */
    int k2i[2] = { 0, 5 };
    logmrgck_scramble_single_test("op6", op6, 6, 2, "k2i", k2i, verbose, ((logsys_op2_t)0xF));
                                                     
  }
    
void logmrgck_scramble_single_test
  ( char *optag, 
    logsys_op_t op_old, 
    int n_old, 
    int n_new, 
    char *ktag, 
    int k[], 
    bool_t verbose,
    logsys_op_t op_new_exp
  )
  {     
    /* Apply scramble: */
    fprintf(stderr, "--- testing scramble/unscramble of %s with n_old = %d  n_new = %d  k = %s---\n", optag, n_old, n_new, ktag);
    if (verbose) { logsys_print_op(stderr, "op_old = \n", n_old, op_old, TRUE, "\n"); }

    /* Apply scramble: */
    logsys_op_t op_new = logsys_scramble_args_in_op(op_old, n_old, n_new, k);
    if (verbose) { logsys_print_op(stderr, "op_new = \n", n_new, op_new, TRUE, "\n"); }
    
    /* Check whether the scramble succeeded: */
    logmrgck_check_scrambled_ops(op_old, n_old, "x", op_new, n_new, "y", k, FALSE);
    assert(op_new == op_new_exp);
    
    /* Build the inverse argument map {unk[0..n_old-1]}: */
    int j_old;
    int unk[n_old];
    for (j_old = 0; j_old < n_old; j_old++)
      { int j_new = n_new - 1;
        while ((j_new >= 0) && (k[j_new] != j_old)) { j_new--; }
        unk[j_old] = j_new;
      }

    /* Check {logsys_unscramble_args_in_op}: */
    logsys_op_t op_uns = logsys_unscramble_args_in_op(op_new, n_old, n_new, k);
    if (verbose) { logsys_print_op(stderr, "op_uns = \n", n_old, op_uns, TRUE, "\n"); }
    
    /* Check whether the unscramble succeeded: */
    logmrgck_check_scrambled_ops(op_new, n_new, "y", op_uns, n_old, "x", unk, TRUE);
    fprintf(stderr, "------------------------------------------------------------\n");
  }
  
void logmrgck_split_tests(void)
  {
    fprintf(stderr, "=== testing equation simplification ================\n");
    int ntries = 100;
    int try;
    srandom(71400417);
    for (try = 0; try < ntries; try++)
      { bool_t verbose = (try < 20);
        /* Choose the total arg count {n}: */
        int n = int32_abrandom(0,logsys_n_MAX);
        /* Generate a random bipartition {ka[0..na],kb[0..nb]} of {0..n-1}: */
        /* Note that either part may be empty. */
        int ka[n], kb[n];
        int na = 0, nb = 0;
        int j;
        for(j = 0; j < n; j++) 
          { if (drandom() < 0.5) { ka[na] = j; na++; } else { kb[nb] = j; nb++; } }
        /* Generate a random equation for each side: */
        logsys_op_t opa = uint64_random() & ((ONE64 << (1 << na)) + ALL64);
        logsys_op_t opb = uint64_random() & ((ONE64 << (1 << nb)) + ALL64);
        /* Test with those ops: */
        logmrgck_split_single_test(n, opa, na, ka, opb, nb, kb, verbose);
      }
  }
    
void logmrgck_split_single_test
  ( int n, 
    logsys_op_t iopa, 
    int na, 
    int ka[], 
    logsys_op_t iopb,  
    int nb,  
    int kb[],
    bool_t verbose
  )
  {
    fprintf(stderr, "--- testing op split with ");
    logsys_print_op(stderr, NULL, na, iopa, FALSE, NULL);
    logsys_print_arg_set(stderr, "(", "x", "u", na, ka, NULL, ",", ")");
    fprintf(stderr, " and ");
    logsys_print_op(stderr, NULL, nb, iopb, FALSE, NULL);
    logsys_print_arg_set(stderr, "(", "x", "v", nb, kb, NULL, ",", ")");
    fprintf(stderr, " ---\n");
    assert((n >= 0) && (n <= logsys_n_MAX));
    assert((na >= 0) && (na <= n));
    assert((nb >= 0) && (nb <= n));
    assert(na + nb == n);
    /* Unscramble the two opcodes and combine them into a single {n}-ary opcode {uop}: */
    logsys_op_t uopa = logsys_unscramble_args_in_op(iopa, n, na, ka);
    logsys_op_t uopb = logsys_unscramble_args_in_op(iopb, n, nb, kb);
    logsys_op_t uop = uopa & uopb;
    /* Now try to split that equation: */
    logsys_op_t sopa, sopb;
    int ska[n], skb[n];
    int sna, snb;
    logsys_split_op(uop, n, &sopa, &sna, ska, &sopb, &snb, skb);
    /* Merge the two factors found: */
    logsys_op_t zopa = logsys_unscramble_args_in_op(sopa, n, sna, ska);
    logsys_op_t zopb = logsys_unscramble_args_in_op(sopb, n, snb, skb);
    logsys_op_t zop = zopa & zopb;
    if (verbose)
      { /* Print factorization: */
        fprintf(stderr, "factorization found:\n");
        fprintf(stderr, "  sopa = ");
        logsys_print_op(stderr, NULL, sna, sopa, FALSE, NULL);
        logsys_print_arg_set(stderr, "(", "x", "u", sna, ska, NULL, ",", ")\n");
        fprintf(stderr, "  sopb = ");
        logsys_print_op(stderr, NULL, snb, sopb, FALSE, NULL);
        logsys_print_arg_set(stderr, "(", "x", "v", snb, skb, NULL, ",", ")\n");
        fprintf(stderr, "operations with unified argument list:\n");
        logsys_op_t prtop[6] = { uopa, uopb, uop, zopa, zopb, zop };
        logsys_print_ops(stderr, "original factors, product, found factors, product:\n", n, 6, prtop, NULL);
      }
    /* Check factor sizes: */
    demand((sna >= 0) && (sna <= n), "invalid sna");
    demand((snb >= 0) && (snb <= n), "invalid snb");
    if (n > 0) { demand(sna > 0, "first factor is empty"); }
    /* If the input was composite, check whether split succeeded: */
    if ((na > 0) && (nb > 0)) { demand((sna > 0) && (snb > 0), "failed to find a split"); }
    /* Check whether {ka[0..na-1]} and {kb[0..nb-1]} are valid and a partition of {0..n-1}: */
    demand(sna + snb == n, "invalid factor sizes");
    bool_t used[n];
    int j;
    for(j = 0; j < n; j++) { used[j] = FALSE; }
    int ia, ib;
    for (ia = 0; ia < sna; ia++) 
      { demand((ska[ia] >= 0) && (ska[ia] < n), "bad {ska[ia]}");
        demand(! used[ska[ia]], "repeated arg index in {ska}");
        used[ska[ia]] = TRUE;
      }
    for (ib = 0; ib < snb; ib++) 
      { demand((skb[ib] >= 0) && (skb[ib] < n), "bad {skb[ib]}");
        demand(! used[skb[ib]], "repeated arg index in {skb}");
        used[skb[ib]] = TRUE;
      }
    /* Now check whether {sopa/\sopb} is equivalent to {iopa/\iopb}: */
    uint32_t N = (1 << n);
    uint32_t Xu;
    for (Xu = 0; Xu < N; Xu++)
      { 
        /* Get the decision of {uop} on {Xu}: */
        bool_t uop_OK = ((uop & (ONE64 << (uint64_t)Xu)) != ZERO64);
        /* Split {Xu} into {Xua,Xub} according to {ka,kb}: */
        uint32_t Xua = 0, Xub = 0;
        for (ia = na-1; ia >= 0; ia--) { Xua = Xua << 1; if ((Xu & (1 << ka[ia])) != 0) { Xua |=  1; } } 
        for (ib = nb-1; ib >= 0; ib--) { Xub = Xub << 1; if ((Xu & (1 << kb[ib])) != 0) { Xub |=  1; } } 
        /* Get the decision of {iopa,iopb} on {Xua,Xub}, must agree  with {uop}: */
        bool_t iopa_OK = ((iopa & (ONE64 << (uint64_t)Xua)) != ZERO64);
        bool_t iopb_OK = ((iopb & (ONE64 << (uint64_t)Xub)) != ZERO64);
        assert(uop_OK == (iopa_OK && iopb_OK));
        
        /* Get the decision of {zop} on {Xu}: */
        bool_t zop_OK = ((zop & (ONE64 << (uint64_t)Xu)) != ZERO64);
        /* Also split {Xu} into {Xsa,Xsb} according to {ska,skb}: */
        uint32_t Xsa = 0, Xsb = 0;
        for (ia = sna-1; ia >= 0; ia--) { Xsa = Xsa << 1; if ((Xu & (1 << ska[ia])) != 0) { Xsa |=  1; } } 
        for (ib = snb-1; ib >= 0; ib--) { Xsb = Xsb << 1; if ((Xu & (1 << skb[ib])) != 0) { Xsb |=  1; } } 
        /* Get the decision of {sopa,sopb} on {Xsa,Xsb}, must also agree  with {zop} and {uop}: */
        bool_t sopa_OK = ((sopa & (ONE64 << (uint64_t)Xsa)) != ZERO64);
        bool_t sopb_OK = ((sopb & (ONE64 << (uint64_t)Xsb)) != ZERO64);
        assert(zop_OK == (sopa_OK && sopb_OK));
        
        /* Check whether found fctorization agrees with original: */
        if(uop_OK != zop_OK)
          { fprintf(stderr, "** bug in factored eq\n");
            fprintf(stderr, "  Xu = %X uop(Xu) = %c\n", Xu, "FT"[uop_OK]);
            fprintf(stderr, "  Xua = %X iopa(Xua) = %c\n", Xua, "FT"[iopa_OK]);
            fprintf(stderr, "  Xub = %X iopb(Xub) = %c\n", Xub, "FT"[iopb_OK]);
            fprintf(stderr, "  Xu = %X zop(Xu) = %c\n", Xu, "FT"[zop_OK]);
            fprintf(stderr, "  Xsa = %X sopa(Xsa) = %c\n", Xsa, "FT"[sopa_OK]);
            fprintf(stderr, "  Xsb = %X sopb(Xsb) = %c\n", Xsb, "FT"[sopb_OK]);
            assert(FALSE);
          }
      }
    fprintf(stderr, "------------------------------------------------------------\n");
  }

void logmrgck_eq_simplify_tests(void)   
  {
    fprintf(stderr, "=== testing equation simplification ================\n");
    /* Pick an operation on 3 variables that depends on all of them: */
    logsys_op_t op3 = logsys_op3_AND;
    logsys_print_op(stderr, "op3 = ", 3, op3, TRUE, "\n");
    /* Create an operation on 5 variables that does not depend on variables {1,3}: */
    int k5[5] = {0, -1, 1, -1, 2};
    logsys_op_t op5 = logsys_scramble_args_in_op(op3, 3, 5, k5);
    logsys_print_op(stderr, "op5 = ", 5, op5, TRUE, "\n");
    logmrgck_check_scrambled_ops(op3, 3, "x", op5, 5, "y", k5, FALSE);
    
    /* Check {logsys_try_to_simplify_equation}: */
    int my_nv = 5;
    int my_ne = 1;

    /* Create two systems with {my_nv} variables: */
    logsys_va_t **my_va[2]; 
    logsys_eq_t **my_eq[2];
    logsys_t *S[2];
    logmrgck_create_empty_systems(my_nv, my_ne, S, my_va, my_eq);

    /* Add to both systems a pair of equivalent equations: */
    assert(my_nv == 5);
    assert(my_ne == 1);
    
    /* Add to System {S[0]} an AND-gate on variables {v0,v2,v4}: */
    logsys_va_t *arg00[3] = { my_va[0][0], my_va[0][2], my_va[0][4] };
    my_eq[0][0] = logsys_add_equation(S[0], op3, 3, arg00, NULL);
    assert(logsys_num_args(my_eq[0][0]) == 3);

    /* Add to System {S[1]} the equivalent 5-ary gate on variables {v0,v1,v2,v3,v4}: */
    logsys_va_t *arg10[5] = { my_va[1][0], my_va[1][1], my_va[1][2], my_va[1][3], my_va[1][4] };
    my_eq[1][0] = logsys_add_equation(S[1], op5, 5, arg10, NULL);
    assert(logsys_num_args(my_eq[1][0]) == 5);
    
    /* Check whether the systems are really equivalent: */
    logmrgck_check_system_equiv(S[0], S[1]);

    /* Simplify the 5-element equation of system 1: */
    fprintf(stderr, "--- simplifying equation 0 of system 1 ---\n");
    logsys_print_equation(stderr, "before = ", my_eq[1][0], "\n");
    bool_t changed = logsys_try_to_simplify_equation(my_eq[1][0]);
    if (changed) { fprintf(stderr, "success claimed\n"); } 
    logsys_print_equation(stderr, "after  = ", my_eq[1][0], "\n");
    assert(logsys_num_args(my_eq[1][0]) == 3);

    /* Check whether the systems are still equivalent: */
    logmrgck_check_system_equiv(S[0], S[1]);
  }

void logmrgck_combine_and_split_tests(bool_t iterative)   
  {
    fprintf(stderr, "=== testing equation combine/split (iterative = %c) ================\n","01"[iterative]);
    /* Check {logsys_try_to_simplify_equation}: */
    int my_nv = 9;
    int my_ne = 4;

    /* Create two systems with {my_nv} variables: */
    logsys_va_t **my_va[2]; 
    logsys_eq_t **my_eq[2];
    logsys_t *S[2];
    logmrgck_create_empty_systems(my_nv, my_ne, S, my_va, my_eq);

    /* Add to each system two equations that can be merged and one that can be split: */
    assert(my_nv == 9);
    assert(my_ne == 4);
    int isys;
    int icomba, icombb, isplit; /* Equation indices for non-iterative tests. */
    for (isys = 0; isys < 2; isys++)
      { 
        logsys_va_t **v = my_va[isys];
        logsys_eq_t **e = my_eq[isys];
        
        int ie = 0;
        /* Add an AND-gate {v2 == v1/\v0}: */
        logsys_va_t *arg0[3] = { v[2], v[1], v[0] };
        e[ie] = logsys_add_equation(S[isys], logsys_op3_AND, 3, arg0, NULL);
        icomba = ie;
        ie++;

        /* Add an OR-gate {v3 == v1\/v0}: */
        logsys_va_t *arg1[3] = { v[3], v[1], v[0] };
        e[ie] = logsys_add_equation(S[isys], logsys_op3_ORR, 3, arg1, NULL);
        icombb = ie;
        ie++;

        /* Add a GEQ-gate {v2 \geq v3}: */
        logsys_va_t *arg2[2] = { v[2], v[3] };
        e[ie] = logsys_add_equation(S[isys], logsys_op2_G01, 2, arg2, NULL);
        ie++;
        
        /* Add a gate {(v5 != v6) /\ (v7 \geq v8)}: */
        logsys_va_t *arg3[4] = { v[5], v[6], v[7], v[8] };
        e[ie] = logsys_add_equation(S[isys], (logsys_op_t)0X6066, 4, arg3, NULL);
        isplit = ie;
        ie++;
        
        assert(ie == my_ne);
      }
    logsys_print_system(stderr, "--- system {S1} before ---", S[1], "--------------------------");
    
    /* Check whether the systems are really equivalent: */
    logmrgck_check_system_equiv(S[0], S[1]);

    if (! iterative)
      { /* Apply single combine step to equations {icomba,icombb} of {S[1]}: */
        fprintf(stderr, "--- trying to merge equations %d and %d of {S[1]} ---\n", icomba, icombb);
        isys = 1;
        logsys_eq_t *a = my_eq[isys][icomba];
        logsys_eq_t *b = my_eq[isys][icombb];
        bool_t amod, bmod;
        logsys_try_to_combine_equations(&a, &amod, &b, &bmod);
        fprintf(stderr, "  modified equations: eq[%d] = %c eq[%d] = %c\n", icomba, "01"[amod], icombb, "01"[bmod]);
        if (a == NULL) { fprintf(stderr, "  eq[%d] deleted\n", icomba); my_eq[isys][icomba] = NULL; }
        if (b == NULL) { fprintf(stderr, "  eq[%d] deleted\n", icombb); my_eq[isys][icombb] = NULL; }
        
        /* Apply single split step to equation {isplit} of {S[1]}: */
        fprintf(stderr, "--- trying to split equation %d of {S[1]} ---\n", isplit);
        isys = 1;
        logsys_eq_t *c = my_eq[isys][isplit];
        bool_t split = logsys_try_to_split_equation(c, NULL, NULL);
        fprintf(stderr, "  split equation: eq[%d] = %c\n", isplit, "01"[split]);
        if (c == NULL) { fprintf(stderr, "  eq[%d] deleted\n", isplit); my_eq[isys][isplit] = NULL; }
      }
    else
      { /* Apply the iterative merge to {S[1]}: */
        logsys_condense_and_split(S[1]);
      }
    logsys_print_system(stderr, "--- system {S1} after ---", S[1], "--------------------------");

    /* Check whether the systems are still equivalent: */
    logmrgck_check_system_equiv(S[0], S[1]);
  }

void logmrgck_create_empty_systems(int nv, int ne, logsys_t *S[], logsys_va_t **va[], logsys_eq_t **eq[])
  {
    int isys, iv, ie;
    for (isys = 0; isys < 2; isys++)
      { va[isys] = malloc(nv*sizeof(logsys_va_t *));
        eq[isys] = malloc(ne*sizeof(logsys_eq_t *));
        S[isys] = logsys_new_system();
        for (iv = 0; iv < nv; iv++) { va[isys][iv] = logsys_add_variable(S[isys]); }
        for (ie = 0; ie < ne; ie++) { eq[isys][ie] = NULL; }
      }
  }    

void logmrgck_check_system_equiv(logsys_t *S0, logsys_t *S1)
  {
    fprintf(stderr, "  checking system equivalence...\n");
    int nv0, nv1;
    logsys_va_t **va0, **va1;
    logsys_get_variables(S0, &nv0, &va0);
    logsys_get_variables(S1, &nv1, &va1);

    /* Check that the systems have the same number of variables with same set of valid IDs: */
    assert(nv0 == nv1);
    int nv = nv0;
    int iv;
    for (iv = 0; iv < nv; iv++) 
      { assert((va0[iv] == NULL) == (va1[iv] == NULL));
        if (va0[iv] != NULL)
          { assert(logsys_variable_id(va0[iv]) == iv); 
            assert(logsys_variable_id(va1[iv]) == iv);
          }
      }
    
    /* Enumerate all assignments of values to those variables and check whether systems agree: */
    bool_t x[nv];
    for (iv = 0; iv < nv; iv++) { x[iv] = FALSE; }
    do
      { bool_t res0 = logsys_is_solution(S0, nv, va0, x);
        bool_t res1 = logsys_is_solution(S1, nv, va1, x);
        if (res0 != res1) 
          { fprintf(stderr, "** systems are not equivalent: S0 = %c  S1 = %c", "01"[res0], "01"[res1]);
            logsys_print_var_assignments(stderr, "  x = ", nv, va0, x, "\n");
          }
        /* Get next assignment: */
        iv = nv - 1;
        while((iv >= 0) && ((va0[iv] == NULL) || x[iv])) { x[iv] = FALSE; iv--; }
        if (iv >= 0) { assert(va0[iv] != NULL); x[iv] = TRUE; }
      }
    while (iv >= 0);
    fprintf(stderr, "  systems are equivalent\n");
  }

void logmrgck_check_scrambled_ops
  ( logsys_op_t op_old,
    int n_old, 
    char *vtag_old,
    logsys_op_t op_new, 
    int n_new, 
    char *vtag_new, 
    int k[], 
    bool_t strict
  )
  {
    bool_t debug = FALSE;
    int NY = (1 << n_new);
    int NX = (1 << n_old);
    /* Enumerate all assignments {Y} to the new variables: */
    int X, Y;
    for (Y = 0; Y < NY; Y++)
      { bool_t y[n_new];
        logsys_get_bits(Y, n_new, y);
        bool_t new_OK = ((op_new & (ONE64 << (uint64_t)Y)) != 0);
        if (debug) 
          { logsys_print_assignments(stderr, "  checking new assignment ", n_new, y, ""); 
            fprintf(stderr, " value = %c\n", "01"[new_OK]);
          }
        /* Check all assignments {X} to old variables that are compatible with {Y}: */
        int num_old_OK = 0; /* Count of old compatible assignments allowed by {op_old}. */
        int num_old_KO = 0; /* Count of old compatible old assignments forbidden by {op_old}. */
        for (X = 0; X < NX; X++)
          { bool_t x[n_old];
            logsys_get_bits(X, n_old, x);
            bool_t old_OK = ((op_old & (ONE64 << (uint64_t)X)) != 0);
            bool_t compatible = TRUE;
            int jy;
            for (jy = 0; compatible && (jy < n_new); jy++)
              { int jx = k[jy];
                if (jx >= 0)
                  { assert(jx < n_old);
                    if (x[jx] != y[jy]) { compatible = FALSE; }
                  }
              }
            if (compatible)
              { if (debug) 
                  { logsys_print_assignments(stderr, "    compatible old assignment ", n_old, x, ""); 
                    fprintf(stderr, " value = %c\n", "01"[old_OK]);
                  }
                if (old_OK) { num_old_OK++; } else { num_old_KO++; }
              }
          }
        bool_t some_old_OK = (num_old_OK > 0); /* TRUE if {old_op} accepted some compar assignment. */
        bool_t bug = FALSE; /* Some problem at the new assignment {Y}. */
        if (strict)
          { /* Check whether the old opcode depended on the lost variables: */
            if ((num_old_OK > 0) && (num_old_KO > 0))
              { fprintf(stderr, "** old opcode depended on omitted variables\n");
                bug = TRUE;
              }
          }
        /* Check whether the new opcode gives the expected result: */
        if (some_old_OK != new_OK)
          { fprintf(stderr, "** scramble/unscramble error\n");
            bug = TRUE;
          }
          
        if (bug)
          { /* Print the state: */
            fprintf(stderr, "  old variables = ");
            logsys_print_arg_set(stderr, "(", vtag_old, vtag_new, n_old, NULL, NULL, ",", ")");
            logsys_print_op(stderr, "op_old = \n", n_old, op_old, TRUE, "\n");
            fprintf(stderr, "  new variables = ");
            logsys_print_arg_set(stderr, "(", vtag_old, vtag_new, n_new, k, NULL, ",", ")");
            logsys_print_op(stderr, "op_new = \n", n_new, op_new, TRUE, "\n");
            fprintf(stderr, "  value assignment of new args: ");
            logsys_print_arg_set(stderr, NULL, vtag_old, vtag_new, n_new, k, y, ", ", NULL);
            fprintf(stderr, "  op_old = %c (%d allowed %d forbidden)", "01"[some_old_OK], num_old_OK, num_old_KO);
            fprintf(stderr, "  op_new = %c", "01"[new_OK]);
            fprintf(stderr, "\n");
            assert(FALSE);
          }
      }
  }

