#ifndef logsys_H
#define logsys_H

/* Last edited on 2012-12-20 17:20:19 by stolfilocal */
/* General system solving. */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>

#include <bool.h>
#include <vec.h>

typedef struct logsys_t logsys_t;
  /* A system of boolean equations. */

typedef struct logsys_va_t logsys_va_t;
  /* A boolean variable belonging to some system. */

typedef struct logsys_eq_t logsys_eq_t;
  /* A logical equation depending on up to three variables. */

vec_typedef(logsys_eq_vec_t,logsys_eq_vec,logsys_eq_t*);
  /* An extensible vector of {logsys_eq_t *}. */

/* BOOLEAN EQUATION OPCODES

    A Boolean equation (predicate) on {n} variables is represented by
    its truth table encoded as an unsigned integer {op} with at least
    {2^n} bits whose '1' bits identify the variable value combinations
    that satisfy the equation.
    
    More precisely, each tuple of {n} boolean values {x[0..n-1]} is
    mapped to the integer {X = SUM{ x[i]*2^i : i IN 0..n-1 }}. Then
    bit number {X} of {op} (whose numeric value is {2^X} is 1 if and
    only if the equation {op(va[0],..va[n-1])} is satisfied when
    {va[i]==x[i]} for {i} in {0..n-1}. */

#define logsys_n_MAX 6
  /* Max number of argument variables in an equation. */

#define logsys_NX_MAX (1<<logsys_n_MAX) 
  /* Max number of assignments for the argument variables in an equation. */

typedef uint64_t logsys_op_t;
  /* Identifies a boolean equation on up to 6 variables. */

typedef logsys_op_t logsys_op0_t;
  /* Identifies a boolean equation on zero variables. Only the
    lowest-order bit is used. Actually there are only two opcodes of
    this type: 0 (never satisfied) and 1 (always satisfied). */

typedef logsys_op_t logsys_op1_t;
  /* Identifies a boolean equation on a single variable {v0}. Only the
    lowest 2 bits are used. There are four such equations: 0 (never
    satisfied, for any value of {v0}), 1 (requires {v0 == 0}), 2
    (requires {v0 == 1}) and 3 (always satisfied, for any value of
    {v0}). */

typedef logsys_op_t logsys_op2_t;
  /* Identifies a boolean equation on two variables {v0,v1}. Only the
    lowest 4 bits are used. There are 16 such equations. */

typedef logsys_op_t logsys_op3_t;
  /* Identifies a boolean equation on three variables {v0,v1,v2}. Only the
    lowest 8 bits are used. There are 256 such equations. */

typedef logsys_op_t logsys_op4_t;
  /* Identifies a boolean equation on four variables {v0,v1,v2,v3}.
    Only the lowest 16 bits are used. There are 65,536 such equations. */
    
typedef logsys_op_t logsys_op5_t;
  /* Identifies a boolean equation on five variables {v0,v1,v2,v3,v4}.
    Only the lowest 32 bits are used. There are 4,294,967,296 such equations. */
    
typedef uint64_t logsys_op6_t;
  /* Identifies a boolean equation on six variables {v0,v1,v2,v3,v4,v5}.
    Only the lowest 64 bits are used. There are 18,446,744,073,709,551,616 
    such equations. */

/* The zero-ary gates: */
#define logsys_op0_TRU ((logsys_op1_t)0x1u) /* Always satisfied. */
#define logsys_op0_FAL ((logsys_op1_t)0x0u) /* Always failed. */

/* Some popular unary gates: */
#define logsys_op1_ZR0 ((logsys_op1_t)0x1u) /* {v0 == 0} (hard-wires 0 to {v0}). */
#define logsys_op1_UN0 ((logsys_op1_t)0x2u) /* {v0 == 1} (hard-wires 1 to {v0}). */

#define logsys_op1_TRU ((logsys_op1_t)0x3u) /* Always satisfied. */
#define logsys_op1_FAL ((logsys_op1_t)0x0u) /* Always failed. */

/* Some popular binary gates: */
#define logsys_op2_ZR0 ((logsys_op2_t)0x5u) /* {v0 == 0} (hard-wires 0 to {v0}). */
#define logsys_op2_UN0 ((logsys_op2_t)0xAu) /* {v0 == 1} (hard-wires 1 to {v0}). */

#define logsys_op2_ZR1 ((logsys_op2_t)0x3u) /* {v1 == 0} (hard-wires 0 to {v1}). */
#define logsys_op2_UN1 ((logsys_op2_t)0xCu) /* {v1 == 1} (hard-wires 1 to {v1}). */

#define logsys_op2_W01 ((logsys_op2_t)0x9u) /* {v0 == v1} (hard-wires {v0} to {v1}). */
#define logsys_op2_N01 ((logsys_op2_t)0x6u) /* {v0 == (! v1)} ({v0} is complement of {v1}). */

#define logsys_op2_G01 ((logsys_op2_t)0xBu) /* {v0 \geq v1} ({v1} implies {v0}). */
#define logsys_op2_L01 ((logsys_op2_t)0xDu) /* {v0 \leq v1} ({v0} implies {v1}). */

/* Some popular ternary gates: */
#define logsys_op3_ZR0 ((logsys_op3_t)0x55u) /* {v0 == 0} (hard-wires 0 to {v0}). */
#define logsys_op3_UN0 ((logsys_op3_t)0xAAu) /* {v0 == 1} (hard-wires 1 to {v0}). */

#define logsys_op3_ZR1 ((logsys_op3_t)0x33u) /* {v1 == 0} (hard-wires 0 to {v1}). */
#define logsys_op3_UN1 ((logsys_op3_t)0xCCu) /* {v1 == 1} (hard-wires 1 to {v1}). */

#define logsys_op3_ZR2 ((logsys_op3_t)0x0Fu) /* {v2 == 0} (hard-wires 0 to {v2}). */
#define logsys_op3_UN2 ((logsys_op3_t)0xF0u) /* {v2 == 1} (hard-wires 1 to {v2}). */

#define logsys_op3_AND ((logsys_op3_t)0x95u) /* {v0 == v1 /\ v2 } (logical AND gate). */
#define logsys_op3_ORR ((logsys_op3_t)0xA9u) /* {v0 == v1 \/ v2} (logical OR operator). */
#define logsys_op3_XOR ((logsys_op3_t)0x69u) /* {v0 == v1 ++ v2} (logical XOR gate). */
#define logsys_op3_EQV ((logsys_op3_t)0x96u) /* {v0 == (v1 == v2)} (logical EQV gate). */

#define logsys_op3_W01 ((logsys_op3_t)0x99u) /* {v0 == v1} (hard-wires {v0} to {v1}). */
#define logsys_op3_N01 ((logsys_op3_t)0x66u) /* {v0 == (! v1)} ({v0} is complement of {v1}). */

#define logsys_op3_W02 ((logsys_op3_t)0xA5u) /* {v0 == v2} (hard-wires {v0} to {v2}). */
#define logsys_op3_N02 ((logsys_op3_t)0x5Au) /* {v0 == (! v2)} ({v0} is complement of {v2). */

#define logsys_op3_W12 ((logsys_op3_t)0xC3u) /* {v1 == v2} (hard-wires {v1} to {v2}). */
#define logsys_op3_N12 ((logsys_op3_t)0x3Cu) /* {v1 == (! v2)} ({v1} is complement of {v2}). */

#define logsys_op3_ANY ((logsys_op3_t)0xFFu) /* Allows any assignment to {v0,v1,v2}. */
#define logsys_op3_NIL ((logsys_op3_t)0x00u) /* Rejects any assignment to {v0,v1,v2}. */

#define logsys_op3_FAN ((logsys_op3_t)0x81u) /* {v0 == v1} and {v1 == v2} (fan-out operator). */

bool_t logsys_op_is_valid(logsys_op_t op, int n);
  /* Returns TRUE iff {op} is a valid operation code for {n} variables. 
     Namely, if all bits are zero beyond the lowest {2^n} . */
     
bool_t logsys_op_depends(logsys_op_t op, int n, int j);
  /* TRUE iff equations on {n} variables with opcode {op} really depend
    on argument {j}. Namely, iff there is at least one combination of
    values for the other arguments such that EXACTLY ONE value for
    argument {j} will satisfy the equation.  */

bool_t logsys_op_is_functional(logsys_op_t op, int n, int j);
  /* TRUE iff equations on {n} variables with opcode {op} are functional
    for argument {j}.  Namely, iff for any assignment of values 
    to the other arguments, there is at EXACTLY ONE assignment for argument
    {j} that satisfies {op}. */

bool_t logsys_op_is_solvable(logsys_op_t op, int n, int j);
  /* TRUE iff equations on {n} variables with opcode {op} are always
    solvable for argument {j}. Namely, if for any assignment of values to the
    other arguments within those ranges, there is AT LEAST ONE
    assignment for argument {j} within its range that satisfies {op}. */

logsys_op_t logsys_build_args_mask(int n, int na, int ka[]);
   /* Build the set of all assignments of values to {n} args where
      only args {ka[0..na-1]} are varied, all other args are set to 0. */
      
logsys_op_t logsys_scramble_args_in_op(logsys_op_t op_old, int n_old, int n_new, int k[]);
  /* Modifies an op code {op_old} for {n_old} argument variables {v_old[0..n_old-1]}
    for a different set of argument variables {v_new[0..n_new-1]}.  
    
    The parameter {k} is a vector of {n_new} indices that defines the
    correspondence beyween the new and old arguments. Specifically, for
    each {j} in {0..n_new-1}, if {k[j]} is negative, the procedure
    assumes that the new argument {v_new[j]} is a variable distinct from
    all the old arguments and other new arguments. In that case the new
    operator will not constrain argument {j}. Otherwise, the procedure
    assumes that the new argument {v_new[j]} is the same variable as the
    old argument {v_old[k[j]]}.
    
    The {k} vector should not contain repeated non-negative values.

    The new operator tries to capture as much information as the old
    one. If some old variable is to be omitted in {op_new} (i.e., some
    integer {i} in the range {0..nold-1} does not appear in
    {k[0..nnew-1]}, then the valid set of assignments allowed by {op} is
    projected along the axis of that variable, which may discard some
    information. If one of the new arguments is a new variable
    ({k[j]<0}) then the new operator will not put any constraint on that
    variable. */
    
logsys_op_t logsys_unscramble_args_in_op(logsys_op_t op_new, int n_old, int n_new, int k[]);
  /* Roughly the inverse of {logsys_scramble_args_in_op}.  Namely,
    given an operator {op_new} on argument {n_new} variables {v_new[0..n_new-1]},
    returns the operator {op_old} for {n_old} argument variables {v_old[0..n_old-1]},
    assuming the correspondence between old and new variables
    is given by the vector {k[0..n_new-1]}, as in that procedure.
    
    If some old variable was omitted in {op_new} (i.e., some integer {i}
    in the range {0..nold-1} does not appear in {k[0..nnew-1]}, then the
    returned op will be indifferent to that variable. If one of the new
    arguments is a new variable ({k[j] < 0}) then {op_new} must not
    depend on that variable. */
    
int logsys_check_arg_dependency(logsys_op_t op, int n, int na, int ka[], int nb, int kb[]);
  /* Checks whether the equation {op} on {n} arguments is the product of two 
    separate equations depending on two given subsets of arguments.  The arguments of {op}
    with indices {ka[0..na-1]} are the first subset, those with indices {kb[0..nb-1]}
    are the second subset.  The lists {ka[0..na-1]} and {kb[0..nb-1]}
    are assumed to be disjoint, free from repetitions, and cover
    {0..n-1}, hence one must have {na+nb == n}.
    
    Returns {-1} {op} splits as specified (in particular, if {na==0} or
    {nb==0}). Otherwise returns an integer {ib} such that the admissible
    values of argument {kb[ib]} depend on the values of arguments
    {ka[0..na-1]}. */

void logsys_split_op
  ( logsys_op_t op, 
    int n, 
    logsys_op_t *opaP, 
    int *naP, 
    int ka[], 
    logsys_op_t *opbP, 
    int *nbP, 
    int kb[]
  );
  /* Checks whether the equation {op} on {n} arguments can be split into
    two equations {opa,opb} respectively with {na,nb} arguments,
    mutually disjoint. 
    
    If the call succeeds, stores the two opcodes into {*opaP,*opbP} and their
    arg counts into {*naP,*nbP}. Also stores into {ka[0..na-1]} the
    indices of the arguments of {op} that are used by {opa}; and
    similarly for {kb[0..nb-1]}.  The set {ka} will be a minimal set
    that has argument 0.
    
    If the call fails, returns a trivial split with {na=n}, {opa=op},
    {ka[i] = i}, {nb=0}, {opb=TRU}. */
    
double logsys_op_entropy_reduction(logsys_op_t op, int n, int j);
  /* Assumes {op} is an {n}-ary operator. Returns the sum of the
    expected reduction in the entropy of the other arguments of {op}
    when argumment {j} is set to either 0 or 1 with equal probability. */

/* SYSTEMS, EQUATIONS, VARIABLES */

logsys_t *logsys_new_system(void);
  /* Creates a new system, initially with no variables and no equations. */

logsys_va_t *logsys_add_variable(logsys_t *S);
  /* Adds a new variable to the system {S}, initially with no users, and returns it. */

logsys_eq_t *logsys_add_equation(logsys_t *S, logsys_op_t op, int n, logsys_va_t *va[], logsys_eq_t *preq);
  /* Adds a new equation to system {S}, consisting of the operation {op} and
    the argument variables {va[0..n-1]}, in that order.  The variables must belong to
    the system {S} and must be non-null and distinct. 
    
    If {preq} is not NULL, it must point to an equation of {S}; the new
    equation is inserted just after it in the list of equations of {S}.
    Otherwise it is inserted after the last equation previously
    added. */

void logsys_delete_equation(logsys_t *S, logsys_eq_t *eq);
  /* Deletes equation {eq} from system {S}, and reclaims its space.
    Do not use {*eq} after calling this procedure! */

int logsys_num_args(logsys_eq_t *eq);
  /* Returns the number of operands (arguments) of equation {eq}. */

int logsys_find_arg(logsys_eq_t *eq, logsys_va_t *va);
  /* Returns the index {k} of the operand of {eq} that is the variable {va},
    or {-1} if {eq} does not use {va}. */

bool_t logsys_narrow_ranges_by_equation(logsys_eq_t *eq, bool_t lo[], bool_t hi[]);
  /* Let{n} {va[0..n-1]} be the argument variables of {eq}.
    Considers the set {A} of all assignments {x[0..n-1]} of Boolean
    values for those variables that will satisfy the equation, and 
    are such that {lo[id[j]] <= x[j] <= hi[j]} for {j} in {0..n-1}. 
    Then narrows each range {lo[j]..hi[j]} as much
    as possible while still enclosing {A}.
    
    The procedure returns FALSE if the equation cannot be satified at all
    with the given range restrictions, and TRUE otherwise. */

logsys_va_t *logsys_replace_arg(logsys_eq_t *eq, int j, logsys_va_t *va);
  /* Replaces argument {j} of equation {eq} by variable {va}. The
    variable {va} may be NULL. Both {eq} and {va} (if not NULL) must
    belong to the same system. */

/* SOLUTIONS 

  For the procedures in this section, the array parameter {va[]} is the 
  list of all variables in the system {S}, indexed by their IDs.
  Namely, {va[0..nv-1]} should be all variables of {S} (and possibly some NULL entries);
  and, for each {i} where {va[i]} is not NULL, the variable {va[i]} must have
  ID number {i}. The array parameter {eq[]} is similarly defined
  for equations. */

int logsys_count_variables(logsys_t *S);
int logsys_count_equations(logsys_t *S);
  /* Count the number of equations and variables in the system. */

int logsys_get_variable_table_size(logsys_t *S);
int logsys_get_equation_table_size(logsys_t *S);
  /* Return the maximum variable or equation ID in the system, plus one.
    If there are none, returns 0. */
  
void logsys_get_variables(logsys_t *S, int *nvP, logsys_va_t ***vaP);
  /* Builds the list {va[0..nv-1]} all variables of {S}.  The entry
    {va[ia]} will be NULL if there is no variable with that ID in the system.
   
    The procedure allocates the vector {va) with sufficient size.
    Its address is returned into {*vaP} and its size in {*nvP}. */

void logsys_get_equations(logsys_t *S, int *neP, logsys_eq_t ***eqP);
  /* Builds the list {eq[0..ne-1]} all equations of {S}.  The entry
    {eq[ie]} will be NULL if there is no equation with that ID in the system.
   
    The procedure allocates the vector {eq) with sufficient size.
    Its address is returned into {*eqP} and its size in {*neP}. */
 
void logsys_sort_variables_for_solver
  ( int heur, 
    logsys_t *S, 
    int nv, 
    logsys_va_t *va[], 
    int vix[], 
    bool_t guess[]
  );
  /* Returns in {vix[0..nv-1]} a permutation of {0..nv-1}, the indices
    of the variables {va[0..nv-1]} in the best order for the
    backtracking solver. Namely the solver should try to fix the values
    of {va[vix[i]]} in order of increasing {i}. In particular, the
    indices {vix[i]} where {va[ix[i]]} is NULL are all at the end of the
    {vix} array. 
    
    If {guess} is not NULL, also stores into {guess[0..nv-1]} a suggested 
    first guess for the variables {va[0..nv-1]}. */

void logsys_sort_variables_by_score(int nv, logsys_va_t *va[], double score[], int sign, int vix[]);
  /* Assumes that {score[ia]} is a numeric score for variable {va[ia]},
    for {ia} n {0..nv-1}. Stores into {vix[0..nv]} a permutation of
    {0..nv} such that {sign*score[vix[k]]} is increasing with increasing
    {k}. If {va[ia]} is NULL, pretends that {sign*score[ia]} is {+oo},
    so that those indices will be moved to the end of the array. */

bool_t logsys_find_solution(logsys_t *S, int nv, logsys_va_t *va[], int vix[], bool_t lo[], bool_t hi[], bool_t guess[]);
  /* Assumes that {va[0..nv-1]} is the list of all variables of {S}.
    Tries to find a solution for the system {S} that is bracketed
    between {lo[0..nv-1]} and {hi[0..nv-1]}.
    
    More precisely, looks for a Boolean vector {x[0..nv-1]} such that
    (1) all equations in {S} are satisfied when each non-null variable
    {va[i]} is set to value {x[i]}, and (2) {lo[i] <= x[i] <= hi[i]} for
    all {i}. If it finds such an assignment, sets {lo[i]=hi[i]=x[i]} for
    all {i}, and returns TRUE. Otherwise, returns FALSE, leaving the
    vectors {lo} and {hi} unchanged.
    
    If {va[i]} is NULL, the entries {lo[i]} and {hi[i]} are irrelevant,
    and are not modified. 
    
    Uses a propagation-with-backtracking method. If {vix} is not NULL,
    it must contain a permutation of the indices {0..nv-1}, which
    defines the order in which variables are tentatively fixed during
    backtracking. If {vix} is null, the variables to be tentatively 
    fixed are picked at random.  If {guess} is not NULL, it must
    be a vector of {nv} booleans; when tentatively fixing an 
    indeterminate variable {va[iv]}, the procedure will try 
    first {guess[iv]} then {!guess[iv]}.  Otherwise it will pick
    the first guess arbitrarily. */

bool_t logsys_equation_is_satisfied(logsys_eq_t *eq, int nv, logsys_va_t *va[], bool_t x[]);
  /* Returns TRUE iff the variable assignment {x[0..nv-1]} for the system
    variables {va[0..nv-1]} satisfies the equation {eq}. 
    Ignores the values of {x[i]} for any {i} such that 
    {va[i]} is NULL or is not an argument of {eq}. */

bool_t logsys_is_solution(logsys_t *S, int nv, logsys_va_t *va[], bool_t x[]);
  /* Returns TRUE iff the variable assignment {x[0..nv-1]} for the system
    variables {va[0..nv-1]} satisfies all the eequations of {S}. 
    Ignores the values of {x[i]} for any {i} such that 
    {va[i]} is NULL or is not an argument of any equation. */

typedef bool_t logsys_solution_proc_t(int nv, logsys_va_t *va[], bool_t x[]);
  /* Type of a procedure suitable to be passed to {logsys_enum_solutions}. */

void logsys_enum_solutions(logsys_t *S, int nv, logsys_va_t *va[], bool_t lo[], bool_t hi[], logsys_solution_proc_t *proc);
  /* Assumes that {va[0..nv-1]} is the list of all variables of {S}. Exhaustively enumerates 
    all assignments {x[0..nv-1} for those variables, and {proc(nv, va, x)} for every
    assignment that satisfies all equations of {S}.  Stops if {proc} return FALSE.
    If {lo,hi} are not NULL, varies each variable {va[i]} 
    only in the range {lo[i]..hi[i]} */

/* SIMPLIFICATION AND MODIFICATION */

bool_t logsys_try_to_simplify_equation(logsys_eq_t *eq);
  /* Tries to modify the equation {eq}, by removing any argument variables
    that the equation does not depend upon. Returns TRUE if the equation
    was modified. */

void logsys_try_to_combine_equations(logsys_eq_t **aP, bool_t *amodP, logsys_eq_t **bP, bool_t *bmodP);
  /* Analyzes equations {**aP} and {**bP} together, and tries to simplify each one 
    in view of constraints imposed by the other.  Also tries to merge them into a single equation. 
    If it succeeds, it may delete one or both of them.  
    
    Assumes that the equations are simplified when taken in isolation.  Does nothing if the 
    arguments of the two equations are disjoint sets of variables.
    
    If {**aP} is either modified or deleted, {*amodP} is set to TRUE. Moreover, if {**aP} gets deleted, 
    {*aP} gets set to NULL.  The same applies to {bmodP} and {bP}. */

bool_t logsys_try_to_split_equation(logsys_eq_t *eq, int *nstP, logsys_eq_vec_t *steq);
  /* Analyzes equation {eq}, and tries to split it into two or more
    separate equations on disjoint non-empty sets of variables. 

    In particular, if the equation does not depend on some argument,
    removes that variable from the equation. If the the equation is
    always true (accepts any value assignment to its args) or trivially
    false (rejects every value assignment to its args), changes it into
    a {logsys_op0_TRU} or {logsys_op0_FAL} with zero arguments.
    
    Any non-trivial factor that is split off is inserted in the system
    just BEFORE {eq} in the list of equations of {S}.
    
    Returns TRUE iff {eq} is modified in any way (split or simplified).
    In that case, if {steq} is not null, stores all factors into
    {steq.e} starting at position {*nstP} and increments {*nstP}
    accordingly. */

void logsys_condense_and_split(logsys_t *S);
  /* Tries to simplify the system {S} by iteratively condensing equations that share
    variables and splitting those that factor into multiple 
    factors with disjoint argument variables. */

void logsys_normalize(logsys_t *S, int *nu, logsys_va_t ***ua);
  /* Modifies the system {S} so that all its variables are
    functionally determined by a certain set of independent variables
    {(*ua)[0..*nu-1]}, new or old, chosen by the procedure.
    The vector {*ua} is allocated by the procedure. */

/* MISCELLANEOUS UTILITIES */

void logsys_get_bits(uint64_t X, int n, bool_t x[]);
  /* Sores into {x[0..n-1]} the bits of the base-2 representation 
    of the unsigned integer {X}.  Bit {x[0]} is the units bit.
    Fails if {X} is {2^n} or greater. */
    
void logsys_sort_ints(int n, int k[]);
  /* Sorts the integers {k[0..n-1]} in increasin order, by insertion sort. */

logsys_op_t logsys_build_args_mask(int n, int na, int ka[]);
  /* Returns an op code with {n} arguments that accepts exactly those
    {2^na} argument value combinations where the arguments with indices
    {ka[0..na-1]} have arbitrary values while the other arguments are
    set to zero. In other words, returns all corners of an
    {(n-na)}-dimensional face of the Boolean {n}-cube, stabbed by axes
    {ka[0..na-1]}. */   

double logsys_entropy(double m0, double m1);
  /* Entropy of a variable that is FALSE with probability {m0/(m0+m1)}, TRUE
   with probability {m1/(m0+m1)}.  Returns 1 if {m0==m1==0}. */

int logsys_ones(uint64_t X);
  /* Number of 1 bits in {X}. */

/* VARIABLE AND EQUATION ID NUMBERS */

typedef uint64_t logsys_eq_id_t;
  /* A serial identifier of an equation. */
  
typedef uint64_t logsys_va_id_t;
  /* A serial identifier of a variable. */

logsys_va_id_t logsys_variable_id(logsys_va_t *va);
  /* Returns an integer (the /ID number/) that uniquely identifies a
    variable in a system. Variable ID numbers are assigned serially,
    starting with 1, as variables are created. */

logsys_eq_id_t logsys_equation_id(logsys_eq_t *eq);
  /* Returns an integer (the /ID number/) that uniquely identifies a
    equation in a system. Equation ID numbers are assigned serially,
    starting with 1, as equations are created. */

/* PRINTOUT */

void logsys_print_stats(FILE *wr, logsys_t *S, int nst, logsys_eq_vec_t *steq);
  /* Prints equation counts by arity for {S} (if not NULL)
    and for {steq.e[0..nst-1]} (if not NULL). */

void logsys_print_equation(FILE *wr, char *pre, logsys_eq_t *eq, char *suf);
  /* Prints the equation {eq} to {wr}, bracketed by the {pre}
    and {suf} strings. (If either string is NULL, it is omitted). The argument variables are printed
    with argument 0 at the RIGHT. */

void logsys_print_equation_id(FILE *wr, char *pre, logsys_eq_id_t id, int d, char *suf);
  /* Prints the equation ID number to {wr} using {d} decimal digits,
    zero-filled, bracketed by the {pre} and {suf} strings. (If either
    string is NULL, it is omitted). */

void logsys_print_variable_id(FILE *wr, char *pre, logsys_va_id_t id, int d, char *suf);
  /* Prints the variable ID number to {wr} using {d} decimal digits,
    zero-filled, bracketed by the {pre} and {suf} strings. (If either
    string is NULL, it is omitted). */

void logsys_print_op(FILE *wr, char *pre, int n, logsys_op_t op, bool_t verbose, char *suf);
  /* Prints to {wr} the opcode {op}, assumed to be on {n} variables,
    bracketed by the {pre} and {suf} strings. (If either string is NULL,
    it is omitted). If {verbose} is true, prints the truth table in the
    format of {logsys_print_ops} (one assignment per line). If {verbose}
    is false, prints only the opcode's conventional name or "o{HEX}"
    where {HEX} is the opcode in hexadecimal. */

void logsys_print_ops(FILE *wr, char *pre, int n, int m, uint64_t op[], char *suf);
  /* Prints to {wr} the truth-tables of the opcodes {op[0..m-1]}, assumed to be all on the same {n} variables,
    bracketed by the {pre} and {suf} strings. (If either string is NULL, it is omitted).
    Prints one assignment {x[0..n-1]} per line (with {x[0]} at the RIGHTMOST
    end), and then the values of {op[0..m-1]} for that assignment. */

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
  );
  /* Prints to {wr} local ersatz names and (optionally) values of {n}
    parameters of an operation (with argument 0 at the RIGHTMOST end).
    The names are prefixed by {pre}, separated by {sep} and terminated
    by {suf}. (If any of these strings is NULL, it is omitted.)
    
    The names are "{vtag_old}{k[i]}" if {k[i]} is non-negative, and
    "{vtag_new}{i}" otherwise, where {i} ranges from {n-1} down to 0. If
    {k} is NULL, assumes {k[i]==i}, so the names are "{vtag_old}{i}".
    
    If {x} is not NULL, assumes {x[i]} is the boolean value of argument
    {i}, and prints "={x[i]}" (as "0" or "1") after the argument name. */

void logsys_print_system(FILE *wr, char *head, logsys_t *S, char *foot);
  /* Prints the system {S} to {wr}, preceded by a line containing the
    string {head} and followed by a line containing {foot}. If either
    string is NULL or empty, the corresponding line is omitted. */


void logsys_print_assignments(FILE *wr, char *pre, int n, bool_t x[], char *suf);
  /* Prints the boolean values {x[0..n-1]} as a binary string (with {x[0]} 
    at the RIGHTMOST end). If {n} is 64 or less, also prints the decimal and hex equivalents. The
    delimiter strings {pre} and {suf} are omitted if NULL. */

void logsys_print_ranges(FILE *wr, char *pre, int n, bool_t lo[], bool_t hi[], char *suf);
  /* Prints the boolean ranges {lo[0..n-1]--hi[0..n-1]} as a binary
    string (with {lo[0],hi[0]} at the RIGHTMOST end), using '*' for bits
    with indeterminate value and '!' for bits with empty range. 
    
    If {n} is 64 or less and all non-null variables have singleton
    ranges (i.e. all bits are definitely 0 or 1), also prints the
    decimal equivalent. The delimiter strings {pre} and {suf} are
    omitted if NULL. */


void logsys_print_var_assignments(FILE *wr, char *pre, int n, logsys_va_t *va[], bool_t x[], char *suf);
  /* Prints the boolean values assigned to variables {va[0..n-1]} in
    RIGHT to LEFT order as in {logsys_print_assignments}, skipping over
    any NULL entries of {va[0..n-1]}. Assumes that {x[k]} is the value
    assigned to the variable whose ID number is {k}, for every variable
    in the system; and that {va[0..nv-1]} is an arbitrary list of {n}
    variables or NULLs, in any order. */

void logsys_print_var_ranges(FILE *wr, char *pre, int n, logsys_va_t *va[], bool_t lo[], bool_t hi[], char *suf);
  /* Prints the bounds of {va[0..n-1]} in RIGHT to LEFT order as in
    {logsys_print_assignments}, skipping over any NULL entries of
    {va[0..n-1]}. Assumes that {lo[k]} and {hi[k]} are the lower and
    upper bounds for the boolean value of the variable whose ID number
    is {k}, for every variable in the system; and that {va[0..nv-1]} is
    an arbitrary list of {n} variables or NULLs, in any order. */


void logsys_print_arg_ranges(FILE *wr, char *pre, logsys_eq_t *eq, bool_t lo[], bool_t hi[], char *suf);
  /* Performs {logsys_print_var_ranges(wr, per, n, va, lo, hi, suf)}
    where {va[0..n-1]} are the argument variables of {eq}. */

#endif

