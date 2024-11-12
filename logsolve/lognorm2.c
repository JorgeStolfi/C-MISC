/* Last edited on 2012-12-20 02:21:59 by stolfilocal */

#include <stdint.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

#include <bool.h>
#include <sign.h>
#include <affirm.h>

#include <logsys.h>

void lognorm2_print_op3_table(void);
  /* Tabulates all the three-variable opcodes and their characteristics: */

void lognorm2_print_op3_op3_factoring_table(void);
  /* Prints a table of all pairs of three-variable equations
    {opA(a,u,v);opB(a,x,y)} and tries to find for each pair an
    equivalent combination {opC(a,p,q);opD(u,x,v,y)} where {p,q} are
    two of the variables {u,x,v,y}. */

bool_t lognorm2_factor_eq_5_3_4(logsys_op5_t op5, int *jpP, int *jqP, logsys_op3_t *op3P, logsys_op4_t *op4P);
  /* Given an equation {op5} on five variables {a,u,x,v,y}, tries to
    find a system of two equations with the form {op3(a,p,q);
    op4(u,x,v,y)} that is equivalent to {op5(a,u,x,v,y)}, where {p,q}
    are two variables among {u,x,v,y}.

    If this search succeeds, stores the relevant opcodes in {*op3P} and
    {*op4P}, and stores in {*jpP} and {*jqP} the positions of {p} and {q}
    in the variable list {u,x,v,y} and returns TRUE. If the search
    fails, returns FALSE, leaving {*op3P,*op5P,*jpP,*jqP} undefined.
    
    The opcodes {op5} and {op4} are interpreted in a way similar to
    {logsys_op3_t}. Namely, {op5 & (1 << x)} is 1 iff the equation
    {op5(a,u,x,v,y)} is satisfied when the bits of {x} are assigned to
    {a,u,x,v,y} (with the lowest bit going to {a}).
    
    The indices {*jpP,*jqP} will be distinct, in the range {0..3},
    with {*jpP < *jqP}, indicating a position in the list {u,x,v,y},
    FROM LFET TO RIGHT. So, for example, {*jpP} will be 0 if {p} is {u},
    and {3} is {p} is {y}. The procedure does not consider solutions
    with {p==NULL} or {q==NULL} or {p==q}, but it may return
    an opcode {op3} that ignores some of its arguments. */

bool_t lognorm2_factor_eq_5_3_4_aux(logsys_op5_t op5, int jp, int jq, logsys_op4_t op4, logsys_op3_t *op3P);
  /* Similar to {lognorm2_factor_eq_5_3_4}, but searchs only for the
    opcode {op3} of the first equation {op3(a,p,q)}, given the opcode
    {op4} of the second equation {op4(u,x,v,y)} and the identities of
    the arguments {p} and {q} (specified by parameters {jp} and {jq},
    each in the range {0..4}, with the same conventions as before).
    
    To maximize the chances of finding a solution, equation {op4} should
    be as restrictive as possible; but it must not forbid any assignment
    to {u,x,v,y} that is allowed by {op5}. */
    
logsys_op_t lognorm2_choose_op(logsys_op_t op_lo, logsys_op_t op_hi);
  /* Chooses an equation opcode {op}, for up to six variables,
    that accepts all assignments that {op_lo} accepts,
    and rejects all assignments that {op_hi} rejects.
    
    The procedure fails if the bounds are incompatible, that is,
    if {op_lo} requires the acceptance of an assignment that
    is forbidden by {op_hi}.  In other words, {op_hi | op_lo} 
    must be equal to {op_lo}. */

void lognorm2_factor_eq_5_3_4_check(logsys_op5_t op5, int jp, int jq, logsys_op3_t op3, logsys_op4_t op4);
  /* Checks whether the results returned by a successful call of {lognorm2_factor_eq_5_3_4} 
     really satisfy the requirements. */

#define ONE64 ((uint64_t)1)
#define ALL64 ((uint64_t)(-1LL))

int main(int argc, char **argv)
  {
    lognorm2_print_op3_table();
    lognorm2_print_op3_op3_factoring_table();
    
    return 0;
  }

void lognorm2_print_op3_table(void)
  {
    fprintf(stderr, "--- opcode properties -----------------------------------------\n");
    uint16_t iop;
    for (iop = 0; iop < 256; iop++)
      { logsys_op3_t op3 = (logsys_op3_t)iop;
        fprintf(stderr, "opcode %03u = %02X", iop, iop);
        logsys_print_op(stderr, " = ", 3, op3, FALSE, "(x,y,z)");
        fprintf(stderr, " depends on ");
        int j;
        for (j = 0; j < 3; j++)
          { char c = (logsys_op_depends(op3, 3, j) ? ("xyz")[j] : '-');
            fputc(c, stderr);
          }
        fprintf(stderr, " is solvable for ");
        for (j = 0; j < 3; j++)
          { char c = (logsys_op_is_solvable(op3, 3,  j) ? ("xyz")[j] : '-');
            fputc(c, stderr);
          }
        fprintf(stderr, " is functional for ");
        for (j = 0; j < 3; j++)
          { char c = (logsys_op_is_functional(op3, 3, j) ? ("xyz")[j] : '-');
            fputc(c, stderr);
          }
        fprintf(stderr, "\n"); 
      }
  }

void lognorm2_print_op3_op3_factoring_table(void)
  { 
    uint16_t iop3A;
    for (iop3A = 0; iop3A < 256; iop3A++)
      { logsys_op3_t op3A = (logsys_op3_t)iop3A;
        uint16_t iop3B;
        for (iop3B = 0; iop3B <= iop3A; iop3B++)
          { logsys_op3_t op3B = (logsys_op3_t)iop3B;
            bool_t solva = logsys_op_is_solvable(op3A, 3, 0);
            bool_t solvb = logsys_op_is_solvable(op3B, 3, 0);
            if (solva && solvb)
              { /* Try to normalize the equations {op3A(a,u,v); op3B(a,x,y)}. */
                fprintf(stderr, "---------------------------------------------------------------\n");
                fprintf(stderr, "normalizing {");
                logsys_print_op(stderr, NULL, 3, op3A, FALSE, "(a,u,v)");
                fprintf(stderr, "; ");
                logsys_print_op(stderr, NULL, 3, op3B, FALSE, "(a,x,y)");
                fprintf(stderr, "}");
                
                /* Find the set {op5} of all 5-bit tuples that satisfy the two equations: */
                logsys_op5_t op5 = 0;
                uint8_t auxvy;  /* An assignment to {a,u,x,v,y}, with {y} as the lowest bit. */
                for (auxvy = 0; auxvy < 32; auxvy++)
                  { /* Check whether the bit combination {auxvy} satsfies both eqs: */
                    int32_t auv = (auxvy & 3LU) | ((auxvy >> 1) & 4LU);
                    int32_t axy = ((auxvy >> 2) & 4LU) | ((auxvy >> 1) & 2LU) | (auxvy & 1LU);
                    if ((((1LU << auv) & op3A) != 0) && (((1LU << axy) & op3B) != 0))
                      { /* Add the valid tuple {auxvy} to the set {op5}. */
                        op5 |= (1LU << auxvy);
                      }
                  }
                fprintf(stderr, " <==> ");
                logsys_print_op(stderr, NULL, 5, op5, FALSE, NULL);
                
                /* Look for a system {op3(a,u,v);op4(u,x,v,y)} equivalent to {op3A(a,u,v); op3B(a,x,y)}. */
                logsys_op3_t op3;
                logsys_op4_t op4;
                int jp; /* Index of {p} among {u,x,v,y}. */
                int jq; /* Index of {q} among {u,x,v,y}. */
                bool_t ok = lognorm2_factor_eq_5_3_4(op5, &jp, &jq, &op3, &op4);
                
                if (ok)
                  { fprintf(stderr, " --> {");
                    logsys_print_op(stderr, NULL, 3, op3, FALSE, NULL);
                    fprintf(stderr, "(a,%c,%c)", ("uxvy")[jp], ("uxvy")[jq]);
                    fprintf(stderr, "; ");
                    logsys_print_op(stderr, NULL, 4, op4, FALSE, NULL);
                    fprintf(stderr, "(u,x,v,y)");
                    fprintf(stderr, "}\n");
                    lognorm2_factor_eq_5_3_4_check(op5, jp, jq, op3, op4);
                    
                  }
                else
                  { fprintf(stderr, "  ** no {a}-factoring found\n"); }
                
                /* Find ternary equations equivalent to {op4(u,x,v,y)}: */
                /* fprintf(stderr, "  ** op4 factoring not implemented yet\n"); */
                fprintf(stderr, "---------------------------------------------------------------\n");
              }
          }
      }
  }

bool_t lognorm2_factor_eq_5_3_4(logsys_op5_t op5, int *jpP, int *jqP, logsys_op3_t *op3P, logsys_op4_t *op4P)
  {
    /* Compute the opcode {op4} of the restrictions on {u,x,v,y}: */
    /* Enumerate all assignments to {a,u,x,v,y} allowed by {op5}: */
    logsys_op4_t op4 = 0x0000; /* Assignments to {u,x,v,y} allowed by {op5}, for now. */
    
    /* !!! There are faster ways to do this !!! */
    uint8_t auxvy;  /* An assignment to {a,u,x,v,y}, with {a} as the lowest bit. */
    for (auxvy = 0; auxvy < 32; auxvy++)
      { if (((1LU << auxvy) & op5) != 0)
          { /* Assignment {auxvy} is allowed by {op5}. */
            uint8_t uxvy = (auxvy >> 1); /* The assignment to {u,x,v,y}, with {u} as the lowest bit. */
            /* Add the sub-assignment {uxvy} to {op4}: */
            op4 |= (1LU << uxvy);
          }
      }

    /* Now choose the two arguments {p.q} of the equation {op3(a,p,q)}: */
    int kp; /* Index of {p} among {u,x,v,y}. */
    int kq; /* Index of {q} among {u,x,v,y}. */
    /* Enumerate all unordered pairs of variables among {u,x,v,y}, in LEFT TO RIGHT order: */
    for (kq = 0; kq < 4; kq++)
      { for (kp = 0; kp < kq; kp++)
          { /* See if {a} can be constrained by just {p} and {q}: */
            logsys_op3_t op3;
            bool_t ok = lognorm2_factor_eq_5_3_4_aux(op5, kp, kq, op4, &op3);
            if (ok) 
              { /* Found a split: */
                (*jpP) = kp; 
                (*jqP) = kq;
                (*op3P) = op3;
                (*op4P) = op4;
                return TRUE;
              }
          }
      }
    /* There is no way to factor {op5} with a single ternary equation on {a}: */
    return FALSE;
  }

bool_t lognorm2_factor_eq_5_3_4_aux(logsys_op5_t op5, int jp, int jq, logsys_op4_t op4, logsys_op3_t *op3P)
  {
    /* Note that {op3} is not unique when {op5} forbids some assignments to {p,q}. */
    /* Therefore we compute two ternary opcodes {op3_lo,op3_hi} that bracket {op3} bitwise. */
    logsys_op3_t op3_lo = 0x00u; /* Assignments to {apq} which MUST be allowed by {op3}, for now. */
    logsys_op3_t op3_hi = 0xFFu; /* Assignments to {apq} which MAY be allowed by {op3}, for now. */
    
    /* Enumerate all possible assignments {auxvy} to {a,u,x,v,y}: */
    uint8_t auxvy;
    for (auxvy = 0; auxvy < 32; auxvy++)
      { /* Extract from {auxvy} the values assigned to {a,p,q}: */
        uint8_t a = auxvy & 1LU;
        uint8_t p = (auxvy >> (1 + jp)) & 1LU;
        uint8_t q = (auxvy >> (1 + jq)) & 1LU;
        uint8_t apq = a | (p << 1) | (q << 2); /* The values of {a,p,q}, with {a} as the lowest bit. */
        /* Extract from {auxvy} the values assigned to {u,x,v,y}: */
        uint8_t uxvy = (auxvy >> 1); /* The assignment to {u,x,v,y}, with {u} as the lowest bit. */
        /* Now update {op3_lo,op3_hi} accordingly: */
        if (((1LU << auxvy) & op5) != 0)
          { /* Assignment {auxvy} is allowed by {op5}. */
            assert(((1LU << uxvy) & op4) != 0); /* Since {op4} is supposed to be compatible with {op5}. */
            /* Then {op3} must allow {apq}: */
            op3_lo |= (1LU << apq);
          }
        else if (((1LU << uxvy) & op4) != 0)
          { /* Assignment {auxvy} is rejected by {op5}, but {uxvy} is allowed by {op4}. */
            /* Then {op3} must reject {apq}: */
            op3_hi &= (~ (1LU << apq));
          }
        else
          { /* Assignment {auxvy} is rejected by {op5}, but {uxvy} is rejected by {op4} anyway. */
            /* In this case, it does not matter what {op3} says about {apq}. */
          }
        
        /* Check whether the implications of {op5(auxvy)} on {op3(apq)} conflict with previous deductions: */
        if ((op3_lo & (~ op3_hi)) != 0)
          { /* The behavior of {op3} for assignment {apq} cannot be defined consistently with {op5,op4}: */
            return FALSE;
          }
      }
      
    /* Now choose an operation: */
    (*op3P) = lognorm2_choose_op(op3_lo, op3_hi);
    return TRUE;
  }
    
uint64_t lognorm2_choose_op(uint64_t op_lo, uint64_t op_hi)
  {
    /* Argument consistency check: */
    demand((op_lo | op_hi) == op_hi, "inconsistent op bounds");
    
    /* Scan the assignments bit by bit: */
    uint64_t op = 0; /* The chosen operation, so far. */
    uint8_t arg;     /* An assignment to up to 6 variables. */
    uint64_t arg_mk; /* Position of the {arg}-acceptance bit in {op}. */
    for (arg = 0, arg_mk = ONE64; (arg < 64) && (arg_mk <= op_hi);  arg++, arg_mk <<= 1)
      { if ((op_lo & arg_mk) != 0)
          { /* There is no choice, {op3} must accept {arg}: */
            op |= arg_mk; 
          }
        else if ((op_hi & arg_mk) == 0)
          { /* There is no choice, {op3} must reject {arg}: */
            op &= (~ arg_mk); 
          }
        else
          { /* We can choose whether to accept or reject {arg}. */
            /* We try to maximize the chances of {op} being functional for {v0}. */
            uint8_t oth = (arg ^ 1U); /* Assignment {arg} with variable {v0} negated. */
            if ((op & (ONE64 << oth)) == 0)
              { /* We do not accept {oth}, so we choose to accept {arg}: */
                op |= arg_mk; 
              }
            /* !!! However, other choices may allow for better {op}. !!! */
          }
      }
    return op;
  }

void lognorm2_factor_eq_5_3_4_check(logsys_op5_t op5, int jp, int jq, logsys_op3_t op3, logsys_op4_t op4)
  {
    /* Enumerate all possible assignments {auxvy} to {a,u,x,v,y}: */
    uint8_t auxvy;
    for (auxvy = 0; auxvy < 32; auxvy++)
      { /* Extract from {auxvy} the values assigned to {a,p,q}: */
        uint8_t a = auxvy & 1LU;
        uint8_t p = (auxvy >> (1 + jp)) & 1LU;
        uint8_t q = (auxvy >> (1 + jq)) & 1LU;
        uint8_t apq = a | (p << 1) | (q << 2); /* The values of {a,p,q}, with {a} as the lowest bit. */
        /* Extract from {auxvy} the values assigned to {u,x,v,y}: */
        uint8_t uxvy = (auxvy >> 1);
        /* Now compare {op5(a,u,x,v,y)} with {op3(a,p,q) /\ op4(u,x,v,y)}: */
        bool_t op5_sat = (((1LU << auxvy) & op5) != 0);
        bool_t op4_sat = (((1LU << uxvy) & op4) != 0);
        bool_t op3_sat = (((1LU << apq) & op3) != 0);
        if (op5_sat != (op4_sat && op3_sat))
          { /* Oops! ... */
            fprintf(stderr, "** bug in result of {lognorm2_factor_eq_5_3_4}\n");
            fprintf(stderr, "op5(%02X) = %d", auxvy, op5_sat);
            fprintf(stderr, "op4(%01X) = %d", uxvy, op4_sat);
            fprintf(stderr, "op3(%01X) = %d", apq, op3_sat);
            fprintf(stderr, "\n");
            demand(FALSE, "aborted");
          }
      }
 }
