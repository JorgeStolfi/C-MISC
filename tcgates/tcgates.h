#ifndef tcgates_H
#define tcgates_H

/* Tools for test-and-complement gates and circuits. */
/* Last edited on 2008-02-12 11:58:02 by stolfi */

#define tcgates_H_COPYRIGHT "Copyright © 2006 by the State University of Campinas (UNICAMP)"
#define tcgates_H_AUTHORS "Created 2006-may-20 by Jorge Stolfi, IC-UNICAMP"

#include <bool.h>

#include <stdio.h>
#include <stdint.h>
#include <values.h>

/* BIT STRINGS. */
  
#define m_max 3
  /* Max number of input bits supported by this implementation. */

#define n_max 8
  /* Max number of TC circuit wires supported in this implementation. */

typedef uint8_t word_t;
  /* A string of up to {n_max} bits. Usually the low-order {n} bits
    are important, for some {n}. */

typedef int32_t word_ct_t; 
  /* A count of words, including -1 and {2^n_max}. */

word_ct_t word_count(int n);
  /* Number of {n}-bit words, namely {2^n}. */

void word_print(FILE *wr, char *pre, int n, word_t x, bool_t bin, bool_t dec, char *suf);
  /* Prints the lowest {n} bits of word {x} to {wr} in the format
    "xxx=X", where "xxx" is the binary version and "X" the decimal
    equivalent. The "xxx=" part is printed only if {bin} is true, and
    the "=X" part is printed only if {dec} is TRUE. The binary version
    is reversed, with bit 0 (the lowest-order one) it at the LEFT end.
    The whole is sandwiched between strings {pre} and {suf}. */

void word_table_print_v(FILE *wr, int ind, int m, int n, word_ct_t ysz, word_t y[]);
  /* Prints the lowest {n} bits of words {y[0..ysz-1]} in vertical
    format. Each line has format "XXX=X -> YYY=Y" where "XXX=X" is the
    index (in reversed binary and decimal), and each "YYY=Y" is the
    value (ditto). Prints the lowest {m} bits of the argument and the 
    lowest {n} bits of the value.  Lines are indented by {ind} spaces. */

void word_table_print_h(FILE *wr, char *pre, int m, int n, word_ct_t ysz, word_t y[], char *suf);
  /* Prints the lowest {n} bits of words {y[0..ysz-1]} in horizontal
    format, "XXX=X->YYY=Y XXX=X->YYY=Y ...", where each "XXX=X" is the
    index (in reversed binary and decimal), and each "YYY=Y" is the
    value (ditto). Prints the lowest {m} bits of the argument and the
    lowest {n} bits of the value. The output is sanwiched
    between {pre} and {suf}. */

/* TC GATES. */
  
typedef uint16_t tcgate_t;
  /* A TC gate with some number {n <= n_max} wires.
    The encoding is defined in {tcgate_apply}. */

typedef int32_t tcgate_ct_t;
  /* A count or index of TC gates, including {-1} and 
    {4^n_max}. */

tcgate_ct_t tcgate_count(int n);
 /* Number of {n}-wire TC gates, including no-ops; namely {4^n}. */

tcgate_ct_t tcgate_proper_count(int n);
 /* Number of {n}-wire TC gates, excluding no-ops; namely, {4^n - 3^n}. */

word_t tcgate_apply(int n, tcgate_t g, word_t x);
  /* Applies the {n}-wire TC gate {g} to the word {x}. */

bool_t tcgate_is_no_op(int n, tcgate_t g);
  /* TRUE iff the TC gate {g} implements the identity function,
    namely, iff it has no inverters. */

void tcgate_print(FILE *wr, char *pre, int n, tcgate_t g, char *suf, bool_t masks);
  /* Writes the first {n} slices of the TC gate {g} to {wr},
    human-readable, between strings {pre} and {suf} If {masks} is
    true, prints also the binary masks of {g}. */

/* BIFUNS (INJECTIONS OF {B^N} TO {B^M}). */
  
typedef uint64_t bifun_t;
  /* An injective map {f} of {B^m} to {B^n}, for some {m <= m_max} and
    some {n <= n_max}. The encoding is defined by {bifun_pack} and
    {bifun_unpack} below. */

typedef int32_t bifun_ct_t;
  /* A count of bifuns, including {-1} and {bifun_count(n,m)} for all 
    supported combinations of {m} and {n}. */

bifun_ct_t bifun_count(int m, int n);
  /* Computes the number of bifuns from {B^m} to {B^n}, namely
    {(2^n_max)!/(2^n_max - 2^m_max)!}, or {0} if {n < m}. 
    Returns {-1} if the number is too large to be represented 
    in a {bifun_ct_t} variable. */

bifun_t bifun_ident(int m, int n);
  /* The canonical injection from {B^m} to {B^n}, that 
    merely pads each argument with {n-m} zeros at the 
    high end. */

word_t bifun_apply(int m, int n, bifun_t f, word_t x);
  /* Applies the bifun {f}, assumed to be from {B^m} to {B^n}, to the
    lowest {m} bits of word {x}, resulting in an {n}-bit word. */

bool_t bifun_is_proper(int m, int n, bifun_t f);
  /* TRUE iff the bifun {f} from {B^m} into {B^n} is actually a
    mapping from {B^m} to {B^m} --- that is, {f(x)} is in {B^m} for
    any {x} in {B^m}. The result is always TRUE if {m == n}. */

bifun_t bifun_tcgate_compose(int m, int n, bifun_t f, tcgate_t g);
  /* Returns the bifun from {B^m} to {B^n} that is equivalent to
    applying {f} then applying the gate {g}. */

typedef uint32_t bifun_hash_t;
  /* A bifun hash index, or size of bifun hash table. */

bifun_hash_t bifun_hash(bifun_t f, bifun_hash_t M);
  /* Maps the bifun {f} to number in {0..M-1},
    hopefully uniformly. */

bifun_hash_t bifun_hash_step(bifun_t f, bifun_hash_t M);
  /* Maps the bifun {f} to a number in {0..M-1} that is 
    relatively prime with {M}, hopefully in a uniform way. */

void bifun_packed_print(FILE *wr, char *pre, int m, int n, bifun_t f, char *suf);
  /* Writes the bifun {f} to {wr} in the format "#d.d....d".
    The whole is sanwiched between {pre} and {suf}. */

void bifun_print_h(FILE *wr, char *pre, int m, int n, bifun_t f, char *suf);
  /* Writes bifun {f} to {wr}, human-readable, in single-line format,
    namely the packd bifun as in {}bifun_packed_print, followed by " = " and the
    word table in horizontal format. The whole is sanwiched between
    {pre} and {suf}. */

void bifun_print_v(FILE *wr, char *head, int ind, int m, int n, bifun_t f, char *foot);
  /* Writes bifun {f} to {wr}, human-readable, one entry per line.
    Namely, writes the packed bifun (as per {bifun_packed_print}) in the
    first line, then the argument-value table in vertical format. Each
    line is indented by {ind} spaces. The whole table is sandwiched
    between lines containing {head} and {foot} (if not NULL). */

#endif
