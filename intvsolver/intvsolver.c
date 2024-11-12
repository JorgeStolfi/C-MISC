/* Last edited on 2011-11-02 03:23:00 by stolfi */

#define PROG_NAME "intvsolver"
#define PROG_DESC "Interval-based solver for algebraic equations and inequations."
#define PROG_VERS "1.0"

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <stdint.h>
#include <string.h>

#include <bool.h>
#include <affirm.h>

#define INF INFINITY

typedef enum { 
    ives_op_ADD, /* Binary '+'. */
    ives_op_MUL, /* Binary '*'. */
    ives_op_NEG, /* Unary '-'. */
    ives_op_IDE, /* Unary identity. */
  } ives_op_t;

typedef struct ives_interval_t 
  { double end[2]; 
    bool_t open[2]; 
    bool_t neg; 
  } ives_interval_t;
/* A generalized interval of real numbers.

  Consider first the case when {neg} is false. If the endpoints {end[0]} and
  {end[1]} are the same value, then the interval is a singleton.  In this case, the value must be finite, and
  {open[0..1]} must be false.  If the endpoints are different, then
  the interval consists of all real numbers between those two
  endpoints. In this casee, each endpoint {end[i]} is included iff
  {open[i]} is false. If {open[0]} is true, {end[0]} may be
  {-INF}; if {open[1]} is true, {end[1]} may be {+INF}; otherwise
  both endpoints must be finite.

  If {neg} is true, the interval defined as above is complemented with
  respect to the whole reals. */

typedef struct ives_variable_t 
  { char *name; 
    int vix; 
    struct ives_variable_t *next; 
  } ives_variable_t;
/* An {ives_variable_t} object represents a variable in the problem; either
  one of the variables mentioned explicitly by the user, or an
  internal variable created by the formula parser. The {name}
  identifies the variable externally, e.g. for printout or debugging;
  the index {ix}, assigned serially, identifies it internally.
  The variables are chained in a list by the {next} field. */

typedef struct ives_equation_t 
  { int eix; 
    struct ives_variable_t *arg[3]; 
    ives_op_t op; 
    struct ives_equation_t *next; 
  } ives_equation_t;
/* An {ives_equation_t} object represents an equation or inequation of the
  form {arg[0] = arg[1] op arg[2]} (if {nex} is false) or {arg[0] !=
  arg[1] op arg[2]} (if {neg} is true).  If the operation {op} is unary 
  (or the identity), {arg[1]} is NULL. The index
  {ix} identifies the equation internally.  The equations are linked
  in a list by their {next} field. */

typedef struct ives_alist_t 
  { int vix; 
    ives_interval_t val; 
    struct ives_alist_t *rest; 
  } ives_alist_t;
/* An {ives_alist_t} is a state of the solver. Each element says that the value of 
  variable {vix} is known (or assumed) to be in the set {val}.    The same variable may appear two or
  more times in the a-list; only the first occurrence is relevant. */

/* PARSING FUNCTIONS

  The procedures {ives_parse_XXX} parse some construct from the input stream {rd}.
  If the construct is not present at the next position of {rd}, they return {NULL}
  without changing {rd}
  In case of error, they print a message to {stderr} and return the special pointer {BADP}.
  Otherwise they consume the construct from {rd}, leaving {rd} positioned
  just before the next character; and return the parsed construct.

  On successful parsing, the new variables, equations, and assignments created by the parsing 
  are prepended to the respective lists {*eqs,*ast,*var}.  Functions that parse a formula or sub-formula
  return as a result the variable which represents that formula's value. If the parsing fails, with or without error,
  the variables {*eqs,*ast,*var} are unchanged, and the returned result is {NULL}. The arguments
  {eqs,arg,vars} nust not be {NULL}. */

#define BADP ((void*)1)
/* Invalid pointer indicating a malformed syntactic element. */

typedef struct ives_parse_stream_t 
  { char *filename;     /* Name of underlying file,or "-" for stdin. */
    FILE *file;         /* Underlying file. */
    int line;           /* Current line. */
    int pos;            /* Current column of cursor (reset to 0 after a newline). */
    bool_t ungot;       /* True if an {ives_ungetchar} was executed with no subsequent {ives_getchar}. */ 
    int prevpos;        /* Right after a newline, this is the length of the prev line except the newline. */
  } ives_parse_stream_t;
/* A character stream that keeps track of the line and column counts. */

ives_variable_t *ives_parse_system(ives_parse_stream_t *rd, ives_equation_t **eqs, ives_alist_t **ast, ives_variable_t **var);
/* Parses a system of equations or inequation from {rd}.  The system
  must be enclosed in '{' and '}', and the equations must be separated
  by semicolons.  The equations {*eqs,*ast,*var} are normally NULL.  The result is
  either {NULL}, {BADP}, or the same as {*var}.  */

ives_variable_t *ives_parse_equation(ives_parse_stream_t *rd, ives_equation_t **eqs, ives_alist_t *ast, ives_variable_t **var);
/* Parses an equation {A rel B} where {A} and {B} are formulas and
  {rel} is a relation symbol: one of '=', '!=", '>', '<', '>=', or
  '<='. 

  Each equation is normally converted to the form {V=A-B}
  where {V} is a new variable, and the appropriate alist item for {V} is added to {ast}.
  The returned result is either {NULL}, {BADP}, or the mismatch variable {V}. */

ives_variable_t *ives_parse_formula(ives_parse_stream_t *rd, ives_equation_t **eqs, ives_alist_t *ast, ives_variable_t **var);
/* Parses an algebraic formula from {rd}.  The formula is broken down into a set of equations between variables,
  either explicitly named in the formula or created by the parser. */

ives_variable_t *ives_parse_term(ives_parse_stream_t *rd, ives_equation_t **eqs, ives_alist_t *ast, ives_variable_t **var);
/* Parses a term of a sum or difference, which may be one or more factors combined with '/' or '*'. */

ives_variable_t *ives_parse_factor(ives_parse_stream_t *rd, ives_equation_t **eqs, ives_alist_t *ast, ives_variable_t **var);
/* Parses a factor of a product or quotient, which may be an atom or a prefix operator, like '-' or function 
  name, applied to a factor. */

ives_variable_t *ives_parse_atom(ives_parse_stream_t *rd, ives_equation_t **eqs, ives_alist_t *ast, ives_variable_t **var); 
/* Parses an atomic formula, which many be a constant, a variable, or a formula in parentheses. */

int ives_getchar(ives_parse_stream_t *rd); 
/* Like {fgetc} for {rd}; updates the line and column count. */

void ives_ungetchar(ives_parse_stream_t *rd, int c); 
/* Like {ungetc} for {rd}. */

int ives_parse_get_non_blank(ives_parse_stream_t *rd);
/* Skip blank characters (ASCII SP, HT, LF, VT, FF, CR) from {rd} 
   and returns the first non-blank character found, possibly {EOF}. */
   
void ives_parse_skip_to_close_par(ives_parse_stream_t *rd);
/* Skips characters from {rd} until EOF, or just after the next unpaired parentheses or bracket;
   or just before the next semicolon, relation symbol, or brace. */

void ives_parse_skip_var(ives_parse_stream_t *rd);
/* Skips characters from {rd} until EOF, or just before a character that 
  is not letter, digit, or underscore. */

ives_variable_t *ives_add_variable(char *name, ives_variable_t **var);
/* Assumes that the list of currently defined variables begins at {*var} 
  and is linked by {.next}. Creates and returns a
  new variable record {v} and prepends it to the list {*var}. Assumes that
  the list is sorted by decreasing index {v.vix}. The index {v.vix} of the
  new variable will be the next free index (or 0 if {*var} is null). 
  The variable's name {v.name} will be the gienv {name} (which must not occur in the list)
  or "@{v.vix}" if {name} is NULL. */

ives_variable_t *ives_get_variable_by_name(char *name, ives_variable_t **var);
/* Assumes that the list of currently defined variables begins at {*var} 
  and is linked by {.next}. If the list contains a variable with the
  given {name} (which must not be NULL), returns that variable. Otherwise calls
  {ives_add_variable(name,var)}. */

ives_op_t ives_parse_unary_op(ives_parse_stream_t *rd);
/* If an unary operator is right upfront in {rd}, consumes it
  and returns the corresponding code.  Otherwise returns {ives_op_IDE}. */

#define ives_parse_error(rd, ...) \
  do \
  { fprintf(stderr, "%s:%d: column %d: ** ", rd->filename,rd->line, rd->pos); \
    fprintf(stderr, __VA_ARGS__); \
    fputc('\n',stderr); \
  } while(FALSE)

/* Prints an error message, reporting the position of {rd}. */

int ives_getchar(ives_parse_stream_t *rd)
{
  int c = fgetc(rd->file);
  rd->ungot = FALSE;
  if (c == '\n')
    { rd->line++; rd->prevpos = rd->pos; rd->pos = 0; }
  else if (c != EOF)
    { rd->pos++; rd->prevpos = -1; }
  return c;
}

void ives_ungetchar(ives_parse_stream_t *rd, int c) 
{
  demand(!rd->ungot, "cannot unget more than once in a row");
  demand(c != EOF, "cannot unget EOF");
  ungetc(c, rd->file);
  rd->ungot = TRUE;
  if (c == '\n')
    { demand((rd->pos == 0) && (rd->prevpos != -1), "invalid ungetc of newline");
      assert(rd->line > 0);
      rd->pos = rd->prevpos; rd->prevpos = -1; rd->line--;
    }
  else 
    { demand(rd->pos > 0, "invalid ungetc of non-newline");
      rd->pos--;
    }
}

int ives_parse_get_non_blank(ives_parse_stream_t *rd)
{
  do 
    { int c = ives_getchar(rd);
      if (c == EOF) 
        { return c; }
      else if ((c != '\011') && (c != '\012') && (c != '\014') && (c != '\015') && (c != 0xA0))
        { ives_ungetchar(rd, c);
          return c;
        }
    }
  while (TRUE);
}

void ives_parse_skip_to_close_par(ives_parse_stream_t *rd)
{
  int level = 1;
  do 
    { int c = ives_getchar(rd);
      if (c == EOF) 
        { return; }
      else if ((c == ';') || (c == '<') || (c == '>') || (c == '=') || (c == '!') || (c == '{') || (c == '}'))
        { ives_parse_error(rd, "missing %d close parentheses",level); 
          ives_ungetchar(rd, c);
          return;
        }
      else if ((c == '(') || (c == '['))
        { level++; }
      else if ((c == ')') || (c == ']'))
        { level--; 
          if (level == 0) { return; }
        }
    }
  while (TRUE);
}

void ives_parse_skip_var(ives_parse_stream_t *rd)
{
  do 
    { int c = ives_getchar(rd);
      if (c == EOF) 
        { return; }
      else if (((c < '0') || (c > '9')) && ((c < 'A') || (c > 'Z')) && ((c < 'a') || (c > 'z')) && (c != '_'))
        { ives_ungetchar(rd, c);
          return;
        }
    }
  while (TRUE);
}

ives_variable_t *ives_get_variable_by_name(char *name, ives_variable_t **var)
{
  demand(name != NULL, "name must not be null");
  ives_variable_t *tail = (*var);
  /* Check whether variable is new: */
  ives_variable_t *p = tail;
  while (p != NULL)
    { if (strcmp(p->name,name) == 0) { return p; }
      p = p->next;
    }
  /* New variable: */
  return ives_add_variable(name, var);
}

ives_variable_t *ives_add_variable(char *name, ives_variable_t **var)
{
  ives_variable_t *tail = (*var);
  int vix = (tail == NULL ? 0 : tail->vix + 1); /* Sequential variable index. */
  if (name == NULL) { asprintf(&name, "@%d", vix); }
  /* New variable: */
  ives_variable_t *vnew = notnull(malloc(sizeof(ives_variable_t)), "no mem");
  (*vnew) = (ives_variable_t)
    { .name = name,
      .vix = vix,
      .next = tail
    };
  (*var) = vnew;
  return vnew;
}

#define ives_parse_MAX_VAR_LENGTH 512
/* Max length of a variable name. Should be enough... */

ives_variable_t *ives_parse_atom(ives_parse_stream_t *rd, ives_equation_t **eqs, ives_alist_t *ast, ives_variable_t **var)
{
  int c = ives_parse_get_non_blank(rd);
  if (c == '(')
    { ives_variable_t *vres = ives_parse_formula(rd, eqs,ast,var);
      if (vres == NULL) 
        { ives_parse_error(rd, "expected formula after '('");
          ives_parse_skip_to_close_par(rd); 
          return BADP;
        }
      if (vres == BADP)
        { ives_parse_skip_to_close_par(rd); 
          return BADP;
       }
      c = ives_parse_get_non_blank(rd);
      if (c != ')') 
        { ives_parse_error(rd, "expected ')', found '%c'", c);
          ives_ungetchar(rd,c); 
          return BADP;
        }
      return vres;
    }
  else if (((c >= '0') && (c <= '9')) || (c == '.'))
    { uint64_t num = 0;
      uint64_t den = 1;
      bool_t dot = FALSE;
      do
        {
          if (c == '.')
            { if (dot) { ives_parse_error(rd, "number has two dots"); return BADP; }
              dot = TRUE;
	    }
	  else
	    { int dig = c - '0';
              if (num > (UINT64_MAX - dig)/10) { ives_parse_error(rd, "overflow in decimal constant"); return BADP; } 
              num = 10*num + dig;
	      if (dot)
		{ if (den > UINT64_MAX/10) { ives_parse_error(rd, "overflow in decimal denominator"); return BADP; } 
		  den = den*10;
		}
	    }
	  c = ives_getchar(rd);
	}
      while (((c >= '0') && (c <= '9')) || (c == '.'));
      ives_ungetchar(rd,c); 
      ives_variable_t *vres = ives_add_variable(NULL, var);
      ives_assign_value(vres,num,den,ast);
      return vres;
    }
  else if (((c >= 'A') && (c <= 'Z')) || ((c >= 'a') && (c <= 'z')))
    {
      char buf[ives_parse_MAX_VAR_LENGTH+1]; 
      int vlen = 1;
      buf[0] = c;
      buf[1] = 0;
      c = ives_getchar(rd);
      while (((c >= 'A') && (c <= 'Z')) || ((c >= 'a') && (c <= 'z')) || ((c >= '0') && (c <= '9')) || (c == '_'))
	{ if (vlen >= ives_parse_MAX_VAR_LENGTH) 
            { ives_parse_error(rd, "variable name '%s...' is too long", buf); 
              ives_parse_skip_var(rd); 
              return BADP;
            } 
          buf[vlen] = c; vlen++; buf[vlen] = 0;
          c = ives_getchar(rd);
	}
      ives_ungetchar(rd,c);
      return ives_get_variable_by_name(buf,var);
    }
  else
    { return BADP; }
}

ives_variable_t *ives_parse_factor(ives_parse_stream_t *rd, ives_equation_t **eqs, ives_alist_t *ast, ives_variable_t **var)
{
  ives_op_t op = ives_parse_unary_op(rd);
  ives_variable_t *varg = NULL;
  if (op != ives_op_IDE)
    { varg = ives_parse_factor(rd, eqs, ast, var);
      if (varg == NULL) { ives_parse_error(rd, "expected factor after unary op"); return BADP; } 
      if (varg == BADP) { return BADP; } 
      ives_variable_t *vres = ives_add_variable(NULL, var); 
      ives_add_equation(op, vres, NULL, varg, eqs);
      return vres;
    }
  else
    { varg = ives_parse_atom(rd, eqs, ast, var);
      return varg;
    }
}

ives_op_t ives_parse_unary_op(ives_parse_stream_t *rd);
{
  int c = ives_parse_get_non_blank(rd);
  if (c == '-')
    { return ives_op_NEG; }
  else
    { ives_ungetchar(rd,c); 
      return ives_op_IDE;
    }
}

ives_variable_t *ives_parse_term(ives_parse_stream_t *rd, ives_equation_t **eqs, ives_alist_t *ast, ives_variable_t **var)
{
  ives_variable_t *vterm = ives_parse_factor(rd, eqs, ast, var);
  if (vterm == BADP) { return BADP; }
  if (vterm == NULL) { return NULL; }
  ives_op_t op = ives_parse_mul_op(rd);
  while(op != ives_op_IDE)
    { ives_variable_t *varg = ives_parse_factor(rd, eqs, ast, var);
      if (varg == BADP) { return BADP; }
      if (varg == NULL) { ives_parse_error(rd, "expected factor after binary multiplicative op"); return BADP; } 
      ives_variable_t *vres = ives_add_variable(NULL, var);
      ives_add_equation(op, vres, vterm, varg, eqs);
      vterm = vres;
    }
  return vterm;
}

ives_variable_t *ives_parse_formula(ives_parse_stream_t *rd, ives_equation_t **eqs, ives_alist_t *ast, ives_variable_t **var)
{
  ives_variable_t *vform = ives_parse_term(rd, eqs, ast, var);
  if (vform == BADP) { return BADP; }
  if (vform == NULL) { return NULL; }
  ives_op_t op;
  while(ives_parse_add_op(rd, &op))
    { ives_variable_t *varg = ives_parse_term(rd, eqs, ast, var);
      if (varg == BADP) { return BADP; }
      if (varg == NULL) { ives_parse_error(rd, "expected term after binary additive op"); return BADP; } 
      ives_variable_t *vres = ives_add_variable(NULL, var);
      ives_add_equation(op, vres, vform, varg, eqs);
      vform = vres;
    }
  return vform;
}

ives_variable_t *ives_parse_equation(ives_parse_stream_t *rd, ives_equation_t **eqs, ives_alist_t *ast, ives_variable_t **var)
{
  ives_variable_t *va = ives_parse_formula(rd, eqs, ast, var);
  if (va == BADP) { return BADP; }
  if (va == NULL) { return NULL; }
  ives_rel_t rel;
  if(ives_parse_rel(rd, &rel))
    { ives_variable_t *vb = ives_parse_formula(rd, eqs, ast, var);
      if (varg == BADP) { return BADP; }
      if (varg == NULL) { ives_parse_error(rd, "expected formula after relational op"); return BADP; } 
      ives_variable_t vdif = ives_add_variable(NULL, var);
      ives_add_equation(ives_op_ADD, va, vdif, vb, eqs);
      ives_assign_value(vdif,0,1,ast);
      return vdif;
    }
  else
    { ives_parse_error(rd, "expected relational operator"); return BADP; } 
}

ives_variable_t *ives_parse_system(ives_parse_stream_t *rd, ives_equation_t **eqs, ives_alist_t **ast, ives_variable_t **var)
{
  char c = ives_getchar(rd);
  bool_t ok = TRUE;
  if (c == '{')
    { do
	{ ives_variable_t *vdif = ives_parse_equation(rd,eqs,ast,var);
          if (vdif != NULL) 
            { if (vdif == BADP) 
                { ives_parse_skip_rest_of_equation(rd); ok = FALSE; }
              c = ives_getchar(rd);
	    }
          if ((c != '}') && (c != ';')) 
            { ives_parse_error(rd, "expected ';' or '}', found '%c'", c);
              ives_parse_skip_rest_of_equation(rd); 
              ok = FALSE;
            }
	}
      while (c == ';');
      if (ok) 
         { assert((*var) != NULL); return (*var); }
      else
	{ return BADP; }
    }
  else
    { return NULL; }
}


