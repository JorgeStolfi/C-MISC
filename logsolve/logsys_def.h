#ifndef logsys_def_H
#define logsys_def_H

/* Last edited on 2012-12-19 19:54:05 by stolfilocal */
/* Reveals the internal structure of {logsys.h} data types. */

#define _GNU_SOURCE
#include <stdint.h>

#include <bool.h>

#include <logsys.h>

#define ONE64 ((uint64_t)1)
#define ALL64 ((uint64_t)(-1LL))
#define ZERO64 ((uint64_t)0)

typedef struct logsys_arg_t  /* A use of a variable in an equation. */
  { logsys_va_t *va;     /* The argument variable. */ 
    logsys_eq_t *nxeq;   /* The next eq in a circular list of eqs that use that variable. */
    logsys_eq_t *preq;   /* The prev eq in a circular list of eqs that use that variable. */
  } logsys_arg_t;

struct logsys_eq_t
  { logsys_t *sys;           /* The system to which this equation belongs. */
    logsys_eq_id_t id;       /* A unique serial id of the equation node in {sys}. */
    uint8_t n;               /* Number of variables in equation (0 to 6). */
    logsys_op_t op;          /* The operation code. */
    /* Links for list of equations in system: */
    logsys_eq_t *nxeq;       /* The next eq in a circular list of all eqs in the system. */
    logsys_eq_t *preq;       /* The prev eq in a circular list of all eqs in the system. */
    /* Links for lists of equations that use each variable: */
    struct logsys_arg_t arg[logsys_n_MAX];  /* The variables used in this equation and their links. */ 
  };
  /* In a finished equation, the first {n} elements of {arg} must
    contain distinct non-null variables with non-null {nxeq,preq} lists.
    The remaining elements of {arg} must be all null. */

struct logsys_va_t
  { logsys_t *sys;           /* The system to which this variable belongs. */
    logsys_va_id_t id;       /* A unique serial id of the variable node in {sys}. */
    /* Links for list of variable users: */
    logsys_va_t *nxva;       /* The next variable in a circular list of all vars in the system. */
    logsys_va_t *prva;       /* The last variable in a circular list of all vars in the system. */
    logsys_eq_t *use;        /* Some eq that uses this variable, or NULL. */
  };

struct logsys_t
  { struct logsys_eq_t *eq;     /* Some equation in  the system. */
    struct logsys_va_t *va;     /* Some variable in the system. */
    logsys_va_id_t next_va_id;  /* Id of next variable to be added to this sys. */
    logsys_eq_id_t next_eq_id;  /* Id of next equation to be added to this sys. */
  };

#endif
