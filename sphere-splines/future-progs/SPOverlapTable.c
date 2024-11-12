/* See SPOverlapTable.h */
/* Last edited on 2002-12-09 11:03:13 by anamaria */

#include <SPOverlapTable.h>
#include <NGet.h>
#include <FGet.h>
#include <stdio.h>
#include <FileFmt.h>
#include <basic.h>
#include <arrays.h>
#include <stdio.h>
#include <stdlib.h>
#include <affirm.h>

#define OverlapTable_FileFormat "2002-12-09"

void Write(FILE *wr, OverlapTable tb)
  { nat NS = tb.ne;
    int i;
    FileFmt.WriteHeader(wr, "overlap table", OverlapTable_FileFormat);
    fprintf(wr, "triangles = %d\n", NS);
    for (i = 0;  i < n; i++)
      { natArray row = tb.e[i];
        nat NT = row.ne;
        int j;
        fprintf(wr, "%d %d ", i, NT);
        for (j = 0; j < NT; j++)
          { fprintf(wr, " %d", row.e[j]); }
        fputc('\n', wr);
      }
    FileFmt.WriteFooter(wr, "overlap table");
    fflush(wr);
  }

OverlapTable Read(FILE *rd)
  { OverlapTable tb;
    nat NS;
    int i;
    FileFmt.ReadHeader(rd, "overlap table", OverlapTable_FileFormat);
    NS = NGet_Int(rd, "triangles"); FGet_EOL(rd);
    tb = natArrayArray_New(NS);
    for (i = 0; i < NS; i++)
      { nat NT;
        int ii, j;
        natArray row;
        ii = FGet_Int(rd);
        affirm(ii == i , "rows out of sequence");
        FGet_SkipSpaces(rd); 
        NT = FGet_Int(rd);
        row = natArray_New(NT);
        for (j = 0; j < NT; j++)
          { FGet_SkipSpaces(rd); row.e[j] = FGet_Int(rd); }
        FGet_EOL(rd);
        tb.e[i] = row;
      }
    FileFmt.ReadFooter(rd, "overlap table");
    return tb;
  }

/* Arrays of unsigned ints ({nat}s) */

natArrayArray natArrayArray_New(nat nel)
{
  natArrayArray arr;
  void *elv = malloc(nel * sizeof(natArray));
  arr.ne = nel;
  arr.e = (natArray *)notnull(elv, "out of mem for array" );
  return arr;
}

void natArrayArray_Expand(natArrayArray *arrp, nat i)
{
  if (i >= arrp->ne)
    { int nelt = 2 * (i > arrp->ne ? i : arrp->ne) + 1;
      void *elv = realloc(arrp->e, nelt * sizeof(natArray));
      arrp->e = (natArray *)notnull(elv, "out of mem for array expansion");
      arrp->ne = nelt;
    }
}

void natArrayArray_Trim(natArrayArray *arrp, nat nel)
{
  if (nel != arrp->ne)
    { void *elv = realloc(arrp->e, nel * sizeof(natArray));
      arrp->e = (natArray *)notnull(elv, "out of mem for array trimming");
      arrp->ne = nel;
    }
}
