/* Last edited on 1997-11-24 00:51:47 by stolfi */
/* See the authorship and copyright notice at the end of this file. */
/* Changes:

*/

#include <jcltree.h>
#include <jclbasic.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <jclerror.h>

/* INTERNAL PROTOTYPES */

/* PROCEDURES */

void SortClusterTree(long *uch, long *vch, long R, long N, double **d, long *kap, long *kzp)
  { 
    long rr = abs(R);
    if (rr < N) 
      { (*kap) = rr; (*kzp) = rr; }
    else
      { long kua, kuz, kva, kvz;
        long ka, kz;
        double daa, daz, dza, dzz, dMin;
        SortClusterTree(uch, vch, uch[rr-N], N, d, &kua, &kuz);
        SortClusterTree(uch, vch, vch[rr-N], N, d, &kva, &kvz);
        daa = d[kua][kva]; daz = d[kua][kvz];
        dza = d[kuz][kva]; dzz = d[kuz][kvz];
        { double d1 = (daa < daz ? daa : daz);
          double d2 = (dza < dzz ? dza : dzz);
          dMin = (d1 < d2 ? d1 : d2);
        }
        /* Reverse direction of subtrees to decrease the gap: */
        if (dza == dMin)
          { /* No swapping necessary */ 
            ka = kua; kz = kvz; 
          }
        else if (daa == dMin)
          { /* Must reverse "u" subtree: */
            uch[rr-N] = -uch[rr-N]; 
            ka = kuz; kz = kvz;
          }
        else if (dzz == dMin)
          { /* Must reverse "v" subtree: */
            vch[rr-N] = -vch[rr-N];
            ka = kua; kz = kva; 
          }
        else if (daz == dMin)
          { /* Must reverse both subtrees: */
            uch[rr-N] = -uch[rr-N]; 
            vch[rr-N] = -vch[rr-N]; 
            ka = kuz; kz = kva; 
          }
        else 
          { Error("program error: dMin"); }
        if (R > 0)
          { (*kap) = ka; (*kzp) = kz; }
        else
          { (*kap) = kz; (*kap) = ka; }
      }
  }
  
void DumpClusterTree(long *uch, long *vch, long R, long N, long *it, long *ip)
  {
    long rr = abs(R);
    if (rr < N)
      { it[*ip] = rr; (*ip)++; }
    else if (R > 0)
      { DumpClusterTree(uch, vch, +uch[rr-N], N, it, ip);
        DumpClusterTree(uch, vch, +vch[rr-N], N, it, ip);
      }
    else
      { DumpClusterTree(uch, vch, -vch[rr-N], N, it, ip);
        DumpClusterTree(uch, vch, -uch[rr-N], N, it, ip);
      }
  }
  
