/* See SPMultiGrid.h */
/* Last edited on 2023-10-15 03:37:58 by stolfi */

#include <SPMultiGrid.h>
#include <SPTriang.h>
#include <SPOverlapTable.h>
#include <SPMatrix.h>
#include <SPPWFunction.h>
#include <SPTriang.h>
#include <r3x3.h>
#include <basic.h>
#include <stdio.h>
#include <stdlib.h>
#include <values.h>
#include <affirm.h>

FaceNum BruteForceLocate(Triangulation *tri, S2Point *p);
  /* Locates the triangle of {tri} that contains {p},
    by brute force search. */

FaceNum SPMultiGrid_ReLocate(MultiGrid *m, Level k, FaceNum fS, S2Point *p)
  { double minDist = INFINITY;
    FaceNum minFace;
    Triangulation *T = m.tri.e[k+1];
    natArray ovT = m.ovl.e[k].e[fn];
    int j;
    for (j = 0; j < ovT.ne; j++)
      { FaceNum fT = ovT.e[j];
        double d = TriDist(p, &(T->c2b.e[fT]));
        if (d <= 0.0) { return fT; }
        if (d < minDist) { minDist = d; minFace = fT; }
      }
    fprintf(stderr, "SPMultiGrid_ReLocate failure:");
    fprintf(stderr, " k = %d fS = %d", k, fS);
    fprintf(stderr, " minFace = %d minDist = %g\n", minFace, minDist);
    return minFace;
  }
  
FaceNum BruteForceLocate(Triangulation *tri, S2Point *p)
  { double minDist = INFINITY;
    FaceNum minFace;
    nat NT = tri->side.ne;
    FaceNum fT;
    for (fT = 0; fT < NT; fT++)
      { double d = TriDist(&p, &(tri->c2b[fT]))
        if (d <= 0.0) { return fT; }
        if (d < minDist) { minDist = d; minFace = fT; }
      }
    fprintf(stderr, "SPMultiGrid_BruteForceLocate failure:");
    fprintf(stderr, "%s",  " minDist == " & Fmt.Double(minDist));
    fprintf(stderr, "\n");
    return minFace;
  }
  
FaceNum Locate(MultiGrid *m, S2Point *p, Level k)
  { FaceNum fn = BruteForceLocate(m->tri.e[0], p);
    int kk;
    for (kk = 0; kk < k; kk++) { fn = ReLocate(m, kk, fn, p); }
    return fn;
  }

#define MultiGrid_FileFormat "2002-12-09"

void SPMultiGrid_Write(FILE *wr, MultiGrid *m)
  { affirm(FALSE , "SPMultiGrid_Write not implemented yet"); } 
    
MultiGrid Read(FILE *rd)
  { MultiGrid m;
    /*
      m.NL = NL;
      m.tri = TriangulationArray_New(NL+1);
      m.dim = natArray_New(NL+1);
      m.ovl = NEW(REFArray SPOverlapTable_T, NL);
      m.GGL = SPMatrixArray_New(NL+1);
      m.FG = SPMatrixArray_New(NL);
      m.SSL = SPMatrixArray_New(NL+1);
      for (k = 0;  k <= NL;  k++) {
        { kTag == Fmt.Pad(Fmt.Int(k), 2, '0'),
          k1Tag == Fmt.Pad(Fmt.Int(k+1), 2, '0')
          rd == open_read(txtcat4(gridName, "-", kTag, ".tri"), TRUE) ) {
          m.tri[k] = SPTriang_Read(rd); Rd.Close(rd);
          rd == open_read(txtcat4(gridName, "-", kTag, "-GGL.mat"), TRUE) ) {
          SPMatrix_Read(rd, m.GGL[k]); Rd.Close(rd);
          m.dim[k] = m.GGL[k].rows;
          rd == open_read(txtcat4(gridName, "-", kTag, "-", txtcat(opTag, "-SSL.mat")), TRUE) )
          SPMatrix_Read(rd, m.SSL[k]); Rd.Close(rd);
          assert(m.dim[k] == m.SSL[k].rows , {??});
          if ((k < NL ))
            { rd == open_read(txtcat4(gridName, "-", kTag, "-", txtcat(k1Tag, ".ovl")), TRUE) ) 
              m.ovl[k] = SPOverlapTable_Read(rd); Rd.Close(rd);
              assert(NUMBER(m.tri[k].side^) == NUMBER(m.ovl[k]^) , {??});
              rd == open_read(txtcat4(gridName, "-", kTag, "-", txtcat(k1Tag, "-FG.mat")), TRUE) ) {
              SPMatrix_Read(rd, m.FG[k]); Rd.Close(rd);
              assert(m.dim[k] == m.FG[k].rows , {??});
              assert(m.dim[k+1] == m.FG[k].cols , {??});
              ???or the other way around???;
            }
          }
        }
      }
    */
    affirm(FALSE , "SPMultiGrid_Read not implemented yet");
    return m;
  }
