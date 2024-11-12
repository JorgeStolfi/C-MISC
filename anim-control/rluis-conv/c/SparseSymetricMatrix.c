
#include <SparseSymetricMatrix.h>



/* Written by R.L.W.Liesenfeld, 1995 */

#include <Rd.h>
#include <Wr.h>
#include <Fmt.h>
#include <Lex.h>
#include <Util.h>
#include <FileFmt.h>
#include <FPut.h>
#include <fget_ NPut.h>
#include <nget_
#include <Thread.h>

T New(nat n)
  {
    { /* with*/ 
      m == NEW(T, n)
    ) {
      for (i = 0;  i < n; i++) {
        { /* with*/ mi == m[i] ) {
          mi.d = 0;
          mi.row = NEW(REF Row, 1);
          mi.row[0] = (MatrixElement){j = i, a = 1.0d0};
        };
      };
    };
  } /* New */;

double GetValue(READONLY T m, i, nat j)
  nat *?k, n;
  {
    { /* with*/ mi == m[i], row == mi.row^ ) {
      if ((j < i )) {
        k = 0;
        n = mi.d;
      } else {
        k = mi.d;
        n = row.nel;
      }
      while (k < n) && (row[k].j < j ) { k++;}
      if ((k < n) && (row[k].j == j )) { return row[k].a } else { return 0.0D0;
    };
  } /* GetValue */;

double GetDiagonalValue(READONLY T m, nat i)
  { 
    return m[i].row[m[i].d].a;
  } /* GetDiagonalValue */;

void SubMul(READONLY T m, READONLY Vector v, VAR Vector x)
  {
    for (i = 0;  i <= ((m.nel - 1)|?|MAX_m);  i++) {
      { /* with*/ r == m[i].row^ ) {
        for (k = 0;  k <= ((r.nel - 1)|?|MAX_r);  k++) {
          x[i] = x[i] - r[k].a * v[r[k].j];
        };
      };
    };
  } /* SubMul */;

void PutRow(VAR T m, nat i, READONLY Vector row)
  VAR k = 0;
  {
    affirm(row[i] != 0.0D0 , "??");
    for (j = 0;  j <= ((row.nel - 1)|?|MAX_row);  j++) {if ((row[j] != 0.0D0 )) { k++;} /*; } */;
    { /* with*/ crow == NEW(REF Row, k), cr == crow^ ) {
      k = 0;
      for (j = 0;  j <= ((row.nel - 1)|?|MAX_row);  j++) {
        if ((row[j] != 0.0D0)) || (j == i )) {
          cr[k].a = row[j];
          cr[k].j = j;
          if ((j == i )) { m[i].d = k; }
          k++;
        };
      }
      m[i].row = crow;
    };
  } /* PutRow */;

PROCEDURE Factor(READONLY T m): REF T == 
  double *?a, s;
      nat p;
  {
    { /* with*/ 
      n == m.nel,
      rrf == NEW(REF Vector, n), rf == rrf^,
      rmf == NEW(REF T, n), mf == rmf^
    ) {
      for (i = 0;  i < n; i++) {
        { /* with*/ 
          ri == m[i].row^,
          jd == m[i].d,
          nz == ri.nel
        ) {
          for (j = 0;  j < i; j++) {rf[j] = GetValue(mf, j, i); }

          s = 0.0D0;
          for (j = 0;  j < i; j++) {
            { /* with*/ rfj == rf[j] ) { s = s + GetDiagonalValue(mf, j)*rfj*rfj;
          }
          rf[i] = ri[jd].a - s;

          p = jd;
          for (j = i+1;  j < n; j++) {
            while (p < nz) && (ri[p].j < j ) { p++;}
            if ((p < nz) && (ri[p].j == j )) { a = ri[p].a } else { a = 0.0D0; }
            s = 0.0D0;
            for (k = 0;  k < i; k++) {
              s = s + GetDiagonalValue(mf, k)*rf[k]*GetValue(mf, k, j);
            }
            rf[j] = (a - s)/rf[i];
          }
          PutRow(mf, i, rf);
        };
      }
      return rmf;
    };
  } /* Factor */;

void Solve(READONLY T mf, READONLY Vector v, VAR Vector x, int i0)
  double *?s;
  {
    { /* with*/ 
      n == mf.nel 
    ) {
      for (i = 0;  i < i0; i++) {x[i] = 0.0D0; }

      /* Solves "LD*y == v", with "y" stored into x */
      for (i = i0;  i < n; i++) {
        { /* with*/ r == mf[i].row^, dj == mf[i].d ) {
          s = 0.0D0;
          for (j = 0;  j < dj; j++) {s = s + r[j].a * x[r[j].j]; }
          x[i] = (v[i] - s)/r[dj].a;
        };
      }

      /* Solves L'x == y */
      for (i = n-1;  i <= 0 BY -1;  i++) {
        { /* with*/ r == mf[i].row^, dj == mf[i].d ) {
          s = 0.0D0;
          for (j = dj+1;  j <= ((r.nel - 1)|?|MAX_r);  j++) {s = s + r[j].a*x[r[j].j]; }
          x[i] = x[i] - s;
        };
      };
    };
  } /* Solve */;

nat NonZeros(READONLY T m)
  VAR nz = 0;
  {
    for (i = 0;  i <= ((m.nel - 1)|?|MAX_m);  i++) {
      { /* with*/ r == m[i].row^ ) {
        nz = nz + r.nel;
        if ((r[m[i].d].a == 0.0d0 )) { nz--;
      };
    }
    return nz;
  } /* NonZeros */;

void Print(FILE *wr, READONLY T m)
  <* FATAL Alerted, Wr.Failure , "??");
  {
    { /* with*/ 
      n == m.nel,
      nr == ((double)n),
      n2 == 100.0D0/(nr*nr),
      nz == NonZeros(m),
      nzr == ((double)nz),
      d == Util.Digits(n)
    ) {
      fprintf(wr, "%s",  
        "Order: " & Fmt.Int(n) & "\n" &
        "Non-zeros: " & Fmt.Int(nz) & " == " & Fmt.Double(nzr*n2, prec = 1) & "%\n"
      );
      for (i = 0;  i < n; i++) {
        fprintf(wr, "%s",  Fmt.Pad(Fmt.Int(i), d) & ":");
        for (j = 0;  j <= m[i].d;  j++) {
          fprintf(wr, "%s",  
            " " & Fmt.Pad(Fmt.Int(m[i].row[j].j), d) & " == " & Fmt.Double(m[i].row[j].a)
          );
        }
        fprintf(wr, "\n");
      };
    }
    fflush(wr);
  } /* Print */;

void Remove(VAR T m, bool_vec sel)
  nat *?nk, newnz;
  {
    affirm(sel.nel == m.nel , "??");
    { /* with*/ 
      RM == ((nat.nel - 1)|?|MAX_CARDINAL)
    ) {
      for (i = 0;  i <= ((m.nel - 1)|?|MAX_m);  i++) {
        { /* with*/ 
          ri == m[i].row^,
          nz == ri.nel
        ) {
          if ((sel[i] )) {
            m[i].row = NEW(REF Row, 1);
            m[i].row[0].j = i;
            m[i].row[0].a = 1.0D0;
            m[i].d = 0;
          } else {
            newnz = nz;
            for (k = 0;  k < nz; k++) {
              if ((ri[k].j == i )) {
                ri[k].a = 1.0d0
              } else if ((sel[ri[k].j] )) { 
                ri[k].j = RM; newnz--;
              };
            }
            { /* with*/ newrow == NEW(REF Row, newnz), nri == newrow^ ) {
              nk = 0;
              for (k = 0;  k < nz; k++) {
                if ((ri[k].j != RM )) {
                  nri[nk] = ri[k];
                  if ((nri[nk].j == i )) { m[i].d = nk; }
                  nk++;
                };
              }
              m[i].row = newrow;
            };
          };
        };
      };
    };
  } /* Remove */;

CONST MatrixFileVersion == "96-12-04";

PROCEDURE Read(FILE *rd, double fp): REF T == 
  <* FATAL Alerted, Rd.Failure , "??");
  {
    FileFmt.ReadHeader(rd, "matrix", MatrixFileVersion);
    { /* with*/ 
      n == nget_Int(rd, "order"),
      fr == nget_Double(rd, "fingerprint"),
      matrix == NEW(REF T, n),
      m == matrix^
    ) {
      affirm(fp == 0.0d0) || (fp == fr , "??");
      for (i = 0;  i < n; i++) {
        Lex.Skip(rd);
        { /* with*/ ir == fget_Int(rd) ) { affirm(ir == i , "??"); }
        fget_Colon(rd);
        { /* with*/ nz == fget_Int(rd) ) { 
          fget_EOL(rd);
          m[i].row = NEW(REF Row, nz);
          m[i].d = ((nat.nel - 1)|?|MAX_CARDINAL);
          { /* with*/ ri == m[i].row^ ) {
            for (j = 0;  j < nz; j++) {
              Lex.Skip(rd);
              ri[j].j = fget_Int(rd);
              ri[j].a = fget_Double(rd);
              affirm(ri[j].j != i) || (ri[j].a != 0.0D0 , "??");
              if ((ri[j].j == i )) { m[i].d = j;
            };
          }
          fget_EOL(rd);
          affirm(m[i].d != ((nat.nel - 1)|?|MAX_CARDINAL) , "??");
        };
      }
      FileFmt.ReadFooter(rd, "matrix");
      return matrix;
    };
  } /* Read */;

void Write(FILE *wr, READONLY T m, double fp)
  <* FATAL Alerted, Wr.Failure , "??");
  {
    FileFmt.WriteHeader(wr, "matrix", MatrixFileVersion);
    fprintf(wr, "order = %d", m.nel);
    NPut.Double(wr, "fingerprint", fp);
    for (i = 0;  i <= ((m.nel - 1)|?|MAX_m);  i++) {
      { /* with*/ ri == m[i].row^, nz == ri.nel ) {
        fprintf(wr, "%d"nz); FPut.Colon(wr);
        for (j = 0;  j < nz; j++) {
          if ((j % 8 == 0 )) {
            fputc('\n', wr); FPut.Space(wr);
          }
          FPut.Space(wr);
          fprintf(wr, "%d"ri[j].j);
          FPut.Space(wr);
          FPut.Double(wr, ri[j].a);
        }
        fputc('\n', wr);
      };
    }
    FileFmt.WriteFooter(wr, "matrix");
    fflush(wr);
  } /* Write */;

{; } /* SparseSymetricMatrix */.
