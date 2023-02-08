/* SOComputeMatrices - matrices of scalar products from function bases. */
/* Last edited on 2008-05-24 12:24:33 by stolfi */

/* Computes the scalar product matrices {M[i,j] = <F[i],G[j]>}
  for two given bases {F[*]} and {G[*]} of functions.
  The dot product {<f,g>} is the integral over the root cell of 
  {f.eval(p)*geval(p)}, or {dot(f.grad(p),g.grad(p))}. */

#include <SOFunction.h> 
#include <SOBasic.h> 
#include <SOIntegral.h> 
#include <SOMatrix.h>
#include <SOBasisMatrices.h>
#include <SOParams.h>
#include <vec.h>
#include <affirm.h>
#include <nat.h>
#include <math.h>
#include <string.h>
#include <limits.h>
#include <values.h>
#include <stdio.h>
#include <stdlib.h>

typedef struct Options
  { char *basisNameF;    /* File of basis {F} */
    char *basisNameG;    /* File of basis {G} */
    char *outName;       /* Prefix for output matrix files. */
    bool_t ignoreTrees;    /* TRUE ignores the grid when computing the dot product. */
    nat gaussOrder;      /* Gaussian quadrature order for dot products. */
    bool_t eval;           /* TRUE to write the {<eval|eval>} matrix. */
    bool_t grad;           /* TRUE to write the {<grad|grad>} matrix. */
    double minVal;       /* Cleanup entries less than {minVal} times elem norm. */ 
    bool_t cholesky;       /* TRUE to write the Cholesky factor of {eval}. */
    bool_t checkAngles;    /* TRUE checks for close pairs in basis. */
    bool_t checkTriang;    /* TRUE checks the triangle inequality.  */
    bool_t verbose;        /* TRUE prints nonzero elements as they are computed. */
  } Options;

Options *GetOptions(int argn, char **argc);

Basis ReadBasis(char *name);
  /* Reads a function basis from the file {name} plus extension ".bas". */
  
void WriteMatrix(SOMatrix M, char *name); 
  /* Writes matrix {M} to the file {name} plus extension ".mat". */
  
void PrintSummary(SOMatrix M);
  /* Prints basic data about matrix {M}. */
  
void CheckAngles(SOMatrix M);
  /* Given the matrix {M} of dot products of a basis {B} against itself,
    checks whether {cos(B[i],B[j])} is close to 1.0 for
    any two distinct basis elements {B[i],B[j]}. */
     
void CheckTriangleInequality(SOMatrix M);
  /* Given the matrix {M} of dot products of a basis {B} against
    itself, checks whether it satisfies the triangle inequality
    {d(B[i],B[j]) \leq d(B[i],B[k])+d(B[k],B[j])}, where {d(f,g)} is
    computed from {M} by the formula 
    {d(f,g) = sqrt(<f|f> + <g|g> - 2*<f|g>)}. */
     
int main(int argn, char **argc)
  { Options *o = GetOptions(argn, argc);
    SOMatrix Meval, Mgrad;
    Basis F = ReadBasis(o->basisNameF);
    Basis G = (strcmp(o->basisNameG, o->basisNameF) == 0 ? F : ReadBasis(o->basisNameG));
    
    SOIntegral_GaussOrder = o->gaussOrder;

    if (o->eval)
      { Meval = SOBasisMatrices_BuildMatrixEval
          ( F, G, NULL, o->ignoreTrees, o->minVal, o->verbose );
        PrintSummary(Meval);
        WriteMatrix(Meval, txtcat(o->outName, "-ev"));
        if (o->cholesky)
          { SOMatrix Mevch = SOMatrix_Cholesky(Meval, o->minVal);
            PrintSummary(Mevch);
            WriteMatrix(Mevch, txtcat(o->outName, "-ev-ch"));
          }
        if (o->checkAngles) { CheckAngles(Meval); }
        if (o->checkTriang) { CheckTriangleInequality(Meval); }
        free(Meval.ents.e);
      }
    if (o->grad)
      { Mgrad = SOBasisMatrices_BuildMatrixGrad
          ( F, G, NULL, o->ignoreTrees, o->minVal, o->verbose );
        PrintSummary(Mgrad);
        WriteMatrix(Mgrad, txtcat(o->outName, "-gr"));
        free(Mgrad.ents.e);
      }
    return 0;
  }
  
Basis ReadBasis(char *name)
  { FILE *rd = open_read(txtcat(name, ".bas"), TRUE);
    Basis F = SOFunction_ReadBasis(rd);
    fclose(rd);
    return F;
  }
  
void WriteMatrix(SOMatrix M, char *name)
  { FILE *wr = open_write(txtcat(name, ".mat"), TRUE);
    SOMatrix_Write(wr, M);
    fclose(wr);
  }
  
void PrintSummary(SOMatrix M)
  { double tel = (double)(M.rows*M.cols);
    double nzel = (double)M.ents.ne;
    fprintf(stderr, "%5d rows\n", M.rows);
    fprintf(stderr, "%5d cols\n", M.cols);
    fprintf(stderr, "%8d entries\n", M.rows*M.cols);
    fprintf(stderr, "%8d nonzero entries\n", M.ents.ne);
    fprintf(stderr, "%8.6f occupancy fraction\n", nzel/tel);
    fprintf(stderr, "%8.6f mean entries per row\n", nzel/((double)M.rows));
    fprintf(stderr, "%8.6f mean entries per col\n", nzel/((double)M.cols));
  }

#define MAXCOS 0.9999
  /* Maximum value of {cos(B[i],B[j])} allowed in {checkAngles}. */

void CheckAngles(SOMatrix M)
  { MatEntry_vec_t Ments = M.ents;
    double_vec_t d = double_vec_new(M.rows);
    int i, k;
    
    affirm(M.rows = M.cols, "matrix is not square");
    
    for (i = 0; i < M.rows; i++) { d.e[i] = 0.0; }
    for (k = 0; k < Ments.ne; k++)
      { MatEntry *Mij = &(Ments.e[k]);
        if (Mij->row == Mij->col) 
          { d.e[Mij->row] = sqrt(Mij->va); }
      }
    for (i = 0; i < M.rows; i++) 
      { affirm(d.e[i] != 0.0, "diagonal element is zero"); }
    for (k = 0; k < Ments.ne; k++)
      { MatEntry *Mij = &(Ments.e[k]);
        nat i = Mij->row, j = Mij->col;
        if (i != j)
          { double c = Mij->va/d.e[i]/d.e[j];
            if (fabs(c) > MAXCOS)
              { fprintf(stderr, 
                  "basis elements %04d and %04d cos = %16.12f\n",
                  Mij->row, Mij->col, c);
              }
          }
      }
  }

void CheckTriangleInequality(SOMatrix M)
  { /* Dumb {O(N^3) log(N)} check - just to make sure. */
    int i, j, k;
    double eps = 1.0e-6;
    affirm(M.rows = M.cols, "matrix is not square");
    
    for (i = 0; i < M.rows; i++)
      { for (j = 0; j < M.cols; j++)
          { for (k = 0; k < M.cols; k++)
              { double mii = SOMatrix_GetValue(M, i,i); 
                double mij = SOMatrix_GetValue(M, i,j);
                double mjj = SOMatrix_GetValue(M, j,j);
                double mik = SOMatrix_GetValue(M, i,k);
                double mjk = SOMatrix_GetValue(M, j,k);
                double mkk = SOMatrix_GetValue(M, k,k);
                double dij2 = mii - 2*mij + mjj;
                double dik2 = mii - 2*mik + mkk;
                double dkj2 = mkk - 2*mjk + mjj;
                double hij = sqrt(fabs(mii)+2*fabs(mij)+fabs(mjj));
                double hik = sqrt(fabs(mii)+2*fabs(mik)+fabs(mkk));
                double hkj = sqrt(fabs(mkk)+2*fabs(mjk)+fabs(mjj));
                double dij = sqrt(fabs(dij2));
                double dik = sqrt(fabs(dik2));
                double dkj = sqrt(fabs(dkj2));
                
                affirm(dij2 >= -eps*hij, "distance dij2 negative");
                affirm(dik2 >= -eps*hik, "distance dik2 negative");
                affirm(dkj2 >= -eps*hkj, "distance dkj2 negative");
                affirm(dij <= dik + dkj, "triangle cond failure");
              }
          }
      }
  }

#define PPUSAGE SOParams_SetUsage

Options *GetOptions(int argn, char **argc)
  {
    Options *o = (Options *)notnull(malloc(sizeof(Options)), "no mem");
    SOParams_T *pp = SOParams_NewT(stderr, argn, argc);

    PPUSAGE(pp, "SOComputeMatrices \\\n");
    PPUSAGE(pp, "  { -basis NAME | -bases FNAME GNAME } \\\n");
    PPUSAGE(pp, "  [ -eval [ -cholesky ] ] [ -grad ] \\\n");
    PPUSAGE(pp, "  -gaussOrder NUM [ -minValue NUM ] [ -ignoreTrees ]\\\n");
    PPUSAGE(pp, "  [ -checkAngles ] [ -checkTriang ] \\\n");
    PPUSAGE(pp, "  -outName NAME \\\n");
    PPUSAGE(pp, "  [ -verbose ]\n");

    if (SOParams_KeywordPresent(pp, "-bases"))
      { o->basisNameF = SOParams_GetNext(pp);  
        o->basisNameG = SOParams_GetNext(pp); 
        o->checkAngles = FALSE;
        o->checkTriang = FALSE;
       }
    else if (SOParams_KeywordPresent(pp, "-basis"))
      { o->basisNameF = SOParams_GetNext(pp);  
        o->basisNameG = o->basisNameF; 
        o->checkAngles = SOParams_KeywordPresent(pp, "-checkAngles");
        o->checkTriang = SOParams_KeywordPresent(pp, "-checkTriang");
      }
    else
      { SOParams_Error(pp, "cannot find neither \"-basis\" nor \"-bases\""); }
        
    o->eval = SOParams_KeywordPresent(pp, "-eval");
    o->grad = SOParams_KeywordPresent(pp, "-grad");
    if (! (o->eval || o->grad))
      { SOParams_Error(pp, "cannot find neither \"-eval\" nor \"-grad\""); }
    if (o->eval)
      { o->cholesky = SOParams_KeywordPresent(pp, "-cholesky"); }
    else
      { o->cholesky = FALSE; }
   
    if (SOParams_KeywordPresent(pp, "-minVal"))
      { o->minVal = SOParams_GetNextDouble(pp, 0.0, DBL_MAX); }
    else 
      { o->minVal = 0.0; }
    
    SOParams_GetKeyword(pp, "-outName");                               
    o->outName = SOParams_GetNext(pp);  

    SOParams_GetKeyword(pp, "-gaussOrder");
    o->gaussOrder = SOParams_GetNextInt(pp, 1, INT_MAX);

    o->ignoreTrees = SOParams_KeywordPresent(pp, "-ignoreTrees");
    
    o->verbose = SOParams_KeywordPresent(pp, "-verbose");

    SOParams_Finish(pp);
    return o;
  }
  









