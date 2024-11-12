/* SPComputeMatrices - matrices of scalar products from function bases. */
/* Last edited on 2008-05-24 12:25:06 by stolfi */

#define PROG_NAME "SPComputeMatrices"

/* Given two bases {F[*]} and {G[*]} of spherical functions,
  Computes the scalar product matrices 

     { SE[i,j] := < F[j] | G[i] > }
     { SV[i,j] := < v ¤ SGrd(F[j]) | G[i] > }
     { SL[i,j] := < SLap(F[j]) | G[i] > = < SGrd(F[j]) | SGrd(G[i]) > }
     
  Here {v(s)} is the velocity vector field corresponding to rigid rotation around
  the Z axis with unit angular speed, and {¤} is scalar product.
  
  The notation { <f|g> } means the integral over the sphere of
  {f(p)*g(p)} for scalar fields, or {dot(f(p),g(p))} for vector
  fields; and {Sgrd(f)} is the spherical gradient of {f}. */

#include <SPFunction.h> 
#include <SPSpline.h> 
#include <SPIntegral.h> 
#include <SPMatrix.h>
#include <SPBasisMatrices.h>
#include <SPOptions.h>
#include <SPBasic.h>

#include <vec.h>
#include <affirm.h>
#include <nat.h>
#include <bool.h>

#include <math.h>
#include <string.h>
#include <limits.h>
#include <values.h>
#include <stdio.h>
#include <stdlib.h>

#define PROG_USAGE \
  PROG_NAME " \\\n" \
  "  { -basis NAME | -bases FNAME GNAME } \\\n" \
  "  [ -evMatrix [ -semi | -ortho ] [ -SVD ] [ -cholesky ] ] \\\n" \
  "  [ -vgMatrix [ -ortho ] ] \\\n" \
  "  [ -slMatrix [ -ortho ] ] \\\n" \
  "  [ -triName NAME ] \\\n" \
  "  -smpOrder NUM [ -minVal NUM ] [ -minSV NUM ] [ -ignoreTriangs ] \\\n" \
  "  [ -checkAngles ] [ -checkTriangles ] \\\n" \
  "  -outName NAME \\\n" \
  "  [ -verbose ]\n"

typedef struct Options_t
  { char *basisNameF;      /* File of basis {F} */
    char *basisNameG;      /* File of basis {G} */
    char *triName;         /* Name of triangulation file (minus ".tri"), or "". */
    char *outName;         /* Prefix for output matrix files. */
    bool_t ignoreTriangs;  /* TRUE ignores the triangulation in the dot product. */
    int smpOrder;          /* Triangle subdivision order for dot product integrals. */
    bool_t evMatrix;       /* TRUE to write the matrix {Mev = <f|g>}. */
    bool_t evSemi;         /* TRUE assumes that basis is semi-orthogonal. */
    bool_t evOrtho;        /* TRUE assumes that {Mev[i,j] = 0} if {i != j}. */
    bool_t vgMatrix;       /* TRUE to write the matrix {Mvg = <v¤SGrd(f)|g>}. */
    bool_t vgOrtho;        /* TRUE assumes that {Mvg[i,j] = 0} if {i != j}. */
    bool_t slMatrix;       /* TRUE to write the matrix {Msl = <SGrd(f)|SGrd(g)>}. */
    bool_t slOrtho;        /* TRUE assumes that {Msl[i,j] = 0} if {i != j}. */
    double minVal;         /* Cleanup entries less than {minVal} times elem norm. */ 
    double minSV;          /* Remove singular values less than {minSV}. */ 
    bool_t SVD;            /* TRUE to write the SVD factors {L,D,R} of {Mev}. */
    bool_t cholesky;       /* TRUE to write the Cholesky factor of {Mev}. */
    bool_t checkAngles;    /* TRUE checks for close pairs in basis. */
    bool_t checkTriangles; /* TRUE checks triangle inequality.  */
    bool_t verbose;        /* TRUE prints nonzero elements as they are computed. */
  } Options_t;

Options_t *GetOptions(int argn, char **argc);

Basis ReadBasis(char *name);
  /* Reads a function basis from the file {name} plus extension ".bas". */
  
Triangulation *GetTriangulation(char *name, int smpOrder, Basis F);
  /* Gets triangulation from file {triName}, or gets
    the common triangulation from {F} if {triName} is empty. */

void WriteMatrix(SPMatrix M, char *name); 
  /* Writes matrix {M} to the file {name} plus extension ".mat". */
  
void WriteVector(SPVector v, char *name);
  /* Writes vector {v} to the file {name} plus extension ".vec". */

void PrintSummary(SPMatrix M);
  /* Prints basic data about matrix {M}. */
  
void CheckAngles(SPMatrix M);
  /* Given the matrix {M} of symmetric dot products of a basis {B}
    against itself, checks whether {cos(B[i],B[j])} is close to 1.0
    for any two distinct basis elements {B[i],B[j]}. */
     
void CheckTriangles(SPMatrix M);
  /* Given the matrix {M} of symmetric dot products of a basis {B}
    against itself, checks whether it satisfies the triangle
    inequality {d(B[i],B[j]) \leq d(B[i],B[k])+d(B[k],B[j])}, where
    {d(f,g)} is computed from {M} by the formula {d(f,g) = sqrt(<f|f>
    + <g|g> - 2*<f|g>)}. */
     
int main(int argn, char **argc)
  { Options_t *o = GetOptions(argn, argc);
    SPIntegral_SetDefaultSamplingOrder(o->smpOrder);
    Basis F = ReadBasis(o->basisNameF);
    bool_t FeqG = (strcmp(o->basisNameG, o->basisNameF) == 0);
    Basis G = (FeqG ? F : ReadBasis(o->basisNameG));
    Triangulation *tri = GetTriangulation(o->triName, o->smpOrder, F);

    if (o->evMatrix)
      { SPMatrix Mev = SPBasisMatrices_BuildMatrixEval
          ( F, G, NULL, 
            tri, o->evSemi, o->evOrtho, o->ignoreTriangs, 
            o->minVal, o->verbose
	  );
        PrintSummary(Mev);
        WriteMatrix(Mev, txtcat(o->outName, "-ev"));
        if (o->SVD)
          { SPMatrix Mevsl, Mevsr;
            SPVector Mevsd;
            SPMatrix_SVD(Mev, o->minSV, &Mevsl, &Mevsd, &Mevsr);
            PrintSummary(Mevsl);
            WriteMatrix(Mevsl, txtcat(o->outName, "-ev-svdL"));
            PrintSummary(Mevsr);
            WriteVector(Mevsd, txtcat(o->outName, "-ev-svdD"));
            PrintSummary(Mevsr);
            WriteMatrix(Mevsr, txtcat(o->outName, "-ev-svdR"));
          }
        if (o->cholesky)
          { SPMatrix Mevch;
            SPMatrix_Cholesky(Mev, o->minVal, &Mevch);
            PrintSummary(Mevch);
            WriteMatrix(Mevch, txtcat(o->outName, "-ev-choL"));
          }
        if (strcmp(o->basisNameF, o->basisNameG) == 0)
          { if (o->checkAngles) { CheckAngles(Mev); }
            if (o->checkTriangles) { CheckTriangles(Mev); }
          }
        else
          { if (o->checkAngles) { fprintf(stderr, "\"-checkAngles\" ignored"); }
            if (o->checkTriangles) { fprintf(stderr, "\"-checkTriangles\" ignored"); }
          }
        free(Mev.ents.e);
      }
    if (o->vgMatrix)
      { SPMatrix Mvg = SPBasisMatrices_BuildMatrixVelSGrd
          ( F, G, NULL, o->vgOrtho, o->ignoreTriangs, o->minVal, o->verbose );
        PrintSummary(Mvg);
        WriteMatrix(Mvg, txtcat(o->outName, "-vg"));
        free(Mvg.ents.e);
      }
    if (o->slMatrix)
      { SPMatrix Msl = SPBasisMatrices_BuildMatrixSLap
          ( F, G, NULL, o->slOrtho, o->ignoreTriangs, o->minVal, o->verbose );
        PrintSummary(Msl);
        WriteMatrix(Msl, txtcat(o->outName, "-sl"));
        free(Msl.ents.e);
      }
    return 0;
  }
  
Basis ReadBasis(char *name)
  { FILE *rd = open_read(txtcat(name, ".bas"), TRUE);
    Basis F = SPFunction_ReadBasis(rd);
    fclose(rd);
    return F;
  }
  
Triangulation *GetTriangulation(char *name, int smpOrder, Basis F)
  { if ((name != NULL) && (*name != '\000'))
      { FILE *rd = open_read(txtcat(name, ".tri"), TRUE);
        Triangulation *tri = SPTriang_Read(rd, smpOrder);
        fclose(rd);
        return tri;
      }
    else
      { return SPSpline_BasisTriangulation(F); }
  }
  
void WriteMatrix(SPMatrix M, char *name)
  { FILE *wr = open_write(txtcat(name, ".mat"), TRUE);
    SPMatrix_Write(wr, M);
    fclose(wr);
  }
  
void WriteVector(SPVector v, char *name)
  { FILE *wr = open_write(txtcat(name, ".vec"), TRUE);
    SPVector_Write(wr, v);
    fclose(wr);
  }
  
void PrintSummary(SPMatrix M)
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

void CheckAngles(SPMatrix M)
  { MatEntry_vec_t Ments = M.ents;
    SPVector d = double_vec_new(M.rows);
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
        int i = Mij->row, j = Mij->col;
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

void CheckTriangles(SPMatrix M)
  { /* Dumb {O(N^3) log(N)} check - just to make sure. */
    int i, j, k;
    double eps = 1.0e-6;
    affirm(M.rows = M.cols, "matrix is not square");
    
    for (i = 0; i < M.rows; i++)
      { for (j = 0; j < M.cols; j++)
          { for (k = 0; k < M.cols; k++)
              { double mii = SPMatrix_GetValue(M, i,i); 
                double mij = SPMatrix_GetValue(M, i,j);
                double mjj = SPMatrix_GetValue(M, j,j);
                double mik = SPMatrix_GetValue(M, i,k);
                double mjk = SPMatrix_GetValue(M, j,k);
                double mkk = SPMatrix_GetValue(M, k,k);
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

Options_t *GetOptions(int argn, char **argc)
  {
    Options_t *o = notnull(malloc(sizeof(Options_t)), "no mem");
    SPOptions_Parser_t *pp = SPOptions_NewParser(stderr, argn, argc);

    SPOptions_SetUsage(pp, PROG_USAGE);

    if (SPOptions_TestKeyword(pp, "-bases"))
      { o->basisNameF = SPOptions_GetNext(pp);  
        o->basisNameG = SPOptions_GetNext(pp); 
        o->checkAngles = FALSE;
        o->checkTriangles = FALSE;
       }
    else if (SPOptions_TestKeyword(pp, "-basis"))
      { o->basisNameF = SPOptions_GetNext(pp);  
        o->basisNameG = o->basisNameF; 
      }
    else
      { SPOptions_Error(pp, "cannot find neither \"-basis\" nor \"-bases\""); }
        
    o->checkAngles = SPOptions_TestKeyword(pp, "-checkAngles");
    o->checkTriangles = SPOptions_TestKeyword(pp, "-checkTriangles");

    if (SPOptions_TestKeyword(pp, "-evMatrix"))
      { o->evMatrix = TRUE; 
        if (SPOptions_TestKeywordNext(pp, "-ortho"))
	  { o->evOrtho = TRUE; o->evSemi = FALSE; }
        else if (SPOptions_TestKeywordNext(pp, "-semi"))
	  { o->evOrtho = FALSE; o->evSemi = TRUE; }
        else
	  { o->evOrtho = FALSE; o->evSemi = FALSE; }
        o->cholesky = o->SVD = FALSE;
        while(TRUE)
          { if (SPOptions_TestKeywordNext(pp, "-cholesky"))
              { o->cholesky = TRUE; }
            else if (SPOptions_TestKeyword(pp, "-SVD"))
              { o->SVD = TRUE; }
            else
              { break; }
          }
      }
    else 
      { o->evMatrix = o->evOrtho = o->cholesky = o->SVD = FALSE; }

    if (SPOptions_TestKeyword(pp, "-vgMatrix"))
      { o->vgMatrix = TRUE; 
        o->vgOrtho = SPOptions_TestKeywordNext(pp, "-ortho");
      }
    else 
      { o->vgMatrix = o->vgOrtho = FALSE; }

    if (SPOptions_TestKeyword(pp, "-slMatrix"))
      { o->slMatrix = TRUE; 
        o->slOrtho = SPOptions_TestKeywordNext(pp, "-ortho");
      }
    else 
      { o->slMatrix = o->slOrtho = FALSE; }

    if (! (o->evMatrix || o->vgMatrix || o->slMatrix))
      { SPOptions_Error(pp, "no output matrix?"); }
   
    if (SPOptions_TestKeyword(pp, "-triName"))                               
      { o->triName = SPOptions_GetNext(pp); }
    else
      { o->triName = ""; }

    if (SPOptions_TestKeyword(pp, "-minVal"))
      { o->minVal = SPOptions_GetNextDouble(pp, 0.0, DBL_MAX); }
    else 
      { o->minVal = 0.0; }
    
    if (SPOptions_TestKeyword(pp, "-minSV"))
      { o->minSV = SPOptions_GetNextDouble(pp, 0.0, DBL_MAX); }
    else 
      { o->minSV = 0.0; }
    
    SPOptions_GetKeyword(pp, "-outName");                               
    o->outName = SPOptions_GetNext(pp);  

    SPOptions_GetKeyword(pp, "-smpOrder");
    o->smpOrder = SPOptions_GetNextInt(pp, 1, INT_MAX);

    o->ignoreTriangs = SPOptions_TestKeyword(pp, "-ignoreTriangs");
    
    o->verbose = SPOptions_TestKeyword(pp, "-verbose");

    SPOptions_Finish(pp);
    return o;
  }
  









