/* SPMatrixDecompTest - tests the Cholesky and SVD decompositions. */
/* Last edited on 2011-07-30 14:18:16 by stolfilocal */

#define PROG_NAME "SPMatrixDecompTest"

/* Tests {SPMatrix_SVD} and related routines. */

#include <SPMatrix.h>
#include <SPVector.h>
#include <SPOptions.h>
#include <SPSys.h>
#include <SPBasic.h>

#include <SPMatrixGSL.h>

#include <vec.h>
#include <affirm.h>
#include <nat.h>
#include <bool.h>
#include <js.h>

#include <gsl/gsl_vector_double.h>
#include <gsl/gsl_matrix_double.h>

#include <math.h>
#include <string.h>
#include <limits.h>
#include <values.h>
#include <stdio.h>
#include <stdlib.h>

#define PROG_USAGE \
  PROG_NAME " \\\n" \
  "  -size M N \\\n" \
  "  [ -cholesky ] [ -SVD ] \\\n" \
  "  [ -minVal NUM ] \\\n" \
  "  [ -verbose ]\n"

typedef struct Options
  { 
    int rows;              /* Number of rows to test. */
    int cols;              /* Number of columns to test. */
    double minVal;         /* Cleanup entries less than {minVal} times elem norm. */ 
    bool_t cholesky;       /* TRUE to write the Cholesky factor of {Mev}. */
    bool_t SVD;            /* TRUE to write the SVD factors {L,D,R} of {Mev}. */
    bool_t verbose;        /* TRUE prints nonzero elements as they are computed. */
  } Options;

Options *GetOptions(int argn, char **argc);

void CheckCholesky(SPMatrix *A, SPMatrix *AL, SPVector *b, SPVector *x);
  /* Verifies if {AL} is the left Cholesky factor of {A},
    and whether the Cholesky solver returns {x} as the solution
    of {A*x = b}. */

void CheckSVD(SPMatrix *A, SPMatrix *AL, SPVector *AD, SPMatrix *AR, SPVector *b, SPVector *x);
  /* Verifies if {AL,AD,AR} are the factors of the SVD decomposition
    {AL*AD*AR^T} of {A}, and whether the SVD solver returns {x} as the
    solution of {A*x = b}. */

void CheckOrtho(SPMatrix *M);
  /* Verifies that {M} is an orthonormal matrix. */
  
void ShowMatrix(SPMatrix M, char *title);
  /* Prints matrix {M} to stderr. */
  
void ShowVector(SPVector v, char *title);
  /* Prints vector {v} to stderr. */

void ShowSolutions(SPVector x, SPVector u);
  /* Prints the true solution {x}, the computed {u},
    and the error {u-x}, side by side, to stderr. */

void MakeSystem(int kind, int rows, int cols, SPMatrix *A, SPVector *b, SPVector *x);
  /* Returns the matrix {A}, the right-hand side {b},
    and the offcial solution {x} of the test system of kind {kind}
    and size {rows} by {cols}. */
  
#define NKINDS 4

int main(int argn, char **argc)
  { 
    Options *o = GetOptions(argn, argc);
    
    int kind;
    for (kind = 0; kind < NKINDS; kind++)
      {
        SPMatrix A; SPVector b, x;
        MakeSystem(kind, o->rows, o->cols, &A, &b, &x);
        ShowMatrix(A, "coeff matrix {A}");
        ShowVector(b, "right-hand side {b}");
        
        if (o->cholesky)
          { SPMatrix AL;
            SPMatrix_Cholesky(A, o->minVal, &AL);
            if (o->verbose)
              { ShowMatrix(AL, "AL"); }
            CheckCholesky(&A, &AL, &b, &x);
          }
        if (o->SVD)
          { SPMatrix AL, AR;
            SPVector AD;
            SPMatrix_SVD(A, o->minVal, &AL, &AD, &AR);
            if (o->verbose)
              { ShowMatrix(AL, "AL");
                ShowVector(AD, "AD");
                ShowMatrix(AR, "AR");
              }
            CheckSVD(&A, &AL, &AD, &AR, &b, &x);
          }
      }
    return 0;
  }
  
void MakeSystem(int kind, int rows, int cols, SPMatrix *A, SPVector *b, SPVector *x)
  { 
    /* Fill the matrix {A}: */
    gsl_matrix *GA = gsl_matrix_calloc(rows, cols);
    srandom(4615);
    int i;
    for (i = 0; i < rows; i++)
      { int j;
        for (j = 0; j < cols; j++)
          { switch(kind)
              { case 0:
                  /* Identity matrix: */
                  gsl_matrix_set(GA, i, j, (i == j ? 1.0 : 0.0)); 
                  break;
                case 1:
                  /* Diagonal matrix: */
                  gsl_matrix_set(GA, i, j, (i == j ? (double)(2*i+1) : 0.0)); 
                  break;
                case 2:
                  /* Positive definite matrix: */
                  gsl_matrix_set(GA, i, j, (i == j ? (double)(rows+cols) : 1.0)); 
                  break;
                case 3:
                  /* Random matrix: */
                  gsl_matrix_set(GA, i, j, drandom()); 
                  break;
                default:
                  affirm(FALSE, "bad matrix kind"); 
              }
          }
      }
    (*A) = SPMatrix_From_GSLMatrix(GA);
    gsl_matrix_free(GA);
        
    /* Fill the official solution vector {x}: */
    (*x) = double_vec_new(cols);
    int j;
    for (j = 0; j < cols; j++)
      { switch(kind)
          { case 0:
            case 1:
            case 2:
            case 3:
              /* Counting vector: */
              x->e[j] = (double)(j+1); 
              break;
            default:
              affirm(FALSE, "bad vector kind"); 
          }
      }

    /* Compute the right-hand side vector {b}: */
    (*b) = double_vec_new(rows);
    SPMatrix_MulCol(*A, *x, *b);
  }

void ShowMatrix(SPMatrix M, char *title)
  { 
    int nM = M.ents.ne;
    
    auto void header(bool_t dashes); 
      /* Prints header or dashed lines. */

    void header(bool_t dashes) 
      { int j;
        if (! dashes) 
          { fprintf(stderr, "%4s ", "row");
            for (j = 0; j < M.cols; j++) { fprintf(stderr, " %12d", j); }
            fprintf(stderr, "\n");
          }
        else
          { fprintf(stderr, "%4s ", "----");
            for (j = 0; j < M.cols; j++) { fprintf(stderr, " %12s", "------------"); }
            fprintf(stderr, "\n");;
          }
      }
    
    /* Print title: */
    fprintf(stderr, "=== %s ===\n", title);
    
    /* Print matrix: */
    header(FALSE);
    header(TRUE);
    int kM = 0;
    int i;
    MatEntry *Mk = &(M.ents.e[0]);
    for (i = 0; i < M.rows; i++)
      {
        fprintf(stderr, "%4d ", i);
        int j;
        for (j = 0; j < M.cols; j++) 
          { double va;
            if ((kM < nM) && (Mk->row == i) && (Mk->col == j))
              { va = Mk->va; kM++; Mk++; }
            else
              { va = 0; }
            fprintf(stderr, " %12.8g", va);
          }
        fprintf(stderr, "\n");
      }
    header(TRUE);
    fprintf(stderr, "\n");      
    affirm(kM == nM, "matrix elements out of order?");
    
    fprintf(stderr, "%5d rows\n", M.rows);
    fprintf(stderr, "%5d cols\n", M.cols);
    fprintf(stderr, "%8d entries\n", M.rows*M.cols);
    fprintf(stderr, "%8d nonzero entries\n", M.ents.ne);

    double tel = (double)(M.rows*M.cols);
    double nzel = (double)nM;
    fprintf(stderr, "%8.6f occupancy fraction\n", nzel/tel);
    fprintf(stderr, "%8.6f mean entries per row\n", nzel/((double)M.rows));
    fprintf(stderr, "%8.6f mean entries per col\n", nzel/((double)M.cols));
  }

void ShowVector(SPVector v, char *title)
  { 
    int nv = v.ne;

    auto void header(bool_t dashes); 
      /* Prints header or dashed lines. */

    void header(bool_t dashes) 
      { if (! dashes) 
          { fprintf(stderr, "%4s ", "pos");
            fprintf(stderr, " %12s", "value");
            fprintf(stderr, "\n");
          }
        else
          { fprintf(stderr, "%4s ", "----");
            fprintf(stderr, " %12s", "------------");
            fprintf(stderr, "\n");
          }
      }
    
    /* Print title and header: */
    fprintf(stderr, "=== %s ===\n\n", title);

    /* Print matrix: */
    header(FALSE);
    header(TRUE);
    int i;
    for (i = 0; i < nv; i++)
      { fprintf(stderr, "%4d ", i);
        fprintf(stderr, " %12.8g", v.e[i]);
        fprintf(stderr, "\n"); 
      }
    fprintf(stderr, "\n"); 
    header(TRUE);
  }
  
void CheckCholesky(SPMatrix *A, SPMatrix *AL, SPVector *b, SPVector *x)
  {
  }

void CheckSVD(SPMatrix *A, SPMatrix *AL, SPVector *AD, SPMatrix *AR, SPVector *b, SPVector *x)
  {
    SPVector u = double_vec_new(x->ne);
    SPVector t = double_vec_new(AD->ne);
    SPSys_SVDSolve(*AL, *AD, *AR, *b, t, u, TRUE);
    ShowSolutions(*x, u);
  }

void CheckOrtho(SPMatrix *M)
  {
  }

void ShowSolutions(SPVector x, SPVector u)
  {
    int nv = x.ne;
    affirm(u.ne == x.ne, "x/u size mismatch");

    auto void header(bool_t dashes); 
      /* Prints header or dashed lines. */

    void header(bool_t dashes) 
      { if (! dashes) 
          { fprintf(stderr, "%4s ", "pos");
            fprintf(stderr, " %12s", "true");
            fprintf(stderr, " %12s", "computed");
            fprintf(stderr, " %12s", "error");
            fprintf(stderr, "\n");
          }
        else
          { fprintf(stderr, "%4s ", "----");
            fprintf(stderr, " %12s", "------------");
            fprintf(stderr, " %12s", "------------");
            fprintf(stderr, " %12s", "------------");
            fprintf(stderr, "\n");
          }
      }
    
    /* Print matrix: */
    header(FALSE);
    header(TRUE);
    int i;
    for (i = 0; i < nv; i++)
      { fprintf(stderr, "%4d ", i);
        fprintf(stderr, " %12.8g", x.e[i]);
        fprintf(stderr, " %12.8g", u.e[i]);
        fprintf(stderr, " %12.8g", u.e[i] - x.e[i]);
        fprintf(stderr, "\n"); 
      }
    fprintf(stderr, "\n"); 
    header(TRUE);
  }

Options *GetOptions(int argn, char **argc)
  {
    Options *o = notnull(malloc(sizeof(Options)), "no mem");
    SPOptions_Parser_t *pp = SPOptions_NewParser(stderr, argn, argc);

    SPOptions_SetUsage(pp, PROG_USAGE);

    SPOptions_GetKeyword(pp, "-size");
    o->rows = SPOptions_GetNextInt(pp, 1, 5000);
    o->cols = SPOptions_GetNextInt(pp, 1, 5000);
    
    o->cholesky = SPOptions_TestKeyword(pp, "-cholesky");
    o->SVD = SPOptions_TestKeyword(pp, "-SVD");
   
    if (SPOptions_TestKeyword(pp, "-minVal"))
      { o->minVal = SPOptions_GetNextDouble(pp, 0.0, DBL_MAX); }
    else 
      { o->minVal = 0.0; }
    
    o->verbose = SPOptions_TestKeyword(pp, "-verbose");

    SPOptions_Finish(pp);
    return o;
  }
  









