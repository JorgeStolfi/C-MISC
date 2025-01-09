#define PROG_NAME "pc_decomp"
#define PROG_DESC "Principal component analysis on a set of data points."
#define PROG_VERS "2023-02-02"

#define pc_decomp_C_COPYRIGHT \
  "Copyright © 2023 by the State University of Campinas (UNICAMP)"
/* Last edited on 2024-12-01 00:23:34 by stolfi */
    
#define PROG_HELP \
  "  " PROG_NAME " \\\n" \
  "    -dataColumns {DATACOL}.. \\\n" \
  "    [ -weightColumn {WTCOL} ] \\\n" \
  "    [ -keyColumn {KEYCOL} ] \\\n" \
  "    -noise {NOISE} \\\n" \
  "    [ -centerFile {CENTERFILE} ] \\\n" \
  "    [ -covarFile {COVARFILE} ] \\\n" \
  "    [ -compFile {COMPFILE} ] \\\n" \
  "    [ -decompFile {DECOMPFILE} ] \\\n" \
  "    [ -format {FORMAT} ] \\\n" \
  "    [ -verbose {VERBOSE} ] \\\n" \
  "    " argparser_help_info_HELP " \\\n" \
  "    < {INFILE}"

#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  "  The program reads from {stdin} a list of {nd} /data points/ each with {nv} coordinates" \
  " or /variables/; that is, an array {D} on {nd} rows and {nv} columns, such" \
  " that D[id,jv]} is coordinate {jv \in 0..nv-1} of data point {id \in 0..nd-1}.\n" \
  "\n" \
  "  The program interprets that data as a cloud of points in the space {\RR^d}. It computes" \
  " the barycenter {d[0..nv-1]} of the cloud, and then finds the {ne} orthogonal" \
  " directions {E[0..ne-1]} of the cloud relative to {d} along which the cloud's" \
  " standard deviation {e[ie]} is successively maximum, in decreasing order of that" \
  " deviation; discarding directions whose {e[i]} is less than a given" \
  " threshold {tol}.\n" \
  "\n" \
  "  More precisely, the directions {E[0..ne-1]} are vectors with {nv} coordinates, unit norm," \
  " and pairwise orthogonal.  The vector {E[0]} is the direction of maximum cloud extent. The" \
  " vector {E[1]} is the maximum-extent direction that is perpendicular to {E[0]}; {E[2]} is" \
  " the maximum-extent direction perpendicular to {E[0]} and {E[1]}; and so on.\n" \
  "\n" \
  "  These directions are the /principal components/ of the data point set.\n" \
  "\n" \
  "  In other words, the program finds the best model of the cloud {D} that consists of" \
  " the barycenter {d} plus the sum of the vectors {E[0..ne-1]}, each {E[ie]} being multiplied" \
  " by {e[ie]} times an independent normal random variable.\n" \
  "\n" \
  "  The program also computes the projection {P[id]} of each data vector {D[id]-d} onto the" \
  " row space of {E}, that is a linear compbination of {E[0..ne-1]} with {ne} coefficients" \
  " {C[id]}; plus the residual {nv}-vector {R[id] = D[id]-d-P[íd]}. These results, as well as the" \
  " barycenter {d} and the component data {E,e}, may be written out to designated files.\n" \
  "\n" \
  "  Optionally the point reads also a weight {w[id]} for each data point (row of {D}). The" \
  " program then assumes that the data is actually {w[id]} copies of row {id}.  The weights" \
  " must be non-negative and may be fractional. Only their ratios matter.\n" \
  "\n" \
  "ALGEBRA\n" \
  "  In the notation of linear algebra, let {W} be the {nd} by {nd} diagonal matrix" \
  " whose element {W[id,id]} is the weight {w[id]}, normalized so that the sum of all" \
  " weights is 1.  Let {u} be the column {nd}-vector whose elements are all 1.\n" \
  "\n" \
  "  The" \
  " barycenter {d} then can be written as the 1 by {nv} (row) vector {u'*W*D}.  The covariance" \
  " matrix is then {A = V'*W*V} where {V} is the {nd} by {nv} matrix of displacements" \
  " from {d} to each {D[id]}, that is, {V = D-u*d}.\n" \
  "\n" \
  "  The pricipal components {E[0..ne-1]} are then the {ne} eigenvectors of the covariance" \
  " matrix {A} whose associated eigenvalues are at least {tol^2}. They are best" \
  " represented as an array {E} with {ne} rows and {nv} columns, whose rows are" \
  " orthonormal (meaning {E*E'} is the {ne} by {ne} identity matrix).\n" \
  "\n" \
  "  The /principal component analysis/ of the data points then consists of" \
  " computing the {nd} by {ne} coefficient matrix {C = E'*V}, the {nd} by {nv} projection" \
  " matrix {P = C*E}, and the {nd} by {nv} residual matrix {R = V - P}.\n" \
  "\n" \
  "INPUTS\n" \
  "  In the input file, blank lines and comments starting with \"#\" to the end" \
  " of the line ignored. Each data line normally has the format\n" \
  "\n" \
  "     {K[id]} {D[id,0]} ... {D[id,nv-1]}\n" \
  "\n" \
  "  or, if the data points are weighted,\n" \
  "\n" \
  "    {K[id]} {w[id]} {D[id,0]} ... {D[id,nv-1]} \n" \
  "\n" \
  "  where {K[id]} is a string with no embedded blanks that serves as an" \
  " identifier of the data point in the decomposition file.  However, the" \
  " order of these fields my be varied with the options \"-dataColumns\"," \
  " \"-weightColumn\", and \"-keyColumn\" (see below).  In general," \
  " apart from '#'-comments each non-blank line" \
  " must have the same number of non-blank fields separated by blanks.\n" \
  "\n" \
  "OUTPUTS\n" \
  "  The program may be asked to write separate files with the results. See" \
  " the options \"-centerFile\", \"-compFile\", and \"decompFile\" below.\n" \
  "OPTIONS\n" \
  "  -dataColumns {DATACOL}.. \n" \
  "    This mandatory argument specifies {nv} indices {col[0..nv-1]} that" \
  " indicate which columns of the input datafile are the data" \
  " values {D[i,jv]}.  Namely {D[id,jv]} will be field {col[jv]} of" \
  " line {id} of the file.  As per GNU/Linux tradition, columns are" \
  " numbered from 1, not 0.  Thus each {col[jv]} must be a number" \
  " in {1..nf}, where {nf} is the number of fileds in every input line.\n" \
  "\n" \
  "  -weightColumn {WTCOL}\n" \
  "    This optional argument specifies the index of the column in" \
  " the input datafile that contains the weight {w[id]}  of the data" \
  " point {D[id]}.  If given, it must be a number" \
  " in {1..nf}.  If omitted, the program assumes there is no weight" \
  " column in the file, and all data points have the same " \
  "weight ({w[id] = 1} for all {id}).\n" \
  "\n" \
  "  -keyColumn {KEYCOL}\n" \
  "    This optional argument specifies the index of the column in" \
  " the input datafile that contains the key {K[id]} identifying eack data" \
  " point {D[id]}.  If given, it must be a number" \
  " in {1..nf}.  If omitted, the program that {K[id]} is the string \"P{NNNNNN}\", where" \
  " {NNNNNN} is the point index {id} formatted as \"%06d\".\n" \
  "\n" \
  "  -noise {NOISE}\n" \
  "    This mandatory argument specifies the standard deviation of the" \
  " noise expected in each variable {D[i,jv]}.  The program will exclude" \
  " from the principal components {E[0..ne-1} any vector whose" \
  " deviation {e[ie]} is less than or equal to {NOISE} (which" \
  " should be non-negative).\n" \
  "\n" \
  "  -centerFile {CENTERFILE}\n" \
  "    This optional argument asks that the barycenter {d[0..nv-1]} of the" \
  " data points be written to the file named {CENTERFILE}. The {nv} coordinates" \
  " of this vector will be written, all on the same line, with the format" \
  " specified by the \"-format\" option.\n" \
  "\n" \
  "  -covarFile {COVARFILE}\n" \
  "    This optional argument asks that the covariance matrix {A[0..nv-1,0..nv-1]} of" \
  " the data points, relative to the barycenter, be written to the file" \
  " named {COVARFILE}. The elements will be written as {nv} lines with {nv} numbers" \
  " per line, in the format \"%24.15e\" \n" \
  "\n" \
  "  -compFile {COMPFILE}\n" \
  "    This optional argument asks that the matrix {E[0..ne-1,0..nv-1]} of principal" \
  " component vectors, and the associated deviations {e[0..ne-1]}, be written to the" \
  " file named {COMPFILE}.  Each row {E[ie]} and the associated deviation {e[ie]} will" \
  " be written, in that order, as a line of the file with {nv+1} numbers. The direction" \
  " vector coordinates will be printed with fixed format \"%+18.15f\", while the" \
  " deviation {e[ie]} will use the format specified by the \"-format\" option.\n" \
  "\n" \
  "  -decompFile {DECOMPFILE}\n" \
  "    This optional argument asks for an output file named {DECOMPFILE} with the" \
  " decomposition of the input data into the computed principal" \
  " components {E[0..ne-1]}.   Each line of the file have {1 + ne + nv + nv} fields: the" \
  " point key {K[id]}, the {ne} coefficients {C[id,0..ne-1]} (projections of {V[id]} along" \
  " each vector {E[ie]}), the {nv} components of the projection {P[id]} of {V[id]} on the" \
  " row space of {E}, and the {nv} components of the residual {R[id]} orthogonal to" \
  " that space. All these numbers will be in the format specified by the \"-format\" option.\n" \
  "\n" \
  "  -format {FORMAT}\n" \
  "    Specifies a format for printing point coordinates and other commensurate" \
  " quantities, including the the coordinates of the barycenter {d} in {CENTERFILE}, the" \
  " deviations {e[ie]} in {COMPFILE}, and the elements of {C[id]}, {P[id]}, and" \
  " {R[id]} in {DECOMPFILE}.  If omitted, defaults to \"-format %+12.5f\".\n" \
  "\n" \
  "  -verbose {VERBOSE}\n" \
  "    If present, this optional argumet asks the program to print various diagnorstic information.\n" \
  "\n" \
  "DOCUMENTATION OPTIONS\n" \
  argparser_help_info_HELP_INFO "\n" \
  "\n" \
  "SEE ALSO\n" \
  "  lineat_fit(JS).\n" \
  "\n" \
  "AUTHOR\n" \
  "  Created 2023-02-07 by Jorge Stolfi, IC-UNICAMP.\n" \
  "\n" \
  "MODIFICATION HISTORY\n" \
  "  2023-02-07 J.Stolfi: Created from {linear_fit.c}.\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " pc_decomp_C_COPYRIGHT ".\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS

#define stringify(x) strngf(x)
#define strngf(x) #x

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include <affirm.h>
#include <argparser.h>
#include <argparser_get_int_list.h>
#include <jsfile.h>
#include <vec.h>
#include <bool.h>
#include <rn.h>
#include <rmxn.h>
#include <jspca.h>
#include <fget.h>
#include <fget_data.h>

vec_typedef(fdtype_vec_t, fdype_vec, fget_data_type_t);

typedef struct jpc_options_t
  { int64_vec_t dataColumns;  /* Columns of input file with the data point coordinates. */
    int32_t keyColumn;        /* Column of input file with the point's key, or -1 if none. */
    int32_t weightColumn;     /* Column of input file with the point's weight, or -1 if none. */
    double noise;             /* Threshold for component magnitude. */
    char *centerFile;         /* Name of barycenter output file, or {NULL}. */
    char *covarFile;          /* Name of  output file, or {NULL}. */
    char *compFile;           /* Name of  output file, or {NULL}. */
    char *decompFile;         /* Name of  output file, or {NULL}. */
    char *format;             /* Format for point coordinates. */
    bool_t verbose;           /* True to print diagnostics. */
    /* Derived column data: */
    fdtype_vec_t field_type;  /* Type of each data field. */
  } jpc_options_t;
  /* Arguments from command line. */
  
#define jpc_MAX_VARS 200
  /* Max number of variables (point coordinates) */
  
#define jpc_MAX_FIELDS 1000
  /* Max number of data fields to parse or skip. */
     
void jpc_read_data
  ( FILE *rd, 
    jpc_options_t *o, 
    int32_t *nd_P,
    char ***K_P,
    double **w_P,
    int32_t *nv_P,
    double **D_P,
    bool_t verbose
  );
  /* Reads from {rd} the tuples {K[id]}, {w[id]}, {D[id,0],...D[id,nv-1]}, 
    for {id} in {0..nd-1}, into newly allocated vectors {K[0..nd-1]}, {w[0..nd-1]},
    and {D[0..nd*nv-1]}.  
    
    The array {D} conceptually has {nd} rows and {nv} columns. It is
    linearized by rows; element {D[id*nv+k]} is conceptually {D[id,jv]}
    for {jv} in {0..nv-1}. Returns the number of samples {nd} in {*NZp}
    and the arrays {ID,Z,w,D} {*IDp,*Zp,*Wp,*Xp}.
    
    The number of variables {nv} is determined by the option
    {o->dataColumns}. The number of data points {nd} is the number of
    lines in the file that are not blank apart from '#'-comments.
    
    The variables {K[id]}, {w[id]}, and {D[id,0..nv-1]} are read from
    the columnd of the file specified by {o->keyColumn},
    {o->weightColumn}, and {o->dataColumns}, respectively. These columns
    should be all distinct and between 1 and {nf} where {nf} is the
    number of non-blank items in each line (which should be the same for
    all lines).
    
    if {o->keyColumn} is {-1}, the key is set to "P{NNNNNN}" where
    {NNNNNN} is {id} in the format "%06d".
    
    If {o->weightColumn} is {-1}, all weights are set to 1.0.
    
    The integers {nv,nd} and the arrays {K,w,D} are returned in
    {*nv_P,*nd_P,*K_P,*w_P,*D_P}, respectively. */
   
void jpc_compute_barycenter(int32_t nd, double w[], int32_t nv, double D[], double d[], bool_t verbose);
  /* Computes the barycenter {d[0..nv-1]}, the average of the rows of
    {D[0..nd-1,0..nv-1]}, with weights {w{0..nd-1]}. */
 
void jpc_compute_principal_components
  ( int32_t nd, 
    double w[], 
    int32_t nv, 
    double D[], 
    double d[], 
    double A[],
    double noise,
    int32_t *ne,
    double E[],
    double e[],
    bool_t verbose
  );
  /* Given a point cloud {D[0..nd-1,0..nv-1]}, the point weights
    {w[0..nd-1]}, and the cloud barycenter {d[0..nv-1]}, computes the
    covariance matrix {A[0..nv-1,0..nv-1]} as well as the unit eigenvectors
    {E[0..nv-1,0..nv-1]} and the respective cloud extents
    (deviations) {e[0..nv-1]}.
    
    The eigenvectors (rows of {E}) will be sorted by decreasing extent {e[iv]}
    The procedure detetermines the number {ne} of eigenvectors
    whose extents {e[0..ne-1]} are at {noise} or more, and returns that
    number in {*ne_P). Rows {0..ne} of {E} are the principal components 
    of the cloud. */
 
void jpc_compute_decomposition
  ( int32_t nd, 
    int32_t nv, 
    double D[], 
    double d[], 
    int32_t ne, 
    double E[], 
    double **C_P, 
    double **P_P, 
    double **R_P,
    bool_t verbose
  );
  /* Given a point cloud {D[0..nd-1,0..nv-1]}, the cloud barycenter {d[0..nv-1]} and
    the principal directions {E[0..ne-1,0..nv-1]}, computes the coefficient matrix 
    {C[0..nd-1,0..ne-1]}, the point projections {P[0..nd-1,0..nv-1]} on the row space of {E},
    and the residuals {R[0..nd-1,0..nv-1]}. 
    
    The matrices {C[0..nd*ne-1]}, {P[0..nd*nv-1]},
    and {R[0..nd*nv-1]} are returned in {*C_P,*P_P,*R_P},
    respectively. */
 
void jpc_write_barycenter(char *fname, char *format, int32_t nv, double d[]);
void jpc_write_covariance_matrix(char *fname, int32_t nv, double A[]);
void jpc_write_principal_components(char *fname, char *format, int32_t ne, int32_t nv, double E[], double e[]);
void jpc_write_decomposition(char *fname, char *format, int32_t nd, int32_t ne, int32_t nv, double C[], double P[], double R[]);
  /* Write the given data tof the file {fname}. */

jpc_options_t *jpc_parse_options(int32_t argc, char **argv);
    
int32_t main (int32_t argc, char **argv)
  {
    jpc_options_t *o = jpc_parse_options(argc, argv);
    bool_t verbose = o->verbose;
    int32_t nv = -1; /* Number of coordinates per data point. */
    int32_t nd = -1; /* Number of data points. */
    
    char **K = NULL; /*  Data point ID key {K[0..nd-1]}. */
    double *w = NULL; /* Weights of data records are {w[0..nd-1]}. */
    double *D = NULL; /* Input independent variable values are {D[0..nd*nv-1]}. */
    jpc_read_data(stdin, o, &nd, &K, &w, &nv, &D, verbose);
    if (verbose) { fprintf(stderr, "read %d data samples {D[id,jv]}\n", nd); }

    /* Find point barycenter: */
    double *d = rn_alloc(nd); /* Point barycenter, {d[0..nv-1]}. */
    jpc_compute_barycenter(nd, w, nv, D, d, verbose);
    if (o->centerFile != NULL)
      { jpc_write_barycenter(o->centerFile, o->format, nv, d); }

    /* Allocate covariance matrix {A[0..nv-1,0..nv-1]}: */
    double *A = rmxn_alloc(nv,nv); 
    /* Allocate eigenvector matrix {E[0..nv-1,0..nv-1]}: */
    double *E = rmxn_alloc(nv,nv);
    /* Allocate extent vactor {e[0..nv-1]}: */
    double *e = rn_alloc(nv);
    int32_t ne = -1;  /* Number of principal components accepted. */
    jpc_compute_principal_components(nd, w, nv, D, d, A, o->noise, &ne, E, e, verbose);
    if (o->covarFile != NULL)
      { jpc_write_covariance_matrix(o->covarFile, nv, A); }
    if (o->compFile != NULL)
      { jpc_write_principal_components(o->compFile, o->format, ne, nv, E, e); }
      
    if (o->decompFile != NULL)
      { /* Compute the PC decomposition: */
        double *C = NULL; /* Coefficients of princ comps, {C[0..dn*ne-1]}. */
        double *P = NULL; /* Projection  {C[0..dn*ne-1]}. */
        double *R = NULL;
        jpc_compute_decomposition(nd, nv, D, d, ne, E, &C, &P, &R, verbose);
        jpc_write_decomposition(o->decompFile, o->format, nd, ne, nv, C, P, R);
      }
    
    return 0;
  }
 
void jpc_write_barycenter(char *fname, char *format, int32_t nv, double d[])
  { FILE *wr = open_write(fname, TRUE);
    for (uint32_t jv = 0;  jv < nv; jv++) 
      { fprintf(wr, format, d[jv]);
        if (jv > 0) { fputc(' ', wr); }
      }
    fputc('\n', wr);
    fclose(wr);
  }
      
void jpc_write_covariance_matrix(char *fname, int32_t nv, double A[])
  { FILE *wr = open_write(fname, TRUE);
    for (uint32_t kv = 0;  kv < nv; kv++) 
      { for (uint32_t jv = 0;  jv < nv; jv++) 
          { fprintf(wr, "%24.14e", A[kv*nv + jv]);
            if (jv > 0) { fputc(' ', wr); }
          }
        fputc('\n', wr);
      }
    fclose(wr);
  }
      
void jpc_write_principal_components(char *fname, char *format, int32_t ne, int32_t nv, double E[], double e[])
  { FILE *wr = open_write(fname, TRUE);
    for (uint32_t ie = 0;  ie < ne; ie++) 
      { for (uint32_t jv = 0;  jv < nv; jv++) 
          { fprintf(wr, "%+15.12f", E[ie*nv + jv]);
            if (jv > 0) { fputc(' ', wr); }
          }
        fputs("  ", wr);
        fprintf(wr, format, e[ie]);
        fputc('\n', wr);
      }
    fclose(wr);
  }

void jpc_write_decomposition(char *fname, char *format, int32_t nd, int32_t nv, int32_t ne, double C[], double P[], double R[])
  { FILE *wr = open_write(fname, TRUE);
    for (uint32_t id = 0;  id < nd; id++) 
      { for (uint32_t ie = 0;  ie < ne; ie++) 
          { fprintf(wr, format, C[id*ne + ie]);
            if (ie > 0) { fputc(' ', wr); }
          }
        fputs("  ", wr);
        for (uint32_t jv = 0;  jv < nv; jv++) 
          { fprintf(wr, format, P[id*nv + jv]);
            if (jv > 0) { fputc(' ', wr); }
          }
        fputs("  ", wr);
        for (uint32_t jv = 0;  jv < nv; jv++) 
          { fprintf(wr, format, R[id*nv + jv]);
            if (jv > 0) { fputc(' ', wr); }
          }
        fputc('\n', wr);
      }
    fclose(wr);
  }

void jpc_read_data
  ( FILE *rd, 
    jpc_options_t *o, 
    int32_t *nd_P,
    char ***K_P,
    double **w_P,
    int32_t *nv_P,
    double **D_P,
    bool_t verbose
  )
  { 
    /* Get the number {nv} and columns  of variables (point coordinates): */
    int32_t nv = o->dataColumns.ne;
    assert(nv > 0); 
    
    /* Get number {nf} of data fields to read/skip, and type table: */
    /* Note that there may be no weight of key. */
    int32_t nf = o->field_type.ne; 
    assert(nf >= nv); 
    fget_data_type_t *type = o->field_type.e;
    
    /* Allocate data arrays: */
    int32_t nd0 = 2048; /* Initial allocation; may be expanded. */
    string_vec_t K = string_vec_new(nd0);
    double_vec_t w = double_vec_new(nd0);
    double_vec_t D = double_vec_new(nv*nd0);

    /* Arrays {alf,num} for {fget_data_fields} results: */
    /* Note that these arrays are indexed {0..nf-1} */
    /* because {fget_data_line} numbers the columns from 0. */
    char *alf[nf];
    double num[nf];

    /* Read the data lines: */
    int32_t nd = 0; /* Number of data lines seen. */
    while (TRUE)
      { bool_t ok = fget_data_fields(rd, '#', nf, type, alf, num);
        if (! ok) { break; }
        int32_t id = nd; /* Row index of {D}, index of {K,w}. */

        /* Expand arrays if needed: */
        string_vec_expand(&K, id);
        double_vec_expand(&w, id);
        double_vec_expand(&D, id*nv + nv-1);
        
        /* Get the data point: */
        double *Di = &(D.e[id*nv]);
        for (uint32_t jv = 0;  jv < nv; jv++) 
          { int32_t kf = (int32_t)(o->dataColumns.e[jv] - 1);
            double Dij = num[kf];
            demand((! isnan(Dij)) && isfinite(Dij), "invalid data point coordinate");
            Di[jv] = Dij;
          }
        /* Get or fake the weight: */
        double wi;
        if (o->weightColumn > 0) 
          { int32_t kf = (int32_t)(o->weightColumn - 1);
            wi = num[kf];
            demand((! isnan(wi)) && isfinite(wi) && (wi >= 0), "invalid weight value");
          }
        else
          { wi = 1.0; }
        w.e[id] = wi;
        
        /* Get or fake the key: */
        char *Ki;
        if (o->keyColumn > 0) 
          { int32_t kf = (int32_t)(o->keyColumn - 1);
            char *Ki = alf[kf];
            assert(Ki != NULL);
          }
        else
          { char *Ki = NULL;
            char *Ki = jsprintf("P%06d", id);
          }
        K.e[id] = Ki;

        /* Skip rest of line: */
        fget_skip_to_eol(rd);
      }
    string_vec_trim(&K,nd);
    double_vec_trim(&w,nd);
    double_vec_trim(&D,nd*nv);
      
    (*nd_P) = nd;
    (*K_P) = K.e;
    (*w_P) = w.e;
    (*nd_P) = nd;
    (*D_P) = D.e;
   }

void jpc_compute_barycenter(int32_t nd, double w[], int32_t nv, double D[], double d[], bool_t verbose)
  {
    double d[nv];
    double sum_w = 0;
    for (uint32_t kv = 0;  kv < nv; kv++)
      { double sum_wD = 1.0e-308;
        for (int32_t id = 0; id < nd; id ++) 
          { double wi = w[id];
            double Dik = D[id*nv + kv];
            if (kv == 0) { sum_w += wi; }
            sum_wD += wi*Dik;
          }
        assert(sum_w > 0);
        d[kv] = sum_wD/sum_w;
      }
  }
 
void jpc_compute_principal_components
  ( int32_t nd, 
    double w[], 
    int32_t nv, 
    double D[], 
    double d[], 
    double A[],
    double noise,
    int32_t *ne_P,
    double E[],
    double e[],
    bool_t verbose
  )
  {
    if (
    jspca_compute_components(nd, nv, D, d, w, E, e, verbose);
    
    /* Determine the number {ne} of significant components: */
    int32_t ne = nv;
    while ((ne > 0) && (e[ne-1] < noise)) { ne--; }
    
    if (A_P != NULL) (*A_P) = A;
    (*ne_P) = ne;
    (*E_P) = E;
    if (A_P != NULL) (*A_P) = A;
    
  
  }
 
void jpc_compute_decomposition
  ( int32_t nd, 
    int32_t nv, 
    double D[], 
    double d[], 
    int32_t ne, 
    double E[], 
    double **C_P, 
    double **P_P, 
    double **R_P,
    bool_t verbose
  )
  {
  }


void jpc_build_model
  ( int32_t nd, 
    double Z[], 
    double w[], 
    int32_t nv, 
    double D[], 
    bool_t unitTerm, 
    int32_t ne, 
    double E[], 
    bool_t verbose
  )
  {
    double Xave[nv];     /* If {unitTerm}, average value of each {D[id]}. */
    double Zave = NAN;   /* If {unitTerm}, average value of {Z}. */
    if (unitTerm)
      { demand(ne == nv + 1, "{ne} should be {nv+1}");
        if (verbose) { fprintf(stderr, "removing variable averages...\n"); }
        /* Compute the weighted sums of {D[k],Z} in {Xave[k],Zave}: */
        for (uint32_t k = 0;  k < nv; k++) { Xave[k] = 0; }
        Zave = 0;
        double sumW = 0;
        for (uint32_t id = 0;  id < nd; id++)
          { double Wi = w[id];
            demand(Wi >= 0, "invalid weight {w[id]}");
            double *Xi = &(D[id*nv]); 
            for (uint32_t k = 0;  k < nv; k++) 
              { Xave[k] += Wi*Xi[k]; }
            Zave += Wi*Z[id];
            sumW += Wi;
          }
        /* Convert weighted sums to weighted averages: */
        if (sumW > 0)
          { for (uint32_t k = 0;  k < nv; k++) { Xave[k] /= sumW; }
            Zave /= sumW;
          }
        /* Subtract averages from all variables: */
        for (uint32_t id = 0;  id < nd; id++)
          { double *Xi = &(D[id*nv]); 
            for (uint32_t k = 0;  k < nv; k++) { Xi[k] -= Xave[k]; }
            Z[id] -= Zave;
          }
      }
    else
      { demand(ne == nv, "{ne} should be {nv}"); }

    /* Build the least squares system: */
    if (verbose) { fprintf(stderr, "building linear system matrices...\n"); }
    double *A = rmxn_alloc(nv, nv);
    double B[nv];
    rn_zero(nv, B);
    rmxn_zero(nv, nv, A);
    
    /* Accumulate the matrix and vector of the least squares system: */
    for (uint32_t id = 0;  id < nd; id++)
      { double *Xi = &(D[id*nv]); 
        jpc_accum_system(Z[id], w[id], nv, Xi, A, B);
      }
    if (verbose) { rmxn_gen_print2(stderr, nv,  nv, A,  1, B,  "%+18.9f", "  ","\n  ","\n", "[ "," "," ]", "  "); }
    
    if (verbose) { fprintf(stderr, "solving system...\n"); }
    double tiny = 1.0e-8;
    uint32_t rank;
    gausol_solve(nv, nv, A, 1, B, E, TRUE,TRUE, tiny, NULL, &rank);
    if (rank < nv) { fprintf(stderr, "!! warning: system rank = %d\n", rank); }
    if (verbose) { rmxn_gen_print(stderr, nv, 1, E, "%+18.9f", "  ","\n  ","\n", "[ "," "," ]"); }
    
    if (unitTerm)
      { if (verbose) { fprintf(stderr, "adding variable averages to indep term...\n"); }
        double Cunit = Zave;
        for (uint32_t k = 0;  k < nv; k++) { Cunit -= E[k]*Xave[k]; }
        E[ne-1] = Cunit;
      }
    
    free(A);
  }
  
void jpc_apply_model(int32_t nd, int32_t nv, double D[], bool_t unitTerm, int32_t ne, double E[], double Y[])
  { 
    assert(ne == nv + (unitTerm ? 1 : 0));
    for (uint32_t id = 0;  id < nd; id++)
      { double *Xi = &(D[id*nv]); 
        double sum = 0;
        for (uint32_t k = 0;  k < nv; k++) { sum += E[k]*Xi[k]; }
        if (unitTerm) { sum += E[nv]; }
        Y[id] = sum;
      }
  }

void jpc_residual_stats(int32_t nd, double Z[], double w[], double Y[], double *avgP, double *devP)
  {
    double sumWD = 0;
    double sumW = 0;
    for (uint32_t id = 0;  id < nd; id++) { sumWD += w[id]*(Z[id] - Y[id]); sumW += w[id]; }
    double avg = sumWD/sumW;
    demand(sumW > 0, "total weight is zero");
    double sumWD2 = 0;
    for (uint32_t id = 0;  id < nd; id++) { double Di = (Z[id] - Y[id]) - avg;  sumWD2 += w[id]*Di*Di; }
    double dev = sqrt(sumWD2/sumW);
    (*avgP) = avg;
    (*devP) = dev;
  }

void jpc_write_model(char *fname, int32_t nv, char *tName[], bool_t unitTerm, int32_t ne, double E[], double avg, double dev)
  {
    assert(ne == nv + (unitTerm ? 1 : 0));
    FILE *wr = open_write(fname, TRUE);
    fprintf(wr, "%d\n", nv);
    for (uint32_t k = 0;  k < nv; k++) 
      { fprintf(wr, "%3d %+24.16e", k, E[k]); 
        if (tName != NULL)
          { fprintf(wr, " # %s", tName[k]); }
        fprintf(wr, "\n");
      }
    if (unitTerm) { fprintf(wr, "%3d %+24.16e # 1\n", nv, E[nv]); }
    fprintf(wr, "    %+24.16e\n", avg);
    fprintf(wr, "    %+24.16e\n", dev);
    fclose(wr);
  }

void jpc_print_model(FILE *wr, int32_t nv, char *tName[], bool_t unitTerm, int32_t ne, double E[], double avg, double dev)
  {
    assert(ne == nv + (unitTerm ? 1 : 0));

    /* Index sort of {E[0..nv-1]} by absolute value of coefficient (leave {E[nv]} at end): */
    int32_t ix[ne];
    for (uint32_t k = 0;  k < ne; k++) { ix[k] = k; }
    for (uint32_t k = 0;  k < nv; k++) 
      { /* Largest coeffs are {E[ix[id]]} for {id} in {0..k-1}. */
        /* Find largest unsorted coeff: */
        int32_t imax = k;
        for (int32_t id = k+1; id < nv; id++)
          { if (fabs(E[ix[id]]) > fabs(E[ix[imax]])) { imax = id; } }
        if (imax != k)
          { int32_t t = ix[k]; ix[k] = ix[imax]; ix[imax] = t; }
      }
  
    fprintf(wr, "fitted model:\n");
    for (uint32_t k = 0;  k < ne; k++) 
      { int32_t id = ix[k];
        fprintf(wr, "  %+14.9f", E[id]); 
        if (id < nv)
          { fprintf(wr, " * D[%2d]", id);
            if (tName != NULL) { fprintf(wr, " # %s", tName[id]); }
          }
        else
          { fprintf(wr, "%8s", "");
            fprintf(wr, " # 1");
          }
        fprintf(wr, "\n"); 
      }
    fprintf(wr, "  %+14.9f         # avg residual\n", avg);
    fprintf(wr, "  %+14.9f * RND() # dev residual\n", dev);
  }
    
void jpc_accum_system(double Zi, double Wi, int32_t nv, double Xi[], double *A, double *B)
  { 
    for (uint32_t k = 0;  k < nv; k++)
      { B[k] += Wi*Zi*Xi[k];
        for (uint32_t jv = 0;  jv < nv; jv++) { A[k*nv + jv] += Wi*Xi[k]*Xi[jv]; }
      }
  }
  
void jpc_write_data(FILE *wr, int32_t nd, char *K[], double Z[], double Y[], char *fmt)
  {
    for (uint32_t id = 0;  id < nd; id++)
      { fprintf(wr, "%s ", K[id]);
        fprintf(wr, fmt, Z[id]);
        fprintf(wr, " ");
        fprintf(wr, fmt, Y[id]);
        fprintf(wr, "\n");
      }
    fflush(wr);
  }

#define jpc_terms_MAX 50
  /* Maximum number of terms (group averages) in model. */

jpc_options_t *jpc_parse_options(int32_t argc, char **argv)
  {
    /* Initialize argument parser: */
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);

    /* Allocate the command line argument record: */
    jpc_options_t *o = notnull(malloc(sizeof(jpc_options_t)), "no mem");

    /* Parse keyword parameters: */

    /* Get indices of data file colums where the variables are: */
    uint32_t nv_max = jpc_MAX_VARS;
    uint32_t nf_max = jpc_MAX_FIELDS;
    o->dataColumns = argparser_get_int_list(pp, 1000, nv_max, "-dataColumns", 1, (int32_t)nf_max);
    uint32_t nv = o->dataColumns.ne;  /* Number of variables. */
    if (nv == 0) { argparser_error(pp, "must specify at least one data column"); }
      
    /* Clear out the data file column types: */
    o->type = fdtype_vec_new(nf_max);
    for (uint32_t kf = 0; kf < nf_max; kf++) { o->type.e[kf] = fget_data_type_NOT; }
    
    /* Set type = numeric for all the columns that contain variables. Remember largest col: */
    uint32_t nf = 0; /* Number of data file cols actually used. */
    for (uint32_t kv = 0; kv < nv; kv++) 
      { int32_t kf = o->dataColumns.e[kv];
        if (kf > 0) { defcol(kf, fget_data_type_NUM, FALSE, o->type.e, &nf); }
      }
   
    /* Get the weight column, if any: */
    if (argparser_keyword_present(pp, "-weightColumn"))
      { o->weightColumn = (int32_t)argparser_get_next_int(pp, 1, jpc_MAX_FIELDS); 
        defcol(o->weightColumn, fget_data_type_NUM, TRUE, o->type.e, &nf);
      }
    else
      { o->weightColumn = -1; };
    
    /* Get the point key column, if any: */
    if (argparser_keyword_present(pp, "-keyColumn"))
      { o->keyColumn = (int32_t)argparser_get_next_int(pp, 1, jpc_terms_MAX);
        defcol(o->keyColumn, fget_data_type_ALF, FALSE, o->type.e, &nf);
      }
    else
      { o->keyColumn = -1; };
    
    argparser_get_keyword(pp, "-noise");
    o->noise = argparser_get_next_double(pp, 0.0, 1.0e200);
    
    if (argparser_keyword_present(pp, "-centerFile"))
      { o->centerFile = argparser_get_next_non_keyword(pp); }
    else
      { o->centerFile = NULL; };
    
    if (argparser_keyword_present(pp, "-covarFile"))
      { o->covarFile = argparser_get_next_non_keyword(pp); }
    else
      { o->covarFile = NULL; };
    
    if (argparser_keyword_present(pp, "-compFile"))
      { o->compFile = argparser_get_next_non_keyword(pp); }
    else
      { o->compFile = NULL; };
    
    if (argparser_keyword_present(pp, "-decompFile"))
      { o->decompFile = argparser_get_next_non_keyword(pp); }
    else
      { o->decompFile = NULL; };
    
    if (argparser_keyword_present(pp, "-format"))
      { o->format = argparser_get_next_non_keyword(pp); }
    else
      { o->format = "%+12.5f"; };

    if (argparser_keyword_present(pp, "-verbose"))
      { o->verbose = argparser_get_next_bool(pp); }
    else
      { o->verbose = FALSE; };
    
    /* Check for spurious arguments: */
    argparser_finish(pp);

    return o;
    
    /* ---------------------------------------------------------------------- */
      
    void defcol(uint64_t col, fget_data_type_t tkf,  fget_data_type_t type[], uint32_t *nf_P)
      { 
        demand(col >= 1, "invalid column");
        uint32_t kf = (uint32_t)(col - 1);
        assert(kf < nf_max);
        fget_data_set_field_type(kf, tkf, rep_ok, nf_max, type);
        if (kf >= (*nf_P)) { (*nf_P) = kf + 1; }
      }
  }

vec_typeimpl(fdtype_vec_t, fdype_vec, fget_data_type_t);
