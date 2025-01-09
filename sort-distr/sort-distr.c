/* Last edited on 2024-12-25 09:54:35 by stolfi */

// .\" Last edited on 2003-09-25 15:50:49 by stolfi
// .\" See the authorship and copyrigh notice at the end of this file.
// .\"
// .nh
// .TH sort-distr 1
// .SH Name
// sort-distr \- sort objects by similarity of distributions
// 
// .SH Syntax
// .B sort\-distr \-\-numValues 
// .I N 
// .B [ \-\-width 
// .I W 
// .B ] [ \-\-skip 
// .I S
// .B ] [ \-\-discrete ] [ \-\-geometric ] [ \-\-farSort | \-\-cluSort ] [ \-\-delReSort ] [ \-\-repeat 
// .I COUNT 
// .B ] [ \-\-showMatrix 
// .I FILE 
// .B ] [ \-\-showPicture 
// .I FILE 
// .B ] [ \-\-verbose ]
// 
// .SH Description
// .B sort\-distr
// sorts a set of probability distributions by similarity.
// Each line of stdin should contain a list of 
// .I N
// numbers, which are interpreted as a probability distribution on
// .I N
// alternatives or events.  The program tries to reorder
// those lines so that lines with similar distributions
// are close together.
// .PP
// This goal is equivalent to the Traveling Salesman problem, so don't 
// expect the output to be the best possible solution.
// .SH Options
// .TP 6
// .BI \-n\  NUM
// .TP 0
// .BI \-\-numValues\  NUM
// The number of probability values in each distribution.
// This parameter is required. 
// 
// .TP 6
// .BI \-s\  NUM
// .TP 0
// .BI \-\-skip\  NUM
// Skip the first 
// .I NUM 
// bytes of every line, before parsing the first probability value.
// The default is zero. 
// 
// .TP 6
// .BI \-w\  NUM
// .TP 6
// .BI \-\-width\  NUM
// Specifies that the probabilities are given in fixed format,
// as
// .I N 
// adjacent fields, each  
// .I NUM 
// bytes wide.
// If this parameter is not specified,
// .B sort\-distr
// assumes that the probabilities are given in free format.
// In that case each field must be terminated by 
// a space, comma, or newline, and may be surrounded 
// by extra spaces.  In either case, fields consisting entirely
// of periods or dashes are taken to be zero.
// 
// .TP 6
// .B \-d
// .TP 6
// .B \-\-discrete
// If this option is given,
// .B sort\-distr
// assumes that the numbers in each line are observed counts
// in some sampling experiment, and converts them to probabilities
// taking into account sampling error.
// 
// .TP 6
// .B \-g
// .TP 6
// .B \-\-geometric
// By default,
// .B sort\-distr
// assumes the events are isolated, and
// compares distributions as vectors with the Euclidean metric.
// If this option is given, however,
// .B sort\-distr
// assumes that the events correspond to equally spaced positions
// along a unidimensional path, and uses the "earth-movers"
// distance to compare them. 
// 
// .TP 6
// .B \-fs
// .TP 6
// .B \-\-farSort
// Specifies that the first sorting pass should be done with 
// the
// .I Farthest-Point 
// alternate sorting heuristic. 
// 
// .TP 6
// .B \-cs
// .TP 6
// .B \-\-cluSort
// Specifies that the first sorting pass should be done with 
// the
// .I Cluster-Tree
// alternate sorting heuristic. 
// 
// .TP 6
// .B \-ds
// .TP 6
// .B \-\-delReSort
// Specifies that the second and additional sorting passes should be done with 
// the
// .I Delete-And-Insert 
// alternate re-sorting heuristic. 
// 
// .TP 6
// .BI \-m\  FILE
// .TP 6
// .BI \-\-showMatrix\  FILE
// Print the distance matrix to the named file 
// (or standard output if 
// .I FILE 
// is "-").
// The distances are scaled to the range [0..99.5] and then
// rounded to the nearest integer. 
// 
// .TP 6
// .BI \-p\  FILE
// .TP 6
// .BI \-\-showPicture\  FILE
// Writes the distance matrix to the named file 
// (or standard output if 
// .I FILE 
// is "-"), as a color PPM 
// image.  Each element becomes a colorful pixel, whose
// brightness decreases with increasing distance (white for distance 0,
// black for distance 1).
// 
// .TP 6
// .BI \-r\  COUNT
// .TP 6
// .B \-\-repeat\  COUNT
// Number of sorting passes (default 2). The first pass finds the 
// two items with largest distance, places them at opposite ends of the 
// list, and sorts all other items based on the ratio of their distances to those
// two. The second pass looks at each item in turn, and tries to move it
// along the list so as to decrease the path length.  If 
// .I COUNT
// is greater than 2, the 
// second pass is performed
// .I COUNT \- 1
// times. If 
// .I COUNT
// is 1, performs only the first pass.
// If 
// .I COUNT 
// is zero, the items are not sorted and output is supressed.
// (This choice is useful for generating a distance matrix 
// or picture of already sorted items.)
// 
// .TP 6
// .B \-v
// .TP 6
// .B \-\-verbose
// Print various progress messages to standard error.
// 
// .SH Bugs
// Likely.
// 
// .SH Challenges for the bored
// Add other distribution metrics.
// .PP
// Improve the sorting algorithm.
// .PP
// Ad dendrogram output option.
// 
// .SH Author
// Jorge Stolfi
// 
// .\" (* (C) Copyright 1992 Universidade Estadual de Campinas (UNICAMP)           *)
// .\" (*                    Campinas, SP, Brazil                                  *)
// .\" (*   Jorge Stolfi        - CS Dept, UNICAMP <stolfi@dcc.unicamp.br>         *)


/* See the authorship and copyright notice at the end of this file. */
/* Changes:

     97-12-07  Changed --binSort heuristic to break the path preferably
               at long edges.
     
     97-11-30  Changed --binSort heuristic to break the path preferably
               at long edges.
     
     97-11-29  --showPicture now precomputes the color table so as
               to use only 216 different colors.
     
     97-11-27  Improved (?) sorting algorithm, with DiaSortItems
               or FarSortItems followed by UniReSortItems or DelReSortItems.
               Changed --showPicture and --showMatrix to take a file name.

     97-11-26  Added --showPicture option.
               Also improved the --showMatrix format.
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <assert.h>

#include <jsrandom.h>
#include <bool.h>

/* LIMITS */

#define MAX_ITEMS           (2048)
#define MAX_VALUES          (2048)
#define MAX_WIDTH           (255)
#define MAX_LINE_LENGTH     (2048 + MAX_VALUES*MAX_WIDTH)
#define MAX_SKIP            MAX_LINE_LENGTH
#define MAX_REPEAT          (10000)
#define NUM_COLORS          (216)

/* INTERNAL PROTOTYPES */

int main(int argc, char **argv);

void* Alloc(long size);
  /*
    Like "malloc", but bombs out if runs out of memory. 
  */
  
char *CopyString(char *p);
  /* 
    Makes a copy of a null-terminated string. 
  */

bool_t ReadItem(char **L, long *N, long *A, long R);
  /* 
    Returns TRUE if "rd" (whose filename is "name") was at
    end-of-file.  Otherwise reads the next line of "rd" into "*L", and
    stores its length into "*N".  Assumes the buffer "*L" is currently
    "*A" bytes long; reallocates a new one if necessary.  Stops with
    error if the last line does not end with newline.
    Reports "R" as line number in case of errors.
  */
  
char* CopyItem(char *L, long N);
  /*
    Allocates an "N"-byte string and copies "L[0..N-1]" into it.
  */

void WriteItem(char *L, long N);
  /*
    Writes bytes "L[0..N-1]" to stdout, followed by newline.
  */
  
void AllocItemEntry(long i, char ***L, long **N, long *A);
  /*
    Makes sure elements "(*L)[i]" and "(*N)[i]" exist,
    reallocating "*L" and "*N" if necessary.
    Assumes their current allocated size is "*A".
  */
 
void ParseFreeVector(char *L, long N, double *v, long M, long R);
  /*
    Parses "L[0..N-1]" as "M" real numbers "v[0..M-1]" in free format.
    The values be terminated by blanks or end-of-buffer.
    Ignores any "L" bytes beyond "v[M-1]". 
    
    Bombs out if there are missing numbers or format
    errors, reporting "R" as the line number.
  */
  
void ParseFixedVector(char *L, long N, long W, double *v, long M, long R);
  /*
    Parses "L[0..N-1]" as "M" real numbers "v[0..M-1]", each "W" bytes wide.
    Ignores any "L" bytes beyond "v[M-1]". 
    
    Bombs out on "N < M*W" or format error, reporting "R" as the line
    number.
  */
  
double ParseNum(char *L, long W, char *fmt, long R, long K);
  /* 
    Parses "L[0..W-1]" as a real number.  
    If the field contains only blanks, tabs, periods, or dashes,
    returns 0; else parses the field using the format "fmt".
    
    Bombs out in case of format error, reporting "R" 
    as the line number and "K" as the element index.
  */

void NormalizeProbs(double *a, long M, long R);
  /*
    Scales the vector "a" so that it has unit "sum".

    Bombs out if "a"'s sum is zero or any "a[i]" is negative,
    reporting "R" as the line number.
  */

void EstimateProbs(double *a, long M, long R);
  /*
    Assumes the "a[i]" are integer event counts, and
    converts them to estimated probabilities by the
    formula "a[i] = (a[i]+1)/(sum + M)" where "sum" 
    is the sum of all "a[i]".

    Bombs out if any "a[i]" is negative,
    reporting "R" as the line number.
  */

double L2ProbDist(double *a, double *b, long M);
  /*
    Computes the Euclidean distance between vectors "a[0..M-1]" and
    b[0..M-1]", asssuming they are normalized to unit sum.  The
    distance is divided by sqrt(2) to give a result in "[0_1]".
  */
    
double EarthMoverDist(double *a, double *b, long M);
  /*
    Computes the "earth mover's" distance between vectors "a[0..M-1]"
    and b[0..M-1]", assuming they have been normalized to unit sum.
    The result is a number in "[0_1]".
  */
  
void DiaSortItems(long *it, long N, double **d);
  /*
    Permutes "it[0..N-1]" so that "it[0]" and "it[N-1]"
    is a diameter, and the remaining items are sorted by 
    the relative distances to those two.
  */
  
void FarSortItems(long *it, long N, double **d);
  /*
    Sorts "it[0..N-1]" by repeatedly taking the element
    furthest from the path and inserting it in the best place.
  */
  
void CluSortItems(long *it, long N, double **d);
  /*
    Builds a dendrogram tree whose leaves are the items in "it[0..N-1]".
    Then swaps and/or reverses the children of every node so that 
    the jump between them is minimized. Finally enumerates the 
    leaves in this sorted tree and stores them back into "it". 
  */
  
void UniReSortItems(long *it, long N, double **d);
  /*
    Scans "it[0..N-1]" and tries to move each element
    so as to reduce the total path length. 
  */
  
void DelReSortItems(long *it, long N, double **d);
  /*
    Removes a big subset of the indices, and inserts them 
    back one by one. 
  */
  
long InsertItemInPath(long *it, long *Np, long ki, double **d);
  /*
    Inserts item "ki" in the best place of the path "it[0..*Np-1]".
    Assumes there is space in "it[]" for one more element.
    Increments *Np.
  */
  
/* ORDERED CLUSTER TREES:
  
  The procedures SortClusterTree and DumpClusterTree below take an
  integer N and two vectors "uch[0..N-2]" and "vch[0..N-2]", and interpret them 
  as an ordered binary trees with leaves "0..N-1" and internal nodes "N..2*N-2".
  
  A further signed parameter "R" specifies a subtree of that tree, and 
  a direction for its traversal, as follows: 
  
    If "R" is in the range "[0 .. N-1]", the subtree consists of the
    single leaf node number "R".

    If "R" is in the range "[N .. 2*N-2]", the subtree's root is the
    internal node number "R", whose left and right subtrees are
    "uch[R-N]" and "vch[R-N]", respectively.

    If "R < 0", the subtree in question is a left-to-right reversal of
    the subtree "-R".  In particular, if "R <= -N", the
    left and right subtrees are "-vch[-R-N]" and "-uch[-R-N]",
    respectively.

  The pointers "uch[]" and "vch[]" may be negative, in which case
  they too are interpreted accoridng to the last rule above.
  
*/
  
void SortClusterTree(long *uch, long *vch, long R, long N, double **d, long *kap, long *kzp);
  /*
    Interprets "N", "uch", "vch" as an ordered binary tree on the
    items "0..N-1", as described above.  Traverses the subtree "R",
    swapping and/or reversing the two children of every internal node,
    if necessary, so as to minimize the length of the jump between
    them.  Also returns in "*kap" and "*kzp" the first and last leaves
    of the sorted subtree.
  */

void DumpClusterTree(long *uch, long *vch, long R, long N, long *it, long *ip);
  /*
    Interprets "N", "uch", "vch" as an ordered binary tree on the
    items "0..N-1", as described above.  Traverses the subtree "R"
    in its internal order, and stores its leaves in "it[*ip...]", 
    incrementing "*ip".
  */

void WriteMatrix(FILE *f, double **d, long *it, long N, char **L, long skip);
  /*
    Prints the distance matrix "d[it[0..N-1]][it[0..N-1]]" to file "f".
    Uses "L[it[0..N-1]][0..skip-1]" as the row and column labels. */
    
int FormatDist(double d);
  /* 
    Converts a distance value from the range "[0 _ 1]" to 
    an integer in the range "[0..99]". */
    
void WritePicture(FILE* f, double **d, long *it, long N);
  /*
    Writes the distance matrix "d[it[0..N-1]][it[0..N-1]]" to 
    file "f", as a raw PPM image.  Distance 1.0 
    is mapped to black, 0.0 to white, and intermediate values to
    intermediate intensities (with varying hues). */
    
FILE *OpenWrite(char *name);
  /* 
    Opens a file for writing; returns "stdout" if "name" is "-".
    Bombs in case of failure.
  */
    
typedef struct { byte R, G, B; } pixel;

void SelectColors(pixel *p, long M);
  /*
    Stores in "p[0..M-1]" a palette of colors
    suitable for displaying values in the range "[0_1]". */
  
typedef struct {
    long numValues;    /* Number of values in each distribution. */
    long width;        /* Width of each value, or 0 for free format. */
    long skip;         /* Number of bytes to skip before first value. */
    long repeat;       /* Number or times to apply "UniReSortItems" */
    bool_t discrete;     /* TRUE if inputs are counts, FALSE if they are probs. */
    bool_t geometric;    /* TRUE for earth mover's distance, FALSE for Euclidean. */
    bool_t cluSort;      /* TRUE to begin with CluSortItems. */
    bool_t farSort;      /* TRUE to begin with FarSortItems. */
    bool_t delReSort;    /* TRUE to use DelReSortItems. */
    char *matName;     /* File name for distance matrix, or NULL if none. */
    char *picName;     /* File name for PPM distance image; or NULL if none. */
    bool_t verbose;      /* TRUE to mumble while working. */
  } Options;
  
void GetOptions(int argc, char **argv, Options *op);
  /* 
    Parses the command line options.
  */

bool_t GetLongOption(
    int argc, char **argv,
    int *ip,
    char *brief,
    char *wordy,
    long min,
    long max,
    long *vp
  );
  /*
    If "argc[i]" matches either the "brief" or the "wordy" keyword, and is
    followed by one more parameter, parses that parameter as a long
    integer, stores it in "*vp", increments "*i" by 2, and returns
    TRUE. Otherwise leaves "*vp" unchanged and returns FALSE. Bombs if
    the value is not an integer between "min" and "max".
  */
     
bool_t GetCharOption(
    int argc, char **argv,
    int *ip,
    char *brief,
    char *wordy,
    char **vp
  );
  /*
    If "argc[i]" matches either the "brief" or the "wordy" keyword, and is
    followed by one more parameter, then makes a copy of that paramter,
    stores its addres in "*vp", increments "*i" by 2, and returns
    TRUE. Otherwise leaves "*vp" unchanged and returns FALSE.
  */
     
bool_t GetBoolOption(
    int argc, char **argv,
    int *ip,
    char *brief,
    char *wordy,
    bool_t *vp
  );
  /*
    If "argc[i]" matches either the "brief" or the "wordy" keyword,
    stores TRUE in "*vp", increments "*i", and returns TRUE. Otherwise
    leaves "*vp" unchanged and returns FALSE.
  */

void Error(char *msg);
void ParamError(char *msg, char *key, long K);
void FileError(char *msg, char *name);
void LineError(char *msg, long R);
void ElementError(char *msg, long R, long K);
  /*
    Print error message about input record "R" and element/parameter "K",
    then die with nonzero status.
  */

/* PROCEDURES */

static char *progName;

int main(int argc, char **argv)
  {
    /* Command line options: */
    Options o;

    /* Read buffer: */
    char *bufL;        /* Last item of input file. */
    long bufN;         /* Used size of "bufL". */
    long bufA;         /* Allocated size of "bufL". */

    /* The input data: */
    long nItems;       /* Number of items stored so far. */
    char **L;          /* "L[0..nItems-1i]" are pointers to the items read so far. */
    long *N;           /* "N[0..nItems-1]" are their lengths (used and allocated). */
    long aItems;       /* Allocated size of "L", "N", "v", and "s". */

    /* The numerical data: */
    double **v;        /* "v[i][0..numValues-1]" is the parsed distribution of item "i". */
    double **d;        /* "d[i,j]" is the distance between item "i" and item "j" */

    /* The output order: */
    long *it;         /* "it[k]" is the "k"th item in output order. */

    /* Misc variables: */ 
    long i, j, k;
    long minItemLength;

    /* Collect options: */
    progName = argv[0];
    GetOptions(argc, argv, &o);
    minItemLength = o.skip + (o.width == 0 ? (2*o.numValues - 1) : (o.numValues * o.width));

    /* Initialize: */

    bufA = 0;
    nItems = 0;
    aItems = 0;

    /* Read data: */
    while(! ReadItem(&bufL, &bufN, &bufA, nItems+1))
      { if (bufN < minItemLength) 
          { LineError("line too short", nItems+1); }
        else
          { AllocItemEntry(nItems, &L, &N, &aItems);
            L[nItems] = CopyItem(bufL, bufN);
            N[nItems] = bufN;
            nItems++;
          }
      }
    if (o.verbose) { fprintf(stderr, "%6ld items read.\n", nItems); }

    /* Extract the numerical vectors: */
    if (o.verbose) { fprintf(stderr, "parsing distributions...\n"); }
    v = (double**)(Alloc(nItems*sizeof(double*)));
    for (i=0; i<nItems; i++)
      { v[i] = (double*)(Alloc(o.numValues*sizeof(double)));
        if (o.width == 0)
          { ParseFreeVector(&(L[i][o.skip]), N[i]-o.skip, v[i], o.numValues, i+1); }
        else
          { ParseFixedVector(&(L[i][o.skip]), N[i]-o.skip, o.width, v[i], o.numValues, i+1); }
       
        if (o.discrete) 
          { EstimateProbs(v[i], o.numValues, i+1); }
        else
          { NormalizeProbs(v[i], o.numValues, i+1); }
      }

    /* Compute distance matrix: */
    /* We should use a triangular matrix, but it is simpler to keep it square. */
    if (o.verbose) { fprintf(stderr, "computing distance matrix...\n"); }
    d = (double**)(Alloc(nItems*sizeof(double*)));
    for (i=0; i<nItems; i++)
      { d[i] = (double*)(Alloc(nItems*sizeof(double)));
        for (j=0; j<i; j++)
          { double dd;
            if (o.geometric) 
              { dd = EarthMoverDist(v[i], v[j], o.numValues); }
            else
              { dd = L2ProbDist(v[i], v[j], o.numValues); }
            d[i][j] = dd;
            d[j][i] = dd;
          }
        d[i][i] = 0;
      }
    /* We could free the "v" vectors now. */

    /* Allocate index vector: */
    it = (long*)(Alloc(nItems*sizeof(long)));
    for (k=0; k<nItems; k++) it[k] = k;
    
    /* Decide what to do: */
    if (o.repeat > 0)
      { /* Sort and print items: */
        long r;
        if (o.verbose) { fprintf(stderr, "initial sorting...\n"); }
        if (o.farSort) 
          { FarSortItems(it, nItems, d); }
        else if (o.cluSort) 
          { CluSortItems(it, nItems, d); }
        else
          { DiaSortItems(it, nItems, d); }
        if (o.verbose) { fprintf(stderr, "refining order... "); }
        for (r = 2; r <= o.repeat; r++) 
          { if (o.verbose) { fprintf(stderr, "+"); }
            if (o.delReSort) 
              { DelReSortItems(it, nItems, d); }
            else
              { UniReSortItems(it, nItems, d); }
          } 
        if (o.verbose) { fprintf(stderr, "\n"); }
        for (k=0; k<nItems; k++) 
          { long i = it[k];
            WriteItem(L[i], N[i]);
          }
      }
    
    if (o.matName != NULL) 
      { /* Print distance matrix: */
        if (o.verbose) { fprintf(stderr, "writing matrix...\n"); }
        WriteMatrix(OpenWrite(o.matName), d, it, nItems, L, o.skip);
      }
    
    if (o.picName != NULL) 
      { /* Print PBM image of distance matrix: */
        if (o.verbose) { fprintf(stderr, "writing picture...\n"); }
        WritePicture(OpenWrite(o.picName), d, it, nItems);
      }
      
    fclose(stdout);
    return(0);
  }

void* Alloc(long size)
  { void *p = malloc(size);
    if (p == NULL) { Error("out of memory"); }
    return(p);
  }
  
char *CopyString(char *p)
  { char *q = (char*)(Alloc((strlen(p)+1)*sizeof(char)));
    char *t = q;
    while(1) { (*t) = (*p); if ((*t)==0) return q; p++; t++; } 
  }

bool_t ReadItem(char **L, long *N, long *A, long R)
  { int c;
    char *curL = (*L);
    long curA = (*A);
    long curN = 0;
    c = fgetc(stdin);
    if (c == EOF) { (*N) = 0; fclose(stdin); return(TRUE); }
    while (c != '\n')
      { if (curN >= curA) 
          { if (curN >= MAX_LINE_LENGTH) 
              { LineError("line too long", R); }
            else
              { long newA = (curA == 0 ? 256 : 2*curA);
                char *newL = (char *)Alloc(newA*sizeof(char));
                long i;
                for (i=0; i<curN; i++) newL[i] = curL[i];
                if (curA > 0) free(curL);
                curA = newA;
                curL = newL;
              }
          }
        curL[curN] = c;
        curN++;
        c = fgetc(stdin);
        if (c == EOF) { Error("missing NL at end of input file"); }
      }
    (*L) = curL;
    (*A) = curA;
    (*N) = curN;
    return(FALSE);
  } 

char* CopyItem(char *L, long N)
  { char *S = (char*)(Alloc(N*sizeof(char)));
    long i;
    for (i=0; i<N; i++) S[i] = L[i];
    return(S);
  }

void WriteItem(char *L, long N)
  { long i;
    for (i=0; i<N; i++) putchar(L[i]);
    putchar('\n');
  }

void AllocItemEntry(long i, char ***L, long **N, long *A)
  {
    if (i >= (*A))
      { 
        if (i >= MAX_ITEMS)
          { LineError("too many lines", i+1); }
        else 
          { char **curL = (*L);
            long *curN = (*N);
            long curA = (*A);
            long newA = (curA == 0 ? 256 : 2*curA);
            char **newL = (char **)Alloc(newA*sizeof(char*));
            long *newN = (long *)Alloc(newA*sizeof(long));
            long j;
            for (j=0; j<curA; j++) { newL[j] = curL[j]; newN[j] = curN[j]; }
            if (curA > 0) { free(curL); free(curN); }
            (*L) = newL;
            (*N) = newN;
            (*A) = newA;
          }
      }
  }

void ParseFreeVector(char *L, long N, double *v, long M, long R)
  { long i;
    char *eb = L;
    long nb = N;
    char fmt[10];
           
    for (i=0; i<M; i++)
      { /* At this point "eb[0..nb-1]" is the unparsed part of "L" */
        char *ef;
        long nf;
        /* Skip blanks and tabs: */
        while ((nb > 0) && (((*eb) == ' ') || ((*eb) == '\011'))) { eb++; nb--; }
        if (nb <= 0) { LineError("too few elements", R); }
        /* Find next blank, tab, comma or end-of-buffer: */
        ef = eb; nf = nb;
        while ((nf > 0) && ((*ef) != ',') && ((*ef) != ' ') && ((*ef) != '\011')) { ef++; nf--; }
        sprintf(fmt, "%%%ldlf", nb-nf);
        /* Parse element: */
        v[i] = ParseNum(eb, nb-nf, fmt, R, i);
        /* Skip it: */
        eb = ef; nb = nf;
        /* skip comma, if any: */
        if ((nb > 0) && ((*eb) == ',')) { eb++; nb--; }
      }
  }
  
void ParseFixedVector(char *L, long N, long W, double *v, long M, long R)
  { long i;
    char *eb = L;
    long nb = N;
    char fmt[10];
    sprintf(fmt, "%%%ldlf", W);
           
    for (i=0; i<M; i++)
      { /* At this point "eb[0..nb-1]" is the unparsed part of "L" */
        if (nb < W) { LineError("line too short", R); }
        v[i] = ParseNum(eb, W, fmt, R, i);
        eb += W; nb -= W;
      }
  }
  
double ParseNum(char *L, long W, char* fmt, long R, long K)
  { long i;
    for (i=0; i<W; i++)
      { if ((L[i] != ' ') && (L[i] != '.') && (L[i] != '-')) 
          { double val; int res;
            res = sscanf(L, fmt, &val); 
            if (res != 1) { ElementError("format error", R, K); }
            return(val);
          }
      }
    return 0;
  }

double L2ProbDist(double *a, double *b, long M)
  {
    double s = 0.;
    long i;
    for (i=0; i<M; i++) { double d = a[i] - b[i]; s += d*d; }
    return sqrt(s/2);
  }
      
double EarthMoverDist(double *a, double *b, long M)
  {
    double s = 0;
    double sa = 0;
    double sb = 0;
    long i;
    if (M==0) return (0);
    for (i=0; i<M; i++) 
      { double d;
        sa += a[i]; sb += b[i];
        /* "sa" is the earth mass in "a[0..i]", ditto for "sb" */
        d = fabs(sa - sb); 
        /* "d" is the mass that must be carried from "[i]" to "[i+1]" */
        s += d;
      }
    return sqrt(s/((double)(M-1)));
  }

void NormalizeProbs(double *a, long M, long R)
  { double s = 0;
    long i;
    for (i=0; i<M; i++) 
      { if (a[i] < 0) { ElementError("negative probability", R, i); }
        s += a[i];
      }
    if (s <= 0) { LineError("null distribution", R); }
    for (i=0; i<M; i++) a[i] /= s;
  }

void EstimateProbs(double *a, long M, long R)
  { double s = 0;
    long i;
    for (i=0; i<M; i++) 
      { if (a[i] < 0) { ElementError("negative count", R, i); }
        s += a[i];
      }
    s += (double)M;
    for (i=0; i<M; i++) a[i] = (a[i]+1)/s;
  }

void DiaSortItems(long *it, long N, double **d)
  {
    /* Sorting keys: */
    double *rd;
    long ka, kb;
    if (N <= 2) return;
    /* Find a diametral pair ka,kb of items, put it in it[0], it[N-1]: */
    { long i, j, iMax, jMax;
      double dMax = -1;
      for (i=1; i<N; i++)
        for (j=0; j<i; j++)
          { double dd = d[it[i]][it[j]]; 
            if (dd > dMax) { dMax = dd; iMax = i; jMax = j; }
          }
      ka = it[iMax]; it[iMax] = it[0];   it[0] = ka;
      kb = it[jMax]; it[jMax] = it[N-1]; it[N-1] = kb;
    }
    /* Define sort keys based on distances to ka and kb: */
    { long i;
      double eps = 0.000001 * d[ka][kb];
      if (eps == 0) return;
      rd = (double *)(Alloc(N*sizeof(double)));
      rd[0] = -1; rd[N-1] = +1;
      for (i=1; i<N-1; i++)
        { long ki = it[i];
          double da = d[ki][ka]; 
          double db = d[ki][kb];
          rd[i] = (da-db)/(da+db+eps);
          if ((rd[i] <= -1) || (rd[i] >= +1)) { Error("program error - rd[i] out of range"); }
        }
    }
    /* Sort remaining items "it[1..N-2]" by ratio of distances to those two: */
    { long i;
      for (i=1; i<N-1; i++)
        { long ki = it[i];
          double ri = rd[i];
          long j = i;
          while (rd[j-1] > ri) { rd[j] = rd[j-1]; it[j] = it[j-1]; j--; }
          if (i != j) { rd[j] = ri; it[j] = ki; }
        }
    }
    if (rd != NULL) free(rd);
  }

void FarSortItems(long *it, long N, double **d)
  {
    /* Distances from path: */
    long ka, kb;
    if (N <= 2) return;
    /* Find a diametral pair of items ka,kb, put them in it[0], it[1]: */
    { long i, j, iMax, jMax;
      double dMax = -1;
      for (i=1; i<N; i++)
        for (j=0; j<i; j++)
          { double dd = d[it[i]][it[j]]; 
            if (dd > dMax) { dMax = dd; iMax = i; jMax = j; }
          }
      ka = it[iMax]; it[iMax] = it[0]; it[0] = ka;
      kb = it[jMax]; it[jMax] = it[1]; it[1] = kb;
    }
    /* Repeately look for the farthest item, insert it in path: */
    { long i, j;
      double *dd = (double *)(Alloc(N*sizeof(double)));
      /* Initialize the distance-from-path and best-place tables: */
      /* If current path if "it[0..i-1]", then */
      /*   The candidates for insertion are "it[i..N-1]" */
      /*   Element "dd[j]" is distance of "it[j]" to the path *vertices* */
      dd[0] = 0;
      dd[1] = 0;
      for (j=2; j<N; j++) 
        { long kj = it[j];
          double da = d[kj][ka]; 
          double db = d[kj][kb];
          dd[j] = (da < db ? da : db);
        }
      /* Now expand that path: */
      i = 2;
      while(i < N)
        { long j, jMax, ki, where;
          double ddMax;
          /* Look for item "it[jMax]" with largest "dd[jMax]": */
          ddMax = -1;
          for (j=i; j<N; j++)
            { double ddj = dd[j]; 
              if (ddj > ddMax) { ddMax = ddj; jMax = j; }
            }
          /* Save that item in "ki" and move the current "it[i]" there: */
          /* We don't need "dd[jMax]" anymore: */
          ki = it[jMax]; 
          it[jMax] = it[i]; dd[jMax] = dd[i]; 
          dd[i] = 0;
          where = InsertItemInPath(it, &i, ki, d);
          for (j=i; j<N; j++)
            { double dij = d[ki][it[j]];
              if (dij < dd[j]) { dd[j] = dij; }
            } 
        }
      if (dd != NULL) free(dd);
    }
  }

void CluSortItems(long *it, long N, double **d)
  {
    /* The first pass builds a tree. */
    /* The nodes of the tree are numbered in the range "[0..2*N-2]". */
    /* If "i" is in "[0..N-1]", node "i" is a leaf, which is item "i". */
    /* If "i" is in "[N..2*N-2]", node "i" is an internal node, whose children */
    /* are nodes "uch[i-N]" and "vch[i-N]". */
    long M = 2*N-1;
    long *uch = (long *)Alloc((M-N)*sizeof(long));
    long *vch = (long *)Alloc((M-N)*sizeof(long));
    
    /* During construction we have K forests whose roots nodes are "root[0..K-1]": */
    long K;
    long *root = (long *)Alloc(N*sizeof(long));
    
    /* The distance between forest "f" and forest "g" is stored in "e[i][j]": */
    double **e = (double **)Alloc(N*sizeof(double*));
    
    /* Initialize root set: */
    { long i; for (i=0; i<N; i++) { root[i] = it[i]; } }
    K = N;
    
    /* Intialize cluster distance matrix: */
    { long i, j;
      for (i=0; i<N; i++) 
        { e[i] = (double*)Alloc(N*sizeof(double));
          for (j=0; j<N; j++) { e[i][j] = d[i][j]; }
        }
    }
    
    /* Repeatedly merge the two nearest clusters: */
    while (K > 1)
      { /* The internal nodes created so far are [N..M-K] */
        /* Find the two nearest clusters */
        long iMin, jMin;
        { double dMin = 999;
          long i, j;
          for (i=0; i<K; i++)
            for (j=i+1; j<K; j++)
              { double dij = e[i][j];
                if (dij < dMin) { dMin = dij; iMin = i; jMin = j; }
              }
        }
        /* Merge them into a new cluster: */
        { long r = M-K+1; long k;
          uch[r-N] = root[iMin]; 
          vch[r-N] = root[jMin];
          /* Store the new root in root[iMin], */
          /* and compute its distances to other subtrees: */
          root[iMin] = r;
          for (k=0; k<K; k++)
            { if ((k != iMin) && (k != jMin))
                { double dik = e[iMin][k]; double djk = e[jMin][k];
                  /* double drk = (dik + djk + (dik < djk ? dik : djk)) / 3; */
                  double drk = 2/(1/dik + 1/djk);
                  e[iMin][k] = drk; e[k][iMin] = drk;

                }
             }
          /* Eliminate root[jMin]: */
          root[jMin] = root[K-1];
          K--;
          for (k=0; k<K; k++)
            { if (k != jMin)
                { double duk = e[K][k]; e[jMin][k] = duk; e[k][jMin] = duk; }
            }

        }
      }

    /* Second pass: sort the subtrees by proximity */
    { long ka, kz; SortClusterTree(uch, vch, +(M-1), N, d, &ka, &kz); }
    
    /* Third pass: dump it out: */
    { long i = 0;
      DumpClusterTree(uch, vch, +(M-1), N, it, &i);
      if (i != N) { Error("program error: DumpClusterTree"); }
    }
    
    free(e); free(uch); free(vch);
  }

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
  
void DelReSortItems(long *it, long N, double **d)
  {
    long i;
    if (N <= 2) return;
    /* Split path into two subsets "A = it[0..i-1]" and "B = it[i..N-1]": */
    srand(62712357);
    { long j;
      bool_t keep = TRUE;
      double db = d[it[0]][it[1]];
      double da = db;
      i = 0;
      for(j=1; j<N; j++)
        { /* The B items are "it[i..j-1]" */
          /* Break preferentially at long edges: */
          if (j == N-1)
            { keep = TRUE; }
          else 
            { double dc = ((j > N-2) ? db : d[it[j]][it[j+1]]);
              double er = db/(da + db + dc + 1.0e-100);
              double prob = 0.125 + 0.375 * er;
              if (rand() < (int)(prob*32767)) 
                { keep = !keep; 
                  /* fprintf(stderr, "keep=%d at j=%ld\n", keep, j); */
                }
              da = db; db = dc;
            }
          if (keep) 
            { /* keep it in the path */
              long k = it[i]; it[i] = it[j]; it[j] = k; i++;
            }
        }
    }
    /* Insert the B items into the A path: */
    while (i<N)
      { long ki = it[i]; long where;
        where = InsertItemInPath(it, &i, ki, d);
      }
  }
  
long InsertItemInPath(long *it, long *Np, long ki, double **d)
  {
    long p, pMin;
    long N = (*Np);
    double ddMin;
    /* Find the best edge "(ka,kb) = (it[pMin-1],it[pMin])" where to insert "ki": */
    ddMin = 999;
    for (p=0; p<=N; p++)
      { long ka = (p > 0 ? it[p-1] : -1);
        long kb = (p < N ? it[p] : -1);
        double da = (ka > 0 ? d[ki][ka] : 0); 
        double db = (kb > 0 ? d[ki][kb] : 0);
        double de = ((ka > 0) && (kb > 0) ? d[ka][kb] : 0);
        double ddp = da + db - de;
        if (ddp < ddMin) { ddMin = ddp; pMin = p; }
      }
    /* Open up that edge: */
    p = N;
    while (p > pMin) { it[p] = it[p-1]; p--; }
    it[pMin] = ki;
    N++;
    /* Update path length: */
    (*Np) = N;
    return (pMin);
  }
  
void UniReSortItems(long *it, long N, double **d)
  {
    long i;
    if (N <= 2) return;
    /* Take each item and see if it would fit better elsewhere */
    i = 0;
    while (i<N)
      { long ki = it[i];
        long jMin;
        double ci, cjMin; 
        /* Compute "ci" = gain in removing node "it[i]" from the path: */
        if (i==0) 
          { ci = d[ki][it[1]]; }
        else if (i==N-1)
          { ci = d[ki][it[N-2]]; }
        else
          { ci = d[ki][it[i-1]]+d[ki][it[i+1]] - d[it[i-1]][it[i+1]]; }
        /* Compute best edge "(it[jMin-1],it[jMin])" */
        /* and the cost "cjMin" of inserting "it[i]" into it: */
        { double da, db, cj; 
          long ja, jb, ka, kb;
          /* Consider inserting "it[i]" before "it[0]": */
          ja = (i == 0 ? 1 : 0);
          ka = it[ja]; 
          da = d[ki][ka];
          cjMin = da; jMin = ja;
          for (jb=ja+1; jb<N; jb++)
            { /* skip over "it[i]": */
              if (jb == i) continue;
              /* At this point "ka" is "it[ja]" where "ja" is */
              /* the previous value of "jb", skipping "i". */
              /* Consider inserting "it[i]" between "it[ja]" and "it[jb]": */
              kb = it[jb]; db = d[ki][kb];
              cj = da + db - d[ka][kb];
              if (cj < cjMin) { cjMin = cj; jMin = jb; }
              /* Prepare for next iteration: */
              ka = kb; da = db; ja = jb;
            }
          /* Consider inserting "it[i]" after end of path: */
          if (da < cjMin) { cjMin = da; jMin = N; }
        }
        /* Now dispose of "it[i]": */
        if (jMin == i) 
          { Error("program error - jMin=i"); }
        else if (ci > cjMin)
          { /* Move "ki = it[i]" to between "it[jMin-1" and "it[jMin]": */
            /* Don't increment "i", because "it[i]" changed. */
            if (jMin < i)
              { long r = i;
                while (r > jMin) { it[r] = it[r-1]; r--; }
                it[jMin] = ki;
              }
            else if (jMin > i+1)
              { long r = i+1;
                while (r < jMin) { it[r-1] = it[r]; r++; }
                it[jMin-1] = ki;
              }
            else
              { Error("program error - jMin = i+1, cjMin < ci"); }
          }
        else
          { i++; }
      }
  }
  
FILE *OpenWrite(char *name)
  { if (strcmp(name, "-") == 0) 
      { return (stdout); }
    else
      { FILE *f = fopen(name, WMODE);
        if (f == NULL) { FileError("can't open file ", name); }
        return(f);
      }
  }
    

void WriteMatrix(FILE *f, double **d, long *it, long N, char **L, long skip)
  {
    long i, j;
    /* Write header: */
    for(i=0;i<skip;i++)
      { for (j=0;j<skip;j++) putc(' ', f);
        putc(' ', f);
        for (j=0;j<N;j++) { putc(' ', f); putc(' ', f); putc(L[it[j]][i], f); }
        putc('\n', f);
      }
    for (j=0;j<skip;j++) putc(' ', f);
    putc(' ', f);
    for (j=0;j<N;j++) { putc(' ', f); putc('-', f); putc('-', f); }
    putc('\n', f);
    
    /* Write body: */
    for(i=0;i<N;i++)
      { for (j=0;j<skip;j++) putc(L[it[i]][j], f);
        putc(' ', f);
        for (j=0;j<N;j++) 
          { fprintf(f, " %2d", FormatDist(d[it[i]][it[j]])); }
        putc('\n', f);
      }
    fflush(f);
  }
  
int FormatDist(double d)
  {
    return (int)(99.4999*(d<1?d:1.0) + 0.5);
  }
    
void WritePicture(FILE *f, double **d, long *it, long N)
  {
    long i, j;
    double dMax;
    pixel p[NUM_COLORS];
    /* Compute color palette: */
    SelectColors(p, NUM_COLORS);
    /* Find maximum distance: */
    dMax = 0.0e-16;
    for (i=0; i<N; i++)
      for (j=0; j<N; j++)
        { double dd= d[i][j];
          if (dd > dMax) dMax = dd;
        }
    /* Write PPM header: */
    fprintf(f, "P6\n%ld %ld\n255\n", N, N);
    /* Write distances as color values: */
    for(i=0;i<N;i++)
      { for (j=0;j<N;j++) 
          { 
            double dd = d[it[i]][it[j]];
            long ip = (long)((dd/dMax)*(NUM_COLORS-0.000001));
            long ipr = (ip < 0 ? 0 : (ip >= NUM_COLORS ? NUM_COLORS-1 : ip));
            pixel pi = p[ipr];
            putc(pi.R, f); putc(pi.G, f), putc(pi.B, f);
          }
      }
    fflush(f);
  }
  
#define Gamma  (2.0)
#define LogGamma (0.69314718055994530942)
  
void SelectColors(pixel *p, long M)
 { 
   long ip;
   for (ip=0; ip<M; ip++)
     { 
       pixel pi;
       double y, R, G, B, tR, tB;
       double d = ((double)ip)/((double)(NUM_COLORS-1));
       d = (d<1?(d>0?d:0.0):1.0);
       y = ( d==0 ? 1 : (exp(d*LogGamma)-Gamma)/(1-Gamma) );

       tR = y*y*(3-2*y);
       tB = y;

       R = 0.03*y + 0.97*(1-cos(3*M_PI*tR))/2;
       B = (1-cos(7*M_PI*tB))/2;
       G = (y - 0.299*R - 0.114*B)/0.588;
       if ((G < 0) || (G > 1)) fprintf(stderr, "y = %f G = %f\n", y, G);
       pi.R = (int)(R*255.0 + 0.5);
       pi.G = (int)(G*255.0 + 0.5);
       pi.B = (int)(B*255.0 + 0.5);
       p[ip] = pi;
     }
 }
    
void GetOptions(int argc, char **argv, Options *op)
  { int i = 1;
    /* Parse option switches: */
    op->cluSort = FALSE;
    op->farSort = FALSE;
    op->delReSort = FALSE;
    op->verbose = FALSE;
    op->discrete = FALSE;
    op->geometric = FALSE;
    op->picName = NULL;
    op->matName = NULL;
    op->width = 0;
    op->skip = 0;
    op->numValues = -1;
    op->repeat = 2;
    while ((i<argc) && (argv[i][0] == '-'))
      { if (
            GetLongOption(argc, argv, &i, "-n",  "--numValues",   0, MAX_VALUES,  &(op->numValues)) ||
            GetLongOption(argc, argv, &i, "-w",  "--width",       1, MAX_WIDTH,   &(op->width)) ||
            GetLongOption(argc, argv, &i, "-s",  "--skip",        0, MAX_SKIP,    &(op->skip)) ||
            GetLongOption(argc, argv, &i, "-r",  "--repeat",      0, MAX_REPEAT,  &(op->repeat)) ||
            GetBoolOption(argc, argv, &i, "-d",  "--discrete",    &(op->discrete)) ||
            GetBoolOption(argc, argv, &i, "-g",  "--geometric",   &(op->geometric)) ||
            GetBoolOption(argc, argv, &i, "-cs", "--cluSort",     &(op->cluSort)) ||
            GetBoolOption(argc, argv, &i, "-fs", "--farSort",     &(op->farSort)) ||
            GetBoolOption(argc, argv, &i, "-ds", "--delReSort",   &(op->delReSort)) ||
            GetCharOption(argc, argv, &i, "-m",  "--showMatrix",  &(op->matName)) ||
            GetCharOption(argc, argv, &i, "-p",  "--showPicture", &(op->picName)) ||
            GetBoolOption(argc, argv, &i, "-v",  "--verbose",     &(op->verbose))
          )
          { /* OK */ }
        else
          { ParamError("invalid option", argv[i], i); }
      }

    if (i < argc) 
      { ParamError("extraneous argument", argv[i], i); }

    if ((op->numValues) < 0) 
      { ParamError("must specify --numValues", "", argc); }

    if (op->cluSort && op->farSort)
      { ParamError("--cluSort and --farSort are mutually exclusive", "", argc); }
  }
  
bool_t GetLongOption(
    int argc, char **argv,
    int *ip,
    char *brief,
    char *wordy,
    long min,
    long max,
    long *vp
  )
  {
    int i = (*ip);
    int res; long val; char dum = '\000';
    if ((i+1<argc) && ((strcmp(argv[i], brief)==0) || (strcmp(argv[i], wordy)==0)))
      { res = sscanf(argv[i+1], "%ld%c", &val, &dum); 
        if ((res < 1) || (dum != '\000')) { ParamError("bad value", argv[i], i); }
        if ((val < min) || (val > max)) { ParamError("value out of range", argv[i], i); }
        (*vp) = val;
        (*ip) += 2;
        return (TRUE);
      }
    else
      { return (FALSE); }
  }
     
bool_t GetCharOption(
    int argc, char **argv,
    int *ip,
    char *brief,
    char *wordy,
    char **vp
  )
  {
    int i = (*ip);
    if ((i+1<argc) && ((strcmp(argv[i], brief)==0) || (strcmp(argv[i], wordy)==0)))
      { (*vp) = CopyString(argv[i+1]); 
        (*ip) += 2;
        return (TRUE);
      }
    else
      { return (FALSE); }
  }
  
bool_t GetBoolOption(
    int argc, char **argv,
    int *ip,
    char *brief,
    char *wordy,
    bool_t *vp
  )
  {
    int i = (*ip);
    if ((i<argc) && ((strcmp(argv[i], brief)==0) || (strcmp(argv[i], wordy)==0)))
      { (*vp) = TRUE;
        (*ip) += 1;
        return (TRUE);
      }
    else
      { return (FALSE); }
  }
     
void ParamError(char *msg, char* key, long K)
  { 
    fprintf(stderr, "%s: ARGV[%ld] = %s %s\n", progName, K, key, msg);
    fprintf(stderr, "usage: %s \\\n", progName);
    fprintf(stderr, "  {-n|--numValues} NNN \\\n");
    fprintf(stderr, "  [ {-w|--width} NNN ] \\\n");
    fprintf(stderr, "  [ {-s|--skip} NNN ] \\\n");
    fprintf(stderr, "  [ {-r|--repeat} NNN ] \\\n");
    fprintf(stderr, "  [ {-d|--discrete} ] \\\n");
    fprintf(stderr, "  [ {-g|--geometric} ] \\\n");
    fprintf(stderr, "  [ {-fs|--farSort} | {-cs|--cluSort} ] \\\n");
    fprintf(stderr, "  [ {-ds|--delReSort} ] \\\n");
    fprintf(stderr, "  [ {-m|--showMatrix} FILE ] \\\n");
    fprintf(stderr, "  [ {-p|--showPicture} FILE ] \\\n");
    fprintf(stderr, "  [ {-v|--verbose} ] \\\n");
    fprintf(stderr, "  < infile  > outfile\n");
    fflush(stderr);
    exit(1);
  }
  
void LineError(char *msg, long R)
  {
    fprintf(stderr, "%s:%ld: %s\n", progName, R, msg);
    fflush(stderr);
    exit(1);
  }

void FileError(char *msg, char *name)
  {
    fprintf(stderr, "%s: %s %s\n", progName, msg, name);
    fflush(stderr);
    exit(1);
  }

void ElementError(char *msg, long R, long K)
  {
    fprintf(stderr, "%s:%ld:value[%ld]: %s\n", progName, R, K, msg);
    fflush(stderr);
    exit(1);
  }

void Error(char *msg)
  {
    fprintf(stderr, "%s: %s\n", progName, msg);
    fflush(stderr);
    exit(1);
  }

/****************************************************************************/
/* (C) Copyright 1992 Universidade Estadual de Campinas (UNICAMP)           */
/*                    Campinas, SP, Brazil                                  */
/*                                                                          */
/* Authors:                                                                 */
/*                                                                          */
/*   Jorge Stolfi        - CS Dept, UNICAMP <stolfi@dcc.unicamp.br>         */
/*                                                                          */
/* This file can be freely distributed, modified, and used for any          */
/*   non-commercial purpose, provided that this copyright and authorship    */
/*   notice be included in any copy or derived version of this file.        */
/*                                                                          */
/* DISCLAIMER: This software is offered ``as is'', without any guarantee    */
/*   as to fitness for any particular purpose.  Neither the copyright       */
/*   holder nor the authors or their employers can be held responsible for  */
/*   any damages that may result from its use.                              */
/****************************************************************************/
