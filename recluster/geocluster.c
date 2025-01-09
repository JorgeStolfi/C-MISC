#define PROG_NAME "geocluster"
#define PROG_DESC "(re)groups vectors into clusters, geometrically"
#define PROG_VERS "2.0"

/* Copyright © 1997 by the State University of Campinas (UNICAMP). */
/* See the authorship, rights and warranty notices in the PROG_INFO below. */
/* Last edited on 2024-12-21 11:54:46 by stolfi */

#define PROG_HELP \
  "  " PROG_NAME " \\\n" \
  "    { -numValues} NNN \\\n" \
  "    [ {-width} NNN ] \\\n" \
  "    [ {-skip} NNN ] \\\n" \
  "    [ {-repeat} NNN ] \\\n" \
  "    [ {-discrete} ] \\\n" \
  "    [ {-geometric} ] \\\n" \
  "    [ {-farSort} | {-cluSort} ] \\\n" \
  "    [ {-delReSort} ] \\\n" \
  "    [ {-showMatrix} FILE ] \\\n" \
  "    [ {-showPicture} FILE ] \\\n" \
  "    [ {-verbose} ] \\\n" \
  "    " argparser_help_info_HELP " \\\n" \
  "    < {INFILE} \\\n" \
  "    > {OUTFILE}"

#define PROG_INFO_DESC \
  "  Reads from stdin a set of items, each a point of {N}-dimensional" \
  " vector space, and writes them to stdout, (re)grouped into a specified" \
  " number of clusters, based on their relative distances.\n" \
  "\n" \
  "  Each item (input or output) must be given in a separate line," \
  " containing a vector of {N} real numbers, possibly preceded and/or" \
  " followed by other information.\n" \
  "\n" \
  "  The program tries to find a partition that minimizes the" \
  " square of the distance from each item to centroid of its" \
  " respective cluster, summed over all items." \
  "\n" \
  "  The program uses an iterative, non-deterministic algorithm that" \
  " improves a random initial partition until a local optimum " \
  " is found.  It can also take an initial partition of the" \
  " items into clusters, as a starting solution.\n" \
  "\n" \
  "  In the input and output files, clusters are separated by " \
  " blank lines."

#define PROG_INFO_OPTS \
  "  -numValues  {NUM}\n" \
  "    The number of data coordinates for each item.  This" \
  " parameter is required.\n" \
  "\n" \
  "  -skip  {NUM}\n" \
  "    Skip the first {NUM} bytes of every line, before parsing the" \
  "  first data value.  The default is zero.\n" \
  "\n" \
  "  -width  {NUM}\n" \
  "    Specifies that the data values are given in fixed format," \
  " as {N} adjacent fields, each {NUM} bytes wide.  If this parameter" \
  " is not specified, {" PROG_NAME "} assumes that the values are given" \
  " in free format.  In that case each field must be terminated by" \
  " a space, comma, or newline, and may be surrounded by extra spaces.  In" \
  " either case, fields consisting entirely of periods or dashes are" \
  " taken to be zero.\n" \
  " \n" \
  "  -discrete\n" \
  "    If this option is given, {" PROG_NAME "} assumes that the numbers in" \
  " each line are observed counts in some sampling experiment, and" \
  " converts them to probabilities taking into account sampling error.\n" \
  "\n" \
  "  -geometric\n" \
  "    By default, {" PROG_NAME "} compares the item data" \
  " vectors with the Euclidean metric.  If this option is given," \
  " however, {" PROG_NAME "} assumes that successive items are associated" \
  " with equally spaced positions along a unidimensional path, and uses" \
  " the \"earth-movers\" distance to compare their data vectors.\n" \
  "\n" \
  "  -farSort\n" \
  "    Specifies that the first sorting pass should be done with" \
  " the {Farthest-Point} alternate sorting heuristic.\n" \
  "\n" \
  "  -cluSort\n" \
  "    Specifies that the first sorting pass should be done with" \
  " the {Cluster-Tree} alternate sorting heuristic.\n" \
  "\n" \
  "  -delReSort\n" \
  "    Specifies that the second and additional sorting passes should" \
  " be done with the {Delete-And-Insert} alternate re-sorting heuristic.\n" \
  "\n" \
  "  -showMatrix {FILE}\n" \
  "    Print the item distance matrix to the named file" \
  " (or to standard output if {FILE} is \"-\"). The distances" \
  " are scaled to the range [0..99.5] and then" \
  " rounded to the nearest integer.\n" \
  "\n" \
  "  -showPicture {FILE}\n" \
  "    Writes the distance matrix to the named file" \
  " (or standard output if {FILE} is \"-\"), as a color PPM" \
  " image.  Each element becomes a colorful pixel, whose" \
  " brightness decreases with increasing distance (white for distance 0," \
  " black for distance 1).\n" \
  "\n" \
  "  -repeat {COUNT}\n" \
  "    Number of sorting passes (default 2). The first pass finds the" \
  " two items with largest distance, places them at opposite ends of the" \
  " list, and sorts all other items based on the ratio of their distances" \
  " to those two. The second pass looks at each item in turn, and tries to" \
  " move it along the list so as to decrease the path length.  If" \
  " {COUNT} is greater than 2, the second pass is performed" \
  " {COUNT-1} times. If {COUNT} is 1, performs only the first pass." \
  " If {COUNT} is zero, the items are not sorted and output is supressed." \
  " (This choice is useful for generating a distance matrix" \
  " or picture of already sorted items.)\n" \
  "\n" \
  "  -verbose\n" \
  "    Print various progress messages to standard error."

#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  PROG_INFO_DESC "\n" \
  "\n" \
  "OPTIONS\n" \
  PROG_INFO_OPTS "\n" \
  "\n" \
  "DOCUMENTATION OPTIONS\n" \
  argparser_help_info_HELP_INFO "\n" \
  "\n" \
  "SEE ALSO\n" \
  "  " ??? "(1).\n" \
  "\n" \
  "CHALLENGES FOR THE BORED\n" \
  "  Add other distribution metrics.\n" \
  "  Improve the sorting algorithm.\n" \
  "  Add dendrogram output option.\n" \
  "\n" \
  "AUTHOR\n" \
  "  Created as {sort-distr.c} on 1997-11-23 (or earlier) by Jorge Stolfi, IC-UNICAMP." \
  "  Originally motivated by Voynich Manuscript word-frequency analysis.\n" \
  "\n" \
  "MODIFICATION HISTORY\n" \
  "   1997-12-07  Changed -binSort heuristic to break the path preferably" \
  "               at long edges.\n" \
  "\n" \
  "   1997-11-30  Changed -binSort heuristic to break the path preferably" \
  "               at long edges.\n" \
  "\n" \
  "   1997-11-29  -showPicture now precomputes the color table so as" \
  "               to use only 216 different colors.\n" \
  "\n" \
  "   1997-11-27  Improved (?) sorting algorithm, with DiaSortItems" \
  "               or FarSortItems followed by UniReSortItems or DelReSortItems." \
  "               Changed -showPicture and -showMatrix to take a file name.\n" \
  "\n" \
  "   1997-11-26  Added -showPicture option." \
  "               Also improved the -showMatrix format.\n" \
  "\n" \
  "   2003-09-25 Edited the documentation?\n" \
  "\n" \
  "   2007-08-15 Revamped code to my current standards, tried to finish it " \
  "              and make it work.\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  Copyright © 1997 by the State University of Campinas (UNICAMP).\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS

#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>

#include <bool.h>
#include <argparser.h>

#include <jclbasic.h>
#include <jcldist.h>
#include <jclerror.h>
#include <jclgeo.h>
#include <jclimage.h>
#include <jcloptions.h>
#include <jclparse.h>
#include <jclsort.h>
#include <jcltree.h>

/* COMMAND-LINE OPTIONS */

typedef struct options_t
  {
    long numValues;    /* Number of values in each distribution. */
    long width;        /* Width of each value, or 0 for free format. */
    long skip;         /* Number of bytes to skip before first value. */
    long repeat;       /* Number or times to apply "UniReSortItems" */
    bool_t discrete;   /* TRUE if inputs are counts, FALSE if they are probs. */
    bool_t geometric;  /* TRUE for earth mover's distance, FALSE for Euclidean. */
    bool_t cluSort;    /* TRUE to begin with CluSortItems. */
    bool_t farSort;    /* TRUE to begin with FarSortItems. */
    bool_t delReSort;  /* TRUE to use DelReSortItems. */
    char *matName;     /* File name for distance matrix, or NULL if none. */
    char *picName;     /* File name for PPM distance image; or NULL if none. */
    bool_t verbose;    /* TRUE to mumble while working. */
  } options_t;

/* INTERNAL PROTOTYPES */

int main(int argc,char** argv);

options_t *parse_options(int argc, char **argv);
  /* Parses the command line arguments and packs them as an {options_t}. */

void AllocItemEntry(long i, byte ***L, long **N, long *A);
  /* Makes sure elements "(*L)[i]" and "(*N)[i]" exist,
    reallocating "*L" and "*N" if necessary.
    Assumes their current allocated size is "*A". */
 
/* IMPLEMENTATIONS */

int main(int argc, char **argv)
  {
    /* Command line options: */
    options_t *o = parse_options(argc, argv);

    /* Read buffer: */
    byte *bufL;        /* Last item of input file. */
    long bufN;         /* Used size of "bufL". */
    long bufA;         /* Allocated size of "bufL". */

    /* The input data: */
    long nItems;       /* Number of items stored so far. */
    byte **L;          /* "L[0..nItems-1i]" are pointers to the items read so far. */
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
    minItemLength = o->skip + (o->width == 0 ? (2*o->numValues - 1) : (o->numValues * o->width));

    /* Initialize: */

    bufA = 0;
    nItems = 0;
    aItems = 0;

    /* Read data: */
    while(! ReadItem(&bufL, &bufN, &bufA, nItems+1))
      { if (bufN < minItemLength) 
          { LineError("line too short", "", nItems+1); }
        else
          { AllocItemEntry(nItems, &L, &N, &aItems);
            L[nItems] = CopyItem(bufL, bufN);
            N[nItems] = bufN;
            nItems++;
          }
      }
    if (o->verbose) { fprintf(stderr, "%6ld items read.\n", nItems); }

    /* Extract the numerical vectors: */
    if (o->verbose) { fprintf(stderr, "parsing distributions...\n"); }
    v = (double**)(Alloc(nItems*sizeof(double*)));
    for (i=0; i<nItems; i++)
      { v[i] = (double*)(Alloc(o->numValues*sizeof(double)));
        if (o->width == 0)
          { ParseFreeVector(&(L[i][o->skip]), N[i]-o->skip, v[i], o->numValues, "", i+1); }
        else
          { ParseFixedVector(&(L[i][o->skip]), N[i]-o->skip, o->width, v[i], o->numValues, "", i+1); }
       
        if (o->discrete) 
          { EstimateProbs(v[i], o->numValues, i+1); }
        else
          { NormalizeProbs(v[i], o->numValues, i+1); }
      }

    /* Compute distance matrix: */
    /* We should use a triangular matrix, but it is simpler to keep it square. */
    if (o->verbose) { fprintf(stderr, "computing distance matrix...\n"); }
    d = (double**)(Alloc(nItems*sizeof(double*)));
    for (i=0; i<nItems; i++)
      { d[i] = (double*)(Alloc(nItems*sizeof(double)));
        for (j=0; j<i; j++)
          { double dd;
            if (o->geometric) 
              { dd = EarthMoverDist(v[i], v[j], o->numValues); }
            else
              { dd = L2ProbDist(v[i], v[j], o->numValues); }
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
    if (o->repeat > 0)
      { /* Sort and print items: */
        long r;
        if (o->verbose) { fprintf(stderr, "initial sorting...\n"); }
        if (o->farSort) 
          { FarSortItems(it, nItems, d); }
        else if (o->cluSort) 
          { CluSortItems(it, nItems, d); }
        else
          { DiaSortItems(it, nItems, d); }
        if (o->verbose) { fprintf(stderr, "refining order... "); }
        for (r = 2; r <= o->repeat; r++) 
          { if (o->verbose) { fprintf(stderr, "+"); }
            if (o->delReSort) 
              { DelReSortItems(it, nItems, d); }
            else
              { UniReSortItems(it, nItems, d); }
          } 
        if (o->verbose) { fprintf(stderr, "\n"); }
        for (k=0; k<nItems; k++) 
          { long i = it[k];
            WriteItem(L[i], N[i]);
          }
      }
    
    if (o->matName != NULL) 
      { /* Print distance matrix: */
        if (o->verbose) { fprintf(stderr, "writing matrix...\n"); }

        /* Compute color palette: */
        SelectColors(map, NUM_COLORS);
        /* Find maximum distance: */
        dMax = 0.0e-16;
        for (i=0; i<N; i++)
          for (j=0; j<N; j++)
            { double dd= d[i][j];
              if (dd > dMax) dMax = dd;
            }

        WriteMatrix(OpenWrite(o->matName), d, it, nItems, L, o->skip);
      }
    
    if (o->picName != NULL) 
      { /* Print PBM image of distance matrix: */
        if (o->verbose) { fprintf(stderr, "writing picture...\n"); }
        WritePicture(OpenWrite(o->picName), d, it, nItems);
      }
      
    fclose(stdout);
    return 0;
  }

void AllocItemEntry(long i, byte ***L, long **N, long *A)
  {
    if (i >= (*A))
      { 
        if (i >= MAX_ITEMS)
          { LineError("too many lines", "", i+1); }
        else 
          { byte **curL = (*L);
            long *curN = (*N);
            long curA = (*A);
            long newA = (curA == 0 ? 256 : 2*curA);
            byte **newL = (byte **)Alloc(newA*sizeof(byte*));
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

options_t *parse_options(int argc, char **argv)
  {
    /* Initialize argument parser: */
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);
    
    /* Allocate the command line argument record: */
    options_t *o = (options_t *)malloc(sizeof(options_t)); 
    
    /* Parse keyword parameters: */
    
    argparser_get_keyword(pp, "-numValues");
    o->numValues = argparser_get_next_int(pp, 0, MAX_VALUES);

    if (argparser_keyword_present(pp, "-width"))
      { o->width = argparser_get_next_int(pp, 1, MAX_WIDTH); }
    else
      { o->width = 0; }

    if (argparser_keyword_present(pp, "-skip"))
      { o->skip = argparser_get_next_int(pp, 0, MAX_SKIP); }
    else
      { o->skip = 0; }

    if (argparser_keyword_present(pp, "-repeat"))
      { o->repeat = argparser_get_next_int(pp, 0, MAX_REPEAT); }
    else
      { o->repeat = 2; }

    o->discrete = argparser_keyword_present(pp, "-discrete");

    o->geometric = argparser_keyword_present(pp, "-geometric");

    o->cluSort = argparser_keyword_present(pp, "-cluSort");

    o->farSort = argparser_keyword_present(pp, "-farSort");

    o->delReSort = argparser_keyword_present(pp, "-delReSort");
 
    if (argparser_keyword_present(pp, "-showMatrix"))
      { o->matName = argparser_get_next(pp); }
    else
      { o->matName = NULL; }

    if (argparser_keyword_present(pp, "-showPicture"))
      { o->picName = argparser_get_next(pp); }
    else
      { o->picName = NULL; }

    o->verbose = argparser_keyword_present(pp, "-verbose");

    /* Parse positional arguments: */
    argparser_skip_parsed(pp);

    /* Check for spurious arguments: */
    argparser_finish(pp);
    
    if (o->cluSort && o->farSort)
      { argparser_error(pp, "\"-cluSort\" and \"-farSort\" are mutually exclusive"); }

    return o;
  }
