#define PROG_NAME "bool"
#define PROG_DESC "boolean op (union, intersection, or difference) of two sorted files"
#define PROG_VERS "2023-02-08"

#define PROG_C_COPYRIGHT \
  "Copyright © 1992 by the State University of Campinas (UNICAMP)"

/* Last edited on 2023-02-08 21:37:28 by stolfi */ 

#define PROG_HELP \
  "  " PROG_NAME " \\\n" \
  "    [ -s | --signedBytes | -u | --unsignedBytes ] \\\n" \
  "    [ -v | --verbose ] \\\n" \
  "    { 1-2 | 2-1 | 1+2 | 1.2 | 1~2 } \\\n" \
  "    {FILE1} {FILE2} \\\n" \
  "    > {OUTFILE} \\\n" \
  "    " argparser_help_info_HELP ""

#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  "  This program computes a specified Boolean operation (intersection, union," \
  " difference, or symmetric difference) of  {FILE1}  and  {FILE2},  viewed as sets of lines.\n" \
  "\n" \
  " Either  {FILE1}  or  {FILE2}  (but not both) may be \"-\", signifying the standard input ({stdin})." \
  " The result is always written to the standard output ({stdout}).\n" \
  "\n" \
  " The Boolean operation is specified by the first argument that is not an option switch" \
  " according to the following code:\n" \
  "\n" \
  "    1-2  difference: only the lines of  {FILE1}  that do not occur in  {FILE2}.\n" \
  "\n" \
  "    2-1  reverse difference: only the lines of  {FILE2}  that do not occur in  {FILE1}.\n" \
  "\n" \
  "    1.2  intersection: only the lines that occur in both files.\n" \
  "\n" \
  "    1+2  union: all the lines that occur in  {FILE1}  or  {FILE2}  or both.\n" \
  "\n" \
  "    1~2  symmetric difference: only the lines that occur in either  {FILE1}  or" \
  " {FILE2}, but not in both.\n" \
  "\n" \
  " The specifications \"2.1\", \"2+1\", \"2~1\" are also accepted" \
  " as equivalent to \"1.2\", \"1+2\", and \"1~2\", respectively.\n" \
  "\n" \
  " In the input files, every line, including the last one," \
  " must be terminated by an end-of-line character (octal \\012, hex 0x0a).  An empty" \
  " (zero length) file is assumed to contain no lines; all other" \
  " files must end with newline.\n" \
  "\n" \
  " The lines in each input file" \
  " must be sorted in strictly increasing collating order. See the" \
  " options \"-signedBytes\" and \"-unsignedBytes\" for details.\n" \
  "\n" \
  "  This program is somewhat similar to {comm}(1). However, unlike the latter," \
  " {bool} properly handles characters in the range \\200--\\377," \
  " lines of arbitrary length, embedded {NUL}s (zero bytes)," \
  " and complains if the input files are not sorted.\n" \
  "\n" \
  " Moreover, {bool} always writes the output in a single column, with no extraneous" \
  " spaces or {TAB}s." \
  "\n" \
  "OPTIONS\n" \
  "  -v\n" \
  "  --verbose\n" \
  "     If this argument is present, {bool} prints the input and output line" \
  " counts to the standard error stream ({stderr})." \
  "\n" \
  "  -s\n" \
  "  --signedBytes\n" \
  "     If this argument is present, {bool} assumes that the sorting of the input lines" \
  " compared bytes as if they were 8-bit signed integers (like the {char} type of" \
  " the C language).  Namely, characers in the range \\200--\\377 (which include" \
  " all the ISO  accented letters) come before the plain ASCII characters" \
  " \\000--\\177 in alphabetical order.  This option may be useful for compatibility with some" \
  "versions of {sort}(1) or some {LOCALE}s.\n" \
  "\n" \
  "  -u\n" \
  "  --unsignedBytes\n" \
  "     This option is the opposite of \"-signedBytes\", and tells {bool} to assume" \
  " that the inputs are sorted assuming bytes are /unsigned/ 8-bit integers; namely" \
  " that bytes \\000--\\177 come before \\200--\\377 in alphabetical" \
  " order.  It is the default behavior." \
  "BUGS\n" \
  " Likely.\n" \
  "\n" \
  "CHALLENGES FOR THE BORED\n" \
  "  Add an option to assume sorting as specified by the current {LC-COLLATE} environment" \
  " variable.  (But please do NOT make it the default behavior!)." \
  "\n" \
  "DOCUMENTATION OPTIONS\n" \
  argparser_help_info_HELP_INFO "\n" \
  "\n" \
  "SEE ALSO\n" \
  "  comm(1), join(1), uniq(1).\n" \
  "\n" \
  "AUTHOR\n" \
  "  Created 1992-02-11 by Jorge Stolfi, IC-UNICAMP.\n" \
  "\n" \
  "MODIFICATION HISTORY\n" \
  "  1992-02-11 J.Stolfi: created.\n" \
  "  2023-02-08 J.Stolfi: use the {artgparser.h} and {bool.h} modules from {JSLIBS/libjs}.\n" \
  "  2023-02-08 J.Stolfi: replaced the \"bool.1\" manpage by internal {PROG_INFO} string.\n" \
  "  2023-02-08 J.Stolfi: replaced legacy C int types by {stdint.h} ones.\n" \
  "  2023-02-08 J.Stolfi: changed the sym diff op from \"#\" to \"~\" for shell compat.\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " PROG_C_COPYRIGHT ".\n"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdint.h>

#include <bool.h>
#include <argparser.h>

typedef uint8_t byte_t;

typedef struct options_t 
  { 
    char *inName[2];     /* Input file names */
    bool_t verbose;      /* {TRUE} to print line counts */
    bool_t signedBytes;  /* {TRUE} to use signed byte comparisons */
    bool_t only0;        /* {TRUE} to write lines that are only in file 0 */
    bool_t both;         /* {TRUE} to write lines that are in both files */
    bool_t only1;        /* {TRUE} to write lines that are only in file 1 */
  } options_t;
  /* Command line arguments. */

#ifdef MSDOS

#include <process.h> 
#define RMODE "rb"
#define WMODE "wb"

#else

#define RMODE "r"
#define WMODE "w"

#endif

/* INTERNAL PROTOTYPES */

int32_t main(int32_t argc, char *argv[]);

FILE *OpenInputFile(char *name);
 /* Opens file "name" for reading; bombs in case of failure.
   If "name" is "-", returns stdin. */
     
bool_t ReadLine(FILE *rd, char *name, byte_t **L, uint64_t *N, uint64_t *A, uint64_t *R);
  /*  Returns TRUE if "rd" (whose filename is "name") was at end-of-file.
    Otherwise reads the next line of "rd" into "*L", and stores its length 
    into "*N", and increments "*R".  Assumes the buffer "*L" is "*A" 
    bytes long; reallocates a new one if necessary. 
    Stops with error if the last line does not end with newline. */

void WriteLine(byte_t L[], uint64_t N, uint64_t *R);
  /* Writes bytes "L[0..N-1]" to stdout, followed by newline.
    Increments the line counter "(*R)". */
 
int32_t CompareLines(byte_t *aL, uint64_t aN, byte_t *bL, uint64_t bN, bool_t signedBytes);
  /* Lexicographic comparison of lines "aL[0..aN-1]" and "bL[0..bN-1]".
    Returns -1 for "<", 0 for "=", +1 for ">".  Uses signed byte comparisons
    iff "signedBytes" is TRUE.  */
  
options_t *GetOptions(int32_t argc, char *argv[]);
  /* Parses the command line options. */

void ErrorOutOfMemory(char *name);
void ErrorCannotOpen(char *name);
void ErrorOutOfOrder(char *name, uint64_t R);
void ErrorNoFinalEOL(char *name);
  /* Each of these procedures prints an error message and aborts the program. */

/* PROCEDURES */

int32_t main(int32_t argc, char *argv[])
  {
    options_t *o = GetOptions(argc, argv);

    /* Internal variables: */
    FILE *inFile[2];      /* Input files */
    byte_t *thisL[2];       /* Last line of each file */
    uint64_t thisN[2];    /* Used size of "thisL[i]". */
    uint64_t thisA[2];    /* Allocated size of "thisL[i]" */
    byte_t *prevL[2];       /* next-to-last line of each file */
    uint64_t prevN[2];    /* Used size of "prev". */
    uint64_t prevA[2];    /* Allocated size of "prevL[i]" */
    bool_t done[2];       /* TRUE if file[i] is exhausted. */
    uint64_t nIn[2];       /* Count input lines by file */

    /* Open files, allocate line buffers, and read first lines */
    for (int32_t i = 0; i < 2; i++)
      { done[i] = FALSE;
        nIn[i] = 0;
        thisA[i] = 250;
        thisL[i] = (byte_t *)malloc(thisA[i]*sizeof(byte_t));
        prevA[i] = 250;
        prevL[i] = (byte_t *)malloc(prevA[i]*sizeof(byte_t));
        inFile[i] = OpenInputFile(o->inName[i]);
        done[i] = ReadLine(
          inFile[i], o->inName[i], 
          &thisL[i], &thisN[i], &thisA[i], 
          &nIn[i]
        );
      }

    uint64_t nOut = 0; /* Counts output lines */
    /* Merge lines and output some as requested: */
    while(! (done[0] & done[1]))
      { int32_t cmp;
        /* Compare next lines, one from each file: */
        if(done[0]) 
          cmp = 1;
        else if (done[1]) 
          cmp = -1;
        else 
          cmp = CompareLines(
            thisL[0], thisN[0], 
            thisL[1], thisN[1], 
            o->signedBytes
          );
        /* Output desired lines: */
        switch (cmp) {
          case -1: if (o->only0) WriteLine(thisL[0], thisN[0], &nOut); break;
          case 00: if (o->both)  WriteLine(thisL[0], thisN[0], &nOut); break;
          case  1: if (o->only1) WriteLine(thisL[1], thisN[1], &nOut); break;
        }
        /* Read next lines: */
        for (int32_t i = 0; i < 2; i++)
          { /* Check if {inFile[i]} needs to be read: */
            if (((cmp <= 0) && (i == 0)) || ((cmp >= 0) && (i == 1)))
              { /* Must read next line from file {inFile[i]} */
                /* Swap buffers {prevL[i],thisL[i]} and their alloc sizes {prevA[i],thisA[i]}: */
                byte_t *tl; uint64_t tn;
                tl = prevL[i]; prevL[i] = thisL[i]; thisL[i] = tl; 
                tn = prevA[i]; prevA[i] = thisA[i]; thisA[i] = tn;
                /* Save the current line as the previous one: */
                prevN[i] = thisN[i];
                /* Since we need to read, we must not have hit end-of-file yet: */
                assert(!done[i]); 
                /* So read the nct line, expanding the buffer as needed: */
                done[i] = ReadLine(
                  inFile[i], o->inName[i], 
                  &thisL[i], &thisN[i], &thisA[i],
                  &nIn[i]
                );
                /* If we did read a new line, check the order: */
                if (!done[i]) 
                  { int32_t ord = CompareLines
                      ( prevL[i], prevN[i], 
                        thisL[i], thisN[i], 
                        o->signedBytes
                      );
                    if (ord >=  0) { ErrorOutOfOrder(o->inName[i], nIn[i]); }
                  }
              }
          }
      }

    /* Close files and print report: */
    for (int32_t i = 0; i < 2; i++)
      { if (inFile[i] != stdin) { fclose(inFile[i]); }
        if (o->verbose)
          { fprintf(stderr, "file %d: %12ld lines\n", i+1, nIn[i]); }
      }
    fclose(stdout);
    if (o->verbose)
      { fprintf(stderr, "output: %12ld lines\n", nOut); }
    return 0;
  }
  
FILE *OpenInputFile(char *name)
  { FILE *f;
    if (strcmp(name, "-") == 0) { return(stdin); }
    f = fopen(name, RMODE);
    if (f == NULL) { ErrorCannotOpen(name); }
    return f;
  }

bool_t ReadLine(FILE *rd, char *name, byte_t **L, uint64_t *N, uint64_t *A, uint64_t *R)
  { byte_t *curL = (*L);
    uint64_t curA = (*A);
    uint64_t curN = 0;
    /* Resut of {fgetc} must be signed and more than 8 bits to detect end-of-file: */
    int32_t c = fgetc(rd); 
    if (c == EOF) { (*N) = 0; return(TRUE); }
    while (c != '\n')
      { if (curN >= curA) 
          { /* Must reallocate the buffer: */
            uint64_t newA = 2*curA;
            curL = (byte_t *)realloc(curL, newA*sizeof(byte_t));
            if (curL == NULL) { ErrorOutOfMemory(name); }
            curA = newA;
          }
        curL[curN] = (byte_t)c;
        curN++;
        c = fgetc(rd);
        if (c == EOF) { ErrorNoFinalEOL(name); }
      }
    (*L) = curL;
    (*A) = curA;
    (*N) = curN;
    (*R)++;
    return(FALSE);
  } 

void WriteLine(byte_t L[], uint64_t N, uint64_t *R)
  { for (uint64_t i = 0; i < N; i++) putchar(L[i]);
    putchar('\n');
    (*R)++;
  }

int32_t CompareLines(byte_t *aL, uint64_t aN, byte_t *bL, uint64_t bN, bool_t signedBytes)
  { uint64_t n = ( aN < bN ? aN : bN);
    if (signedBytes)
      { for (uint64_t i = 0; i < n; i++)
          { signed char ac = (signed char)(aL[i]);
            signed char bc = (signed char)(bL[i]);
            if (ac < bc) return -1;
            if (ac > bc) return +1;
          }
      }
    else
      { for (uint64_t i = 0; i < n; i++)
          { byte_t ac = aL[i];
            byte_t bc = bL[i];
            if (ac < bc) return -1;
            if (ac > bc) return +1;
          }
      }
    if (aN < bN) return -1;
    if (aN > bN) return +1;
    return 0;
  }

options_t *GetOptions(int32_t argc, char *argv[])
  {
    /* Initialize argument parser: */
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);

    /* Allocate the command line argument record: */
    options_t *o = (options_t*)malloc(sizeof(options_t));

    /* Parse keyword parameters: */

    if 
      ( argparser_keyword_present(pp, "--verbose") ||
        argparser_keyword_present(pp, "-verbose") ||
        argparser_keyword_present(pp, "-v")
      )
      { o->verbose = TRUE; }
    else
      { o->verbose = FALSE; }

    if 
      ( argparser_keyword_present(pp, "--signedBytes") ||
        argparser_keyword_present(pp, "-signedBytes") ||
        argparser_keyword_present(pp, "-s")
      )
      { o->signedBytes = TRUE; }
    else if 
      ( argparser_keyword_present(pp, "--unsignedBytes") ||
        argparser_keyword_present(pp, "-unsignedBytes") ||
        argparser_keyword_present(pp, "-u")
      )
      { o->signedBytes = FALSE; }
    else
      { o->signedBytes = FALSE; }
      
    /* Skip the command line switches, if any: */
    argparser_skip_parsed(pp);
      
    char *op = argparser_get_next_non_keyword(pp);
    if (strcmp(op, "1-2") == 0)
      { o->only0 = TRUE;  o->both = FALSE; o->only1 = FALSE; }
    else if (strcmp(op, "2-1") == 0)
      { o->only0 = FALSE; o->both = FALSE; o->only1 = TRUE;  }
    else if ((strcmp(op, "1+2") == 0) || (strcmp(op, "2+1") == 0))
      { o->only0 = TRUE;  o->both = TRUE;  o->only1 = TRUE;  }
    else if ((strcmp(op, "1.2") == 0) || (strcmp(op, "2.1") == 0))
      { o->only0 = FALSE; o->both = TRUE;  o->only1 = FALSE; }
    else if ((strcmp(op, "1~2") == 0) || (strcmp(op, "2~1") == 0))
      { o->only0 = TRUE;  o->both = FALSE; o->only1 = TRUE;  }
    else
      { argparser_error(pp, "bad operation"); }
 
    /* Parse filenames: */
    o->inName[0] = argparser_get_next_non_keyword(pp);
    o->inName[1] = argparser_get_next_non_keyword(pp);
    if((strcmp(o->inName[0], "-") == 0) && (strcmp(o->inName[1], "-") == 0))
      { argparser_error(pp, "standard input can be used only once"); }
    
    /* Check for spurious parameters: */
    argparser_finish(pp);
    
    return o;
  }
  
void ErrorCannotOpen(char *name)
  {
    fprintf(stderr, "%s: cannot open file \"%s\"\n", PROG_NAME, name);
    fflush(stderr);
    exit(1);
  }

void ErrorOutOfOrder(char *name, uint64_t R)
  {
    fprintf(stderr, 
      "%s: ** lines %lu and %lu of file \"%s\" are out of order\n", 
      PROG_NAME, R-1, R, name
    );
    fflush(stderr);
    exit(1);
  }

void ErrorNoFinalEOL(char *name)
  {
    fprintf(stderr, 
      "%s: ** missing final newline in file \"%s\"\n", 
      PROG_NAME, name
    );
    fflush(stderr);
    exit(1);
  }


void ErrorOutOfMemory(char *name)
  {
    fprintf(stderr, 
      "%s: ** out of memory while reading a line of \"%s\"\n", 
      PROG_NAME, name
    );
    fflush(stderr);
    exit(1);
  }


