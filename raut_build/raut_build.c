#define PROG_NAME "dicio_build"
#define PROG_DESC "build a finite automaton from a set of strings"
#define PROG_VERS "4.0.0-1"

#define PROG_C_COPYRIGHT "Copyright © 2009 Universidade Estadual de Campinas (UNICAMP)."

/* Last edited on 2009-10-31 00:21:29 by stolfi */

#define PROG_HELP \
  "  " PROG_NAME " \\\n" \
  "    {ACTION}.. \\\n" \
  "  where {ACTION} is any of \n" \
  "    -load {DMP_FILE} \n" \
  "    -dump {DMP_FILE} \n" \
  "    -list {STRINGS_FILE} \n" \
  "    -check {STRINGS_FILE} {BAD_FILE} \n" \
  "    -add {STRINGS_FILE} \n" \
  "    -remove {STRINGS_FILE} \n" \
  "    -modify {STRINGS_FILE} \n" \
  "    -expand {NUM_NODES} \n" \
  "    -stats {TXT_FILE} \n" \
  "    -debug {TXT_FILE} \n" \
  "    -comment {STRING} \n" \
  "    -description {DESCR_FILE} \n" \
  "    -report {NUM_STRINGS} \n" \
  "    -crunch {NUM_STRINGS} \n" \
  "    -plot {PLT_FILE} {NUM_STRINGS} \n" \
  "    -forgive"

#define MAX_STRING_SYMBOLS 255
  /* Maximum length (symbols) of a language string, in internal encoding. */

#define MAX_LINE_BYTES 600
  /* Maximum length of a line (bytes) in string files, excluding EOL.
    Allow for leading '+'/'-' code, some leading and trailing whitespace,
    and multiple embedded blanks. Must be at least {MAX_STRING_SYMBOLS + 2} */

#define FIRST_EXPAND 100000
  /* Argument of first implicit expand when adding to an empty automaton. */

#define DEFAULT_REPORT_STEP 10000
  /* Print a report report message every this many input strings. */

#define DEFAULT_CRUNCH_STEP 10000
  /* Print a report report message every this many input strings. */

#define MAX_STRING_STEP (~ 0LLU)
  /* Max argument for the "-report", "-crunch", and "-plot" options. */

#define MAX_EXPAND_ARG 1073741823
  /* Maximum argument of "-expand" command.
    Should not exceed {2^{rdag_nn_MAX}-1}. */

#define STRINGIFY(x) #x

#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  "  Maintains finite-state automata representing finite sets of " \
  " strings.\n" \
  "\n" \
  "  In normal usage, this program reads or builds an initial" \
  " automaton, performs some modifications on it, prints some" \
  " reports and auxiliary listings, and optionally writes" \
  " out the resulting automaton.\n" \
  "\n" \
  "  The {ACTION}s specified on the command line are executed" \
  " in the order given. Some actions have parameters.  Each {OPTION}" \
  " applies to all subsequent actions.\n" \
  "\n" \
  "ACTIONS\n" \
  PROG_INFO_ACTIONS1 "\n" \
  "\n" \
  PROG_INFO_ACTIONS2 "\n" \
  "\n" \
  "DOCUMENTATION OPTIONS\n" \
  argparser_help_info_HELP_INFO "\n" \
  "\n" \
  "BUGS\n" \
  "  This is a free program, so there is no extra charge for the bugs.\n" \
  "\n" \
  "SEE ALSO\n" \
  "  cat(1), sort(1), uniq(1).\n" \
  "\n" \
  "AUTHOR\n" \
  "  Created on 2009-10-27 by Jorge Stolfi, UNICAMP.\n" \
  "  Based on MaintainAutomaton.m3 by J.Stolfi, C.Lucchesi and T.Kowaltowski (early 1992).\n" \
  "\n" \
  "MODIFICATION HISTORY\n" \
  PROG_INFO_HISTORY "\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " PROG_C_COPYRIGHT ".\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS
  
#define PROG_INFO_ACTIONS1 \
  "  -load {DMP_FILE}\n" \
  "    Read the initial automaton from {DMP_FILE}, which should be in" \
  " the \".dmp\" format (see {raut.h}).  If {DMP_FILE} is \"-\", reads" \
  " from {stdin}.\n" \
  "\n" \
  "  -dump {DMP_FILE}\n" \
  "    Write the automaton to {DMP_FILE}, in the \".dmp\" format.  If {DMP_FILE}" \
  " is \"-\", writes to {stdout}." \
  "\n" \
  "  -add {STRINGS_FILE}\n" \
  "    Read strings from {STRINGS_FILE} and add them to the" \
  " automaton's language.  If {STRINGS_FILE} is \"-\", reads from" \
  " {stdin}.\n" \
  "\n" \
  "  -remove {STRINGS_FILE}\n" \
  "    Read strings from {STRINGS_FILE} and remove them from the" \
  " automaton's language.  If {STRINGS_FILE} is \"-\", reads from" \
  " {stdin}.\n" \
  "\n" \
  "  -modify {STRINGS_FILE}" \
  "     Reads strings from {STRINGS_FILE}, each preceded by '+' or '-' and optional" \
  " whitespace. Adds the '+' strings and removes the '-' strings from" \
  " the automaton's language.  If {STRINGS_FILE} is \"-\", reads from" \
  " {stdin}.\n" \
  "\n" \
  "  -check {STRINGS_FILE} {BAD_FILE}\n" \
  "    Read strings from {STRINGS_FILE} and write to {BAD_FILE} those that are" \
  " not accepted by the automaton.  If {STRINGS_FILE} is \"-\", reads the" \
  " query words from {stdin}.  If {BAD_FILE} is \"-\", writes the bad words" \
  " to {stdout}.\n" \
  "\n" \
  "  -list {STRINGS_FILE}\n" \
  "    Write to {STRINGS_FILE} all words accepted by the automaton.  If {STRINGS_FILE}" \
  " is \"-\", writes them to {stdout}.\n" \
  "\n" \
  "  -expand {NUM_NODES}\n" \
  "    Expands the automaton to allow for {NUM_NODES} proper nodes" \
  " (automaton transitions).  This action only affects the" \
  " running time, since the automaton is automatically" \
  " expanded (even beyond the {NUM_NODES}) as needed.\n" \
  "\n" \
  "  -forgive" \
  "    In subsequent \"-add\", \"-remove\", and \"-modify\" actions, silently ignore" \
  " strings that are redundant (strings to" \
  " add are already in the automaton, or strings to remove which aren't there).\n" \
  "\n" \
  "  -stats {TXT_FILE}\n" \
  "    Writes statistics about the current automaton to the file" \
  " {TXT_FILE}.  If {TXT_FILE} is \"-\", writes them to {stderr}.\n" \
  "\n" \
  "  -debug {TXT_FILE}\n" \
  "    Write to {TXT_FILE} a readable description of the automaton.  If {TXT_FILE}" \
  " is \"-\", writes to {stderr}.\n" \
  "\n" \
  "  -comment {STRING}\n" \
  "    Include {STRING} in the next \"-dump\" file, as a" \
  " comment.  Overrides any preceding \"-comment\" or \"-description\" options.\n" \
  "\n" \
  "  -description {DOC_FILE}\n" \
  "    Copy the content of {DOC_FILE} into the next \"-dump\" file, as" \
  " a comment.  Overrides any preceding \"-comment\" or \"-description\" options.\n" \
  "\n" \
  "  -report {NUM_STRINGS}\n" \
  "    During all subsequent \"-add\", \"-remove\", \"-modify\", and \"-check\" actions," \
  " print a progress report line after processing every {NUM_STRINGS} input strings." \
  " The default is to report every " STRINGIFY(DEFAULT_REPORT_STEP) " strings.  If" \
  " {NUM_STRINGS} is zero, reporting is turned off.\n" \
  "\n" \
  "  -crunch {NUM_STRINGS}\n" \
  "    During all subsequent \"-add\", \"-remove\", and \"-modify\" actions," \
  " crunch the automaton after processing every {NUM_STRINGS} input strings." \
  " The default is to crunch every " STRINGIFY(DEFAULT_CRUNCH_STEP) " strings.  If" \
  " {NUM_STRINGS} is zero, crunching is never done.\n" \
  "\n" \
  "  -plot {PLT_FILE} {NUM_STRINGS}\n" \
  "    During the next \"-add\", \"-remove\", or \"-modify\" action," \
  " write a line to {PLT_FILE} after processing every {NUM_STRINGS} input" \
  " strings.  Each line will show:\n" \
  "      * the number of strings read so far,\n" \
  "      * the number of input symbols in all those strings,\n" \
  "      * the maximum allocated node number,\n" \
  "      * the maximum node number present in the dag,\n" \
  "      * the root state (accept bit and dag node), and\n" \
  "      * the number of nodes reachable from the root node.\n" \
  "    If {NUM_STRINGS} is zero, this command is ignored.  Plotting only" \
  " applies to the next one of the above commands."
  
#define PROG_INFO_ACTIONS2 \
  "  The files read by \"-add\", \"-remove\", \"-modify\", and \"-check\" must" \
  " contain one string per line.  Valid word characters are all printing" \
  " ISO-Latin-1 characters ('\041' to '\176' and '\241' to '\377') plus" \
  " the plain ASCII space ('\040').  Leading and trailing blanks are" \
  " discarded; any embedded whitespace is squeezed to a single space and treated as" \
  " a letter.  In particular, a blank line is treated as the empty" \
  " word.  For safety, input lines are currently limited" \
  " to " STRINGIFY(MAX_LINE_BYTES) " characters, including" \
  " all whitespace and '+'/'-' op codes; and to " STRINGIFY(MAX_STRING_SYMBOLS) " significant" \
  " symbols, including any embedded blanks.\n" \
  "\n" \
  "  The \"-list\" command writes the words in the same format, without" \
  " leading or trailing blanks, in lexicographic order.\n" \
  "\n" \
  "  Usually the first action is \"-load\", \"add\", or \"-modify\"; until" \
  " then the automaton is empty (accepts no strings).  The first \"-add\"" \
  " or \"-modify\" action on an empty automaton will expand it to a" \
  " large initial size, as if \"-expand " STRINGIFY(FIRST_EXPAND) "\" had" \
  " been performed first." \

#define PROG_INFO_HISTORY \
  "  ca. jan/1992 MaintainAutomaton.m3 written by J.Stolfi at DCC-IMECC-UNICAMP.\n" \
  "\n" \
  "  1992-04-13: Modified for compatibility with version 2.04 of Modula-3 and the" \
  " new interfaces [TK].\n" \
  "\n" \
  "  1992-11-02: \"-noRedundant\" option added [TK].\n" \
  "\n" \
  "  1995-08-06: Converted to Modula-3 release 3.5.3.  Substantial" \
  " rewrite of command line processing.  Fixed comments to better" \
  " match the code, etc. [JS]\n" \
  "\n" \
  "  1995-09-24: Removed the \"-build\" command (use \"-add\" without \"-load\"). Changed" \
  " command syntax to require {FILENAME} argument for the" \
  " commands \"-check\", \"-add\", \"-remove\", \"-addRemove\" (which formerly" \
  " read from stdin) and \"-spell\" (which formerly wrote to stdout). [JS]\n" \
  "\n" \
  "  1996-11-01: Added a quick fix to the word length histogram code so" \
  " that words longer than {MaxWordLength} are counted too.  [JS]. This" \
  " was the last version of MaintainAutomaton.m3 to be effectively used (dicio-v353-1).\n" \
  "\n" \
  "  1997-01-19: Removed the historical compatibility" \
  " hacks (\"-dumpin\", \"-dumpout\", \"-mess\", etc.).  Upgraded" \
  " to use the new {libm3reduced-2}. [JS] This was the last edit" \
  " in the Modula-3 version. Apparently, this new version of" \
  " MaintainAutomaton.m3 (dicio-v353-2) was never used.\n" \
  "\n" \
  "  2009-10-27 Converted from Modula-3 to C and to the" \
  " new dag structure {libdicio.a}.  Renamed some actions" \
  " and options: \"-spell\"=\"-list\", \"-addRemove\"=\"-modify\"," \
  " \"-noRedundant\"=\"-forgive\", \"-print\"=\"debug\"," \
  " \"-commentFile\"=\"-description\".  Assigned a new version" \
  " number, '4.0.0-1'. [JS]\n" \
  "\n" \
  "  2009-10-30 Added \"-plot\" and \"-crunch\" options. [JS]" \

#define int64_NONE (~ 0LL)
  /* A 64-bit integer value used to represent "no number". */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <values.h>
#include <assert.h>
#include <string.h>

#include <bool.h>
#include <jsfile.h>
#include <affirm.h>
#include <argparser.h>
#include <vec.h>
#include <jsstring.h>
#include <jstime.h>

#include <rdag.h>
#include <raut.h>
#include <raut_io.h>
#include <rdag_io.h>

typedef struct command_t
  { char *op;    /* Command ("-load", "-add", etc.). */
    char *str1;  /* First filename/string argument of command, or NULL. */
    char *str2;  /* Second filename/string argument of command, or NULL. */
    int64_t num; /* Numeric argument of command, or NULL. */
  } command_t;
  /* A parsed command-line action or option, with its arguments if any. */

vec_typedef(command_vec_t,command_vec,command_t);

typedef struct prog_state_t
  { command_vec_t cmd;    /* Command/option list. */
    int icmd;             /* Commands {cmd[0..icmd-1]} have been finished. */
    raut_t *A;            /* The current automaton. */
    char *cmt;            /* Comment string for next "-dump" command, or NULL. */
    bool_t forgive;       /* TRUE to silently ignore redundant add/remove strings. */
    uint64_t report_step; /* Report progress every this many strings. */
    uint64_t crunch_step; /* Crunch the automaton every this many strings. */
    uint64_t plot_step;   /* Write a plot data line every this many strings. */
    char *plot_fname;     /* Name of plot data file. */
  } prog_state_t;
  /* The program's atate. */

typedef enum
  { update_op_ADD, /* Add strings. */
    update_op_SUB, /* Remove strings. */
    update_op_MOD, /* Add/remove strings, depending on first char. */
    update_op_CHK  /* Check strings. */
  } update_op_t;
  /* Automaton modification mode. */

typedef struct clocks_t { double rtime; double utime; double stime; } clocks_t;

int main (int argc, char **argv);
void do_show_args(FILE *wr, command_vec_t *cmd);
void do_load(prog_state_t *st, char *rd_fname);
void do_dump(prog_state_t *st, char *wr_fname);
void do_process_strings
  ( prog_state_t *st,
    char *rd_fname,
    char *wr_fname,
    update_op_t op
  );
void do_list(prog_state_t *st, char *wr_fname);
void do_debug(prog_state_t *st, char *wr_fname);

command_vec_t parse_commands(int argc, char **argv);
char *parse_filename(argparser_t *pp);
void print_doc(FILE *wr, char *cmt, char *prefix);

void clocks_save(clocks_t *c);
  /* Stores the current user, system, and real time readings into {c}. */

char *clocks_diff_string(clocks_t *c);
  /* Returns a string with the differences between the current
    user, system, and real time readings and those saved in {c}. */

void command_start_message(char *msg);
  /* Prints "started {msg}" to {stderr}, with separating line. */

void command_finish_message(char *msg);
  /* Prints "finished {msg}" to {stderr}, with CPU time and separating line. */

void dashed_message(char *msg1, char *msg2, char *msg3);
  /* Prints "--- {msg1} {msg2} {msg3} ---------\n" to {stderr}.
    Any arg is omitted if NULL or empty. */

void data_error_message(char *fname, uint32_t nlin, char *msg, int nbuf, char buf[]);
  /* Prints to {stderr} a message {msg} about an error on or after line {nlin} of the
    data file {fname}, whose contents (possibly truncated) is {buf[0..nbuf-1]}. */

FILE *do_open_read(char *fname, char *what, FILE *def_file, char *def_fname);
  /* Opens {fname} for reading; bombs out if the operation failed for
    any reason. However, if {fname} is NULL or empty, returns
    {def_file} which is assumed to be already open for reading. In
    either case, prints a warning "reading {what} from {name}\n" where {name}
    is either {fname} or {def_fname}. */

FILE *do_open_write(char *fname, char *what, FILE *def_file, char *def_fname);
  /* Opens {fname} for writing; bombs out ofthe operation failed for
    any reason. However, if {fname} is NULL or empty, returns
    {def_file} which is assumed to be already open for writing. In
    either case, prints a warning "writing {what} to {name}\n" where {name} is
    either {fname} or {def_fname}. */

void do_close_file(FILE *f);
  /* If {f} is NULL or {stdin}, does nothing.  If {f} is {stdout} or {stderr},
    flushes any pendng output of {f}, but leaves it open.  Otherwise closes
    the file {f}. */

void write_plot_line(FILE *wr, uint64_t nstrings, uint64_t nsymbols, raut_t *A);
  /* Writes to {wr} a line containing {nstrings}, {nsmbols}, and some vital statistics
    of the automaton {A}. */

bool_t get_string_from_line
  ( char * rd_fname, 
    uint64_t line_count, 
    uint32_t nbuf, 
    char buf[], 
    update_op_t op, 
    uint32_t *nstr, 
    rdag_symbol_t str[],
    update_op_t *str_op
  );
  /* Extracts an input symbol string argument from a data line
    {buf[0..nbuf-1]}, and decides the operation {*str_op} to be
    performed on it. The number of symbols is placed in {*nstr} and
    the symbols in {str[0..*nstr-1]}.
    
    Assumes that the data line was read from a file that is argument
    to "-add", "-remove", "modify", or "-check". Namely, if {op} is
    {update_op_MOD}, the line must beging with '+' or '-'; the
    procedure discards that character and sets {*str_op} to {update_op_ADD}
    or {update_op_SUB}. In all other cases {*str_op} is set to {op}
    and the first character is retained. In any case, the remaining
    characters are converted to an input symbol string with
    {rdag_string_iso_latin_1_cleanup} and
    {rdag_string_iso_latin_1_encode}. 
    
    If the string contains invalid characters or is too long, reports
    the error to {stderr}, sets {*nstr} to 0, and returns FALSE.
    Otherwise returns TRUE. */

static char *update_op_message[] =
  { 
    [update_op_ADD] = "adding words", 
    [update_op_SUB] = "removing words", 
    [update_op_MOD] = "adding or removing words", 
    [update_op_CHK] = "checking words"
  };

int main (int argc, char **argv)
  {
    prog_state_t *st = (prog_state_t *)notnull(malloc(sizeof(prog_state_t)), "no mem");

    /* Parse the command line: */
    st->cmd = parse_commands(argc, argv);
    st->icmd = 0;
    st->report_step = DEFAULT_REPORT_STEP;
    st->plot_step = 0;
    st->plot_fname = NULL;
    st->forgive = FALSE;
    

    /* Choose the default encoding. */
    /* st->enc = rdag_encoding_iso_latin_1(); */

    /* Create the initial (empty) automaton: */
    /* The input symbols are the encoded string characters. */
    /* The output symbols are 0 = reject, 1 = accept; only the last symbol matters. */
    int nn = 26; /* Number of bits in node number. */
    int ni = 8;  /* Number of bits in input symbols. */
    st->A = raut_new(nn, ni, 0);

    command_start_message(txtcat(PROG_NAME " " PROG_VERS " ", today()));

    do_show_args(stderr, &(st->cmd));
    
    char *old_doc = NULL;

    /* Command loop: */
    while (st->icmd < st->cmd.ne)
      {
        command_t *ci = &(st->cmd.e[st->icmd]);
        if (strcmp(ci->op, "-load") == 0)
          { do_load(st, ci->str1);
            old_doc = raut_doc_get(st->A);
          }
        else if (strcmp(ci->op, "-add") == 0)
          { do_process_strings(st, ci->str1, NULL, update_op_ADD); }
        else if (strcmp(ci->op, "-remove") == 0)
          { do_process_strings(st, ci->str1, NULL, update_op_SUB); }
        else if (strcmp(ci->op, "-modify") == 0)
          { do_process_strings(st, ci->str1, NULL, update_op_MOD); }
        else if (strcmp(ci->op, "-check") == 0)
          { do_process_strings(st, ci->str1, ci->str2, update_op_CHK); }
        else if (strcmp(ci->op, "-dump") == 0)
          { do_dump(st, ci->str1); }
        else if (strcmp(ci->op, "-list") == 0)
          { do_list(st, ci->str1); }
        else if (strcmp(ci->op, "-report") == 0)
          { st->report_step = ci->num; }
        else if (strcmp(ci->op, "-crunch") == 0)
          { st->crunch_step = ci->num; }
        else if (strcmp(ci->op, "-plot") == 0)
          { st->plot_fname = ci->str1;
            st->plot_step = ci->num;
          }
        else
          { fprintf(stderr, "invalid or unimplemented command \"%s\"\n", ci->op);
            assert(FALSE);
          }
        st->icmd++;
      }

    /* Create the automaton: */
    command_finish_message(PROG_NAME " " PROG_VERS);
    
    return 0;
  }

void do_show_args(FILE *wr, command_vec_t *cmd)
  {
    int i;
    for (i = 0; i < cmd->ne; i++)
      { command_t *ci = &(cmd->e[i]);
        fprintf(wr, "  %s", ci->op);
        if (ci->str1 != NULL) { fprintf(wr, " '%s'", ci->str1); }
        if (ci->str2 != NULL) { fprintf(wr, " '%s'", ci->str2); }
        if (ci->num != int64_NONE) { fprintf(wr, " %lld", ci->num); }
        if (i != cmd->ne-1) { fprintf(wr, " \\"); }
        fputc('\n', wr);
      }
  }

FILE *do_open_read(char *fname, char *what, FILE *def_file, char *def_fname)
  {
    if ((fname == NULL) || (strlen(fname) == 0))
      { fprintf(stderr, "reading %s from %s\n", what, def_fname); fflush(stderr);
        return def_file;
      }
   else
      { fprintf(stderr, "reading %s from %s\n", what, fname); fflush(stderr);
        return open_read(fname, FALSE);
      }
  }

FILE *do_open_write(char *fname, char *what, FILE *def_file, char *def_fname)
  {
    if ((fname == NULL) || (strlen(fname) == 0))
      { fprintf(stderr, "writing %s to %s\n", what, def_fname); fflush(stderr);
        return def_file;
      }
   else
      { fprintf(stderr, "writing %s to %s\n", what, fname); fflush(stderr);
        return open_write(fname, FALSE);
      }
  }

void print_doc(FILE *wr, char *cmt, char *prefix)
  { 
    if (cmt == NULL) { return; }
    while (*cmt != '\000')
      { /* Write a new line, advance {cmt} to start of next one or to '\000' */
        char c = (*cmt);
        fputs(prefix, wr);
        fputc(' ', wr);
        fputc(c, wr);
        while (c != '\n')
          { cmt++; c = (*cmt); 
            if (c == '\000') { c = '\n'; }
            fputc(c, wr);
          }
        if (*cmt == '\n') { cmt++; }
      }
  }

void do_load(prog_state_t *st, char *rd_fname)
  {
    command_start_message("loading automaton");
    FILE *rd = do_open_read(rd_fname, "automaton", stdin, "stdin");
    raut_free(st->A);
    st->A = raut_read(rd);
    if (rd != stdin) { fclose(rd); }
    print_doc(stderr, raut_doc_get(st->A), "  |");
    command_finish_message("loading automaton");
  }

void do_dump(prog_state_t *st, char *wr_fname)
  {
    command_start_message("dumping automaton");
    FILE *wr = do_open_write(wr_fname, "automaton", stdout, "stdout");
    raut_write(wr, st->A);
    if ((wr == stdout) || (wr == stderr)) { fflush(wr); } else { fclose(wr); }
    command_finish_message("dumping automaton");
  }

void do_process_strings
  ( prog_state_t *st,
    char *rd_fname,
    char *wr_fname,
    update_op_t op
  )
  {
    command_start_message(update_op_message[op]);

    /* Open the files: */
    FILE *rd = do_open_read(rd_fname, "strings", stdin, "stdin");
    FILE *flag_wr = do_open_write(wr_fname, "bad strings", stdout, "stdout");
    FILE *plot_wr = NULL;
    if ((st->plot_step != 0) && (op != update_op_CHK))
      { plot_wr = do_open_write(st->plot_fname, "plot data", stdout, "stdout"); }

    uint64_t line_count = 0;          /* Number of current line (from 1). */
    uint64_t bad_count = 0;           /* Number of illegal strings read. */
    uint64_t flagged_count = 0;       /* Number of flagged strings. */
    uint64_t input_string_count = 0;  /* Number of input strings processed. */
    uint64_t input_symbol_count = 0;  /* Number of input symbols in those strings. */
    
    /* Allocate a raw line buffer, with one extra byte just in case: */
    char *buf = (char *)notnull(malloc((MAX_LINE_BYTES+1)*sizeof(char)), "no mem");

    uint32_t nbuf; /* Significant bytes in {buf} (excluding final EOL). */

    /* Allocate a symbol buffer: */
    rdag_symbol_t *str = (rdag_symbol_t *)notnull(malloc(MAX_STRING_SYMBOLS*sizeof(rdag_symbol_t)), "no mem");
    uint32_t nstr;

    /* !!! Bring the {raut_update_strings} loop to here !!! */
    /* raut_update_strings(&read_raw_string, stderr, report, !forgive); */
    while (TRUE)
      { int rok = rdag_line_read(rd, &nbuf, buf, MAX_LINE_BYTES);

        if (rok == -1) { /* End of file: */ break; }

        /* We do have a line: */
        line_count++;
        if ((flag_wr != stderr) && (st->report_step != 0) && (line_count % st->report_step == 0))
          { fprintf(stderr, "%llu\n", line_count); }

        if (rok == -2)
          { /* Line too long: */
            data_error_message(rd_fname, line_count, "line too long", nbuf, buf);
            if (op != update_op_CHK) { demand(FALSE, "aborted"); }
            bad_count++;
          }
        else
          {
            /* Convert line to i-symbols {str[0..str_count-1]} and a specific op code {str_op}: */
            update_op_t str_op; /* Operation to apply to this particular string. */
            bool_t eok = get_string_from_line(rd_fname, line_count, nbuf, buf, op, &nstr, str, &str_op);
            if (! eok)
              { /* Something was wrong with the string: */
                if (op != update_op_CHK) { demand(FALSE, "aborted"); }
                bad_count++;
              }
            else
              { /* Crunch if requested (before plotting, so the plot starts crunched): */
                if ((st->crunch_step != 0) && ((input_string_count % st->crunch_step) == 0))
                  { raut_crunch(st->A); }
                  
                /* Plot if requested (before processing this string, so the plot starts at 0): */
                if ((plot_wr != NULL) && (st->plot_step != 0) && ((input_string_count % st->plot_step) == 0))
                  { write_plot_line(plot_wr, input_string_count, input_symbol_count, st->A); }
                  
                /* Process the string, set {flag} if exceptional: */
                raut_state_t u = raut_root_get(st->A); /* Current root state. */
                bool_t flag;
                switch(str_op)
                  {
                    case update_op_ADD:
                      { /* Resize automaton if empty: */
                        if (raut_node_max_alloc(st->A) <= 1) { raut_expand(st->A, FIRST_EXPAND); }
                        /* Add string, flag it if already there: */
                        raut_state_t v = raut_string_add(st->A, u, nstr, str);
                        flag = ((u.ac == v.ac) && (u.nd == v.nd));
                        if (! flag) { raut_root_set(st->A, v); }
                      }
                      break;
                    case update_op_SUB:
                      { /* Remove string, flag it if not there: */
                        raut_state_t v = raut_string_remove(st->A, u, nstr, str);
                        flag = ((u.ac == v.ac) && (u.nd == v.nd));
                        if (! flag) { raut_root_set(st->A, v); }
                      }
                      break;
                    case update_op_CHK:
                      { /* Check string, flag it if not there: */
                        flag = ! raut_state_accepts(st->A, u, nstr, str);
                      }
                      break;
                    case update_op_MOD:
                    default:
                      /* Should not happen: */
                      assert(FALSE);
                  }
                if (flag)
                  { /* Count the flagged string: */
                    flagged_count++;
                    /* Write it out, if appropriate: */
                    if ((! st->forgive) && (flag_wr != NULL)) { fprintf(flag_wr, "%.*s\n", nbuf, buf); }
                  }
                /* Count input symbols: */
                input_string_count++;
                input_symbol_count += nstr;
              }
          }
      }
    /* Final crunch in any case (before plotting, so the final point is crunched): */
    raut_crunch(st->A);

    /* Plot last data point if requested, close the file, and reset plotting: */
    if ((plot_wr != NULL) && (st->plot_step != 0))
      { write_plot_line(plot_wr, input_string_count, input_symbol_count, st->A);
        do_close_file(plot_wr);
        /* Turn off plotting until reset explicitly: */
        st->plot_step = 0;
        st->plot_fname = NULL;
      }
    
    /* Close the files if appropriate: */
    do_close_file(flag_wr);
    do_close_file(rd);

    fputc('\n', stderr);
    fprintf(stderr, "%12llu lines read\n", line_count);
    fprintf(stderr, "%12llu illegal strings\n", bad_count);
    fprintf(stderr, "%12llu input strings processed\n", input_string_count);
    fprintf(stderr, "%12llu input symbols processed\n", input_symbol_count);
    fprintf(stderr, "%12llu flagged strings\n", flagged_count);
    fprintf(stderr, "%12u reachable nodes\n", raut_node_max(st->A) + 1);

    command_finish_message(update_op_message[op]);
  }
 
void write_plot_line(FILE *wr, uint64_t nstrings, uint64_t nsymbols, raut_t *A)
  {
    uint64_t node_max_alloc = raut_node_max_alloc(A);
    uint64_t node_max = raut_node_max(A);
    raut_state_t root = raut_root_get(A);
    uint64_t reach_count = raut_reachable_node_count(A, root);
    
    fprintf(wr, " %12llu", nstrings);
    fprintf(wr, " %12llu", nsymbols);
    fprintf(wr, " %12llu", node_max_alloc);
    fprintf(wr, " %12llu", node_max);
    fprintf(wr, " %u %12u", root.ac, root.nd);
    fprintf(wr, " %12llu", reach_count);
    fprintf(wr, "\n");
  }

void do_close_file(FILE *f)
  { 
    if ((f == stdin) || (f == NULL))
      { return; }
    else if ((f == stdout) || (f == stderr))
      { fflush(f); }
    else
      { fclose(f); }
  }

bool_t get_string_from_line
  ( char * rd_fname, 
    uint64_t line_count, 
    uint32_t nbuf, 
    char buf[], 
    update_op_t op, 
    uint32_t *nstr, 
    rdag_symbol_t str[],
    update_op_t *str_op
  ) 
  {
    /* Get the raw string (minus leading code '+'/'-' if any) and string operation: */
    char *arg; /* The raw argument string. */
    uint32_t narg; /* Count of bytes in raw argument string. */
    switch(op)
      {
        case update_op_ADD:
        case update_op_SUB:
        case update_op_CHK:
          arg = buf;
          narg = nbuf;
          (*str_op) = op;
          break;
        case update_op_MOD:
          /* Require a leading '+' or '-', without leading whitespace: */
          if ((buf[0] != '+') && (buf[0] != '-'))
            { data_error_message(rd_fname, line_count, "invalid add/remove line", nbuf, buf);
              (*nstr) = 0;
              return FALSE;
            }
          (*str_op) = (buf[0] == '+' ? update_op_ADD : update_op_SUB);
          arg = buf + 1;
          narg = nbuf - 1;
          break;
        default:
          assert(FALSE);
      }

    /* Cleanup whitespace and NUL chars: */
    rdag_string_iso_latin_1_cleanup(&narg, arg);

    /* Convert the line to a symbol sequence {str[0..str_count-1]}: */
    bool_t eok = rdag_string_iso_latin_1_encode(narg, arg, nstr, str, MAX_STRING_SYMBOLS);
    if (! eok)
      { /* Something wrong with the string: */
        if ((*nstr) >= MAX_STRING_SYMBOLS)
          { data_error_message(rd_fname, line_count, "string too long", narg, arg); }
        else
          { data_error_message(rd_fname, line_count, "invalid char in string", narg, arg); }
        (*nstr) = 0;
        return FALSE;
      }
    return TRUE;
  }

void do_list(prog_state_t *st, char *wr_fname)
  {
    command_start_message("listing the accepted strings");

    FILE *wr = do_open_write(wr_fname, "accepted strings", stdout, "stdout");

    uint64_t string_count = 0;

    auto rdag_disp_t string_print(uint32_t nstr, rdag_symbol_t str[]);
      /* Prints the string {str[0..nstr-1]} to {wr} in the iso-latin-1 encoding. */

    rdag_disp_t string_print(uint32_t nstr, rdag_symbol_t str[])
      {
        rdag_string_iso_latin_1_write(wr, nstr, str);
        fputc('\n', wr);
        string_count++;
        if ((st->report_step != 0) && (wr != stderr) && (string_count % st->report_step == 0))
          { fprintf(stderr, "%llu\n", string_count); }
        return rdag_disp_FINE;
      }

    raut_enum_suffs(st->A, raut_root_get(st->A), &string_print);
    if ((wr == stdout) || (wr == stderr)) { fflush(wr); } else { fclose(wr); }
    fputc('\n', stderr);
    fprintf(stderr, "listed %llu strings\n", string_count);
    command_finish_message("listing the accepted strings");
  }

void do_debug(prog_state_t *st, char *wr_fname)
  {
    command_start_message("printing automaton");
    fprintf(stderr, "not implemented\n");
//     Util.PrintDoc(wr, A.doc, "  |");
//     raut_Print(wr, A, encoding);
//     fclose(wr);
    command_finish_message("printing automaton");
  }

command_vec_t parse_commands(int argc, char **argv)
  {
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);

    uint32_t ncmd = 0;
    command_vec_t cmd = command_vec_new(10); /* The commands are {cmd.e[0..ncmd-1]}. */

    while (argparser_next_is_keyword(pp))
      {
        /* Next command: */
        command_t ci = (command_t){ .op = NULL, .str1 = NULL, .str2 = NULL, .num = int64_NONE };
        ci.op = argparser_get_next(pp);

        if (strcmp(ci.op, "-load") == 0)
          { ci.str1 = parse_filename(pp); }
        else if (strcmp(ci.op, "-dump") == 0)
          { ci.str1 = parse_filename(pp); }
        else if (strcmp(ci.op, "-add") == 0)
          { ci.str1 = parse_filename(pp); }
        else if (strcmp(ci.op, "-remove") == 0)
          { ci.str1 = parse_filename(pp); }
        else if (strcmp(ci.op, "-modify") == 0)
          { ci.str1 = parse_filename(pp); }
        else if (strcmp(ci.op, "-check") == 0)
          { ci.str1 = parse_filename(pp);
            ci.str2 = parse_filename(pp);
          }
        else if (strcmp(ci.op, "-list") == 0)
          { ci.str1 = parse_filename(pp); }
        else if (strcmp(ci.op, "-debug") == 0)
          { ci.str1 = parse_filename(pp); }
        else if (strcmp(ci.op, "-report") == 0)
          { ci.num = argparser_get_next_uint(pp, 0, MAX_STRING_STEP); }
        else if (strcmp(ci.op, "-crunch") == 0)
          { ci.num = argparser_get_next_uint(pp, 0, MAX_STRING_STEP); }
        else if (strcmp(ci.op, "-plot") == 0)
          { ci.str1 = parse_filename(pp);
            ci.num = argparser_get_next_uint(pp, 0, MAX_STRING_STEP);
          }
        else if (strcmp(ci.op, "-expand") == 0)
          { ci.num = argparser_get_next_uint(pp, 0, MAX_EXPAND_ARG); }
        else
          { argparser_error(pp, "unrecognized or unimplemented option"); }
        command_vec_expand(&cmd, ncmd);
        cmd.e[ncmd] = ci;
        ncmd++;
      }
    command_vec_trim(&cmd, ncmd);
    return cmd;
  }
  
char *parse_filename(argparser_t *pp)
  {
    char *fname = argparser_get_next_non_keyword(pp);
    if (strcmp(fname,"-") == 0) { fname = NULL; }
    return fname;
  }

static clocks_t clocks;

void clocks_save(clocks_t *c)
  { c->rtime = real_time_usec(); 
    c->utime = user_cpu_time_usec();
    c->stime = system_cpu_time_usec(); 
  }

char *clocks_diff_string(clocks_t *c)
  { double r = real_time_usec() - c->rtime; 
    double u = user_cpu_time_usec() - c->utime;
    double s = system_cpu_time_usec() - c->stime;
    double M = 1.0e6; /* One second. */
    char *tim = NULL;
    asprintf(&tim, "(%.3f usr, %.3f sys, %.3f tot)", u/M, s/M, r/M);
    return tim;
  }

void command_start_message(char *msg)
  {
    fputc('\n', stderr);
    dashed_message("started", msg, NULL);
    clocks_save(&clocks);
  }

void command_finish_message(char *msg)
  {
    char *tim = clocks_diff_string(&clocks);
    dashed_message("finished", msg, tim);
    fputc('\n', stderr);
    free(tim);
  }

#define DASHED_MESSAGE_WIDTH 76

void dashed_message(char *msg1, char *msg2, char *msg3)
  {
    int len = 0; /* Number of characters printed. */
    bool_t some = FALSE; /* Some argument was printed. */

    /* Start dashes: */
    fprintf(stderr,"---"); len += 3;

    /* The three arguments, each with a leading blank: */
    if ((msg1 != NULL) && (strlen(msg1) !=  0))
      { fprintf(stderr," %s", msg1); len += 1 + strlen(msg1); some = TRUE; }
    if ((msg2 != NULL) && (strlen(msg2) !=  0))
      { fprintf(stderr," %s", msg2); len += 1 + strlen(msg2); some = TRUE; }
    if ((msg3 != NULL) && (strlen(msg3) !=  0))
      { fprintf(stderr," %s", msg3); len += 1 + strlen(msg3); some = TRUE; }

    /* A final blank if needed: */
    if (some) { fprintf(stderr," "); len += 1; }

    /* At least three more dahses: */
    fprintf(stderr,"---"); len += 3;

    /* Finish off the line: */
    int i;
    for (i = len; i < DASHED_MESSAGE_WIDTH; i++) { fputc('-', stderr); }
    fputc('\n', stderr);
  }

void data_error_message(char *fname, uint32_t nlin, char *msg, int nbuf, char buf[])
  {
    fprintf(stderr, "%s:%ul: ** error: %s\n", fname, nlin, msg);
    fprintf(stderr, "  line = «%.*s»\n", nbuf, buf);
    fflush(stderr);
  }

vec_typeimpl(command_vec_t,command_vec,command_t);

// void MakeDoc(old_doc, cmd, newDoc: char *): char *==
//   
//   {
//     with (wr == TextWr.New()){
//       Wr.PutText(wr, old_doc);
//       Wr.PutText(wr, "\n");
//       Wr.PutText(wr, "--- MaintainAutomaton ");
//       Wr.PutText(wr, PROG_VERSION);
//       Wr.PutText(wr, " --- ");
//       Wr.PutText(wr, Util.FmtDate(Util.GetDate()));
//       Wr.PutText(wr, " ---\n");
//       Wr.PutText(wr, "\n");
//       Wr.PutText(wr, cmd);
//       if ((NOT Text.Empty(newDoc))){
//         Wr.PutText(wr, "\n");
//         Wr.PutText(wr, newDoc);
//         Wr.PutText(wr, "\n");
//       ;};
//       return TextWr.ToText(wr)
//     ;};
//   ;}
//
// void SlurpText(FILE *rd): char *==
//   <* FATAL Rd.Failure, Wr.Failure, Rd.EndOfFile );
//   {
//     with (wr == TextWr.New()){
//       while (NOT Rd.EOF(rd)){
//         Wr.PutChar(wr, fgetc(rd))
//       ;};
//       return TextWr.ToText(wr)
//     ;};
//   ;}
//
// void DoWriteStats(prog_state_t *st) ==
//
//   void Print_counts()  ==
//     
//     {
//       fprintf(stderr,"_counts:\n");
//       fprintf(stderr,"-------\n");
//       with (counts == A._count(ARRAY OF raut_State{raut_root_get(A)})){
//         raut_Print_counts(st->log, counts)
//       ;};
//       fputc('\n', wr);
//     ;}
//
//   void PrintDegreeDistr()  ==
//     NAT *accum;
//     
//     {
//       with (
//         root == raut_root_get(A),
//         otrans ==  NEW(REF ARRAY OF NAT, NUMBER(Letter)),
//         itrans  == NEW(REF ARRAY OF NAT, root+1),
//         nStates == A.NStates(ARRAY OF raut_State{root}),
//         rStates == FLOAT(nStates)/100.0
//      ){
//         for (i = 0 TO root){ itrans^[i] = 0  ;};
//         for (i = 0 TO LAST(Letter)){ otrans^[i] = 0 ;};
//         for (s = raut_UnitState TO root){
//           if ((A.NPrefs(s)>0)){
//             INC(itrans^[A.InDeg(s)]);
//             INC(otrans[A.OutDeg(s)])
//           ;}
//         ;};
//         accum = 0;
//         fprintf(stderr,
//           "\nDistribution of indegrees:" &
//           "\n--------------------------\n\n" &
//           "In-Trans States   %     %acc\n\n"
//         );
//         for (i = 0 TO root){
//           with (it == itrans^[i]){
//             if ((it > 0)){
//               INC(accum, it);
//               fprintf(stderr,
//                 Fmt.Pad(Fmt.Int(i), 7) &
//                 Fmt.Pad(Fmt.Int(it), 7) &
//                 Fmt.Pad(Fmt.Real(FLOAT(it)/rStates, Fmt.Style.Fix, 2), 7) &
//                 Fmt.Pad(Fmt.Real(FLOAT(accum)/rStates, Fmt.Style.Fix, 2), 7)
//               );
//               fputc('\n', wr)
//             ;}
//           ;}
//         ;};
//         assert(accum==nStates );
//         accum = 0;
//         fprintf(stderr,
//           "\n\nDistribution of outdegrees:" &
//           "\n---------------------------\n\n" &
//           "Out-Trans States  %     %acc\n\n"
//         );
//         for (i = 0 TO LAST(Letter)){
//           with (ot==otrans^[i]){
//             if ((ot>0)){
//               INC(accum,ot);
//               fprintf(stderr,
//                 Fmt.Pad(Fmt.Int(i),7) &
//                 Fmt.Pad(Fmt.Int(ot),7) &
//                 Fmt.Pad(Fmt.Real(FLOAT(ot)/rStates, Fmt.Style.Fix, 2), 7) &
//                 Fmt.Pad(Fmt.Real(FLOAT(accum)/rStates, Fmt.Style.Fix, 2), 7)
//               );
//               fputc('\n', wr)
//             ;}
//           ;}
//         ;};
//         assert(accum == nStates );
//       ;};
//     ;}
//
//   void PrintLengthDistr()  ==
//     VAR
//       LenDis = ARRAY [0..MAX_STRING_SYMBOLS+1] OF NAT{0, ..};
//       uint32_t accum = 0;
//
//     void SuffixAction(READONLY sw: String) RAISES {} ==
//       {
//         with (n == NUMBER(sw)){
//           if ((n > MAX_STRING_SYMBOLS)){
//             INC(LenDis[MAX_STRING_SYMBOLS+1])
//           }else{
//             INC(LenDis[n])
//           ;}
//         ;}
//       ;}
//
//     <* FATAL Wr.Failure, Abort );
//     {
//       A.EnumSuffs(raut_root_get(A), SuffixAction);
//       fprintf(stderr,
//         "\nDistribution of string lengths:" &
//         "\n-----------------------------\n\n" &
//         "Lengths   Strings   %     %acc\n\n"
//       );
//       with (
//         nStrings == A.NSuffs(raut_root_get(A)),
//         rStrings == FLOAT(nStrings)/100.0,
//         hack == ARRAY BOOLEAN OF char *{" ", "+"}
//      ){
//         for (i = 0 TO  MAX_STRING_SYMBOLS+1){
//           with (ldi == LenDis[i]){
//             if ((ldi > 0)){
//               INC(accum, ldi);
//               fprintf(stderr,
//                 Fmt.Pad(Fmt.Int(i), 6) & hack[(i > MAX_STRING_SYMBOLS)] &
//                 Fmt.Pad(Fmt.Int(ldi), 7) &
//                 Fmt.Pad(Fmt.Real(FLOAT(ldi)/rStrings, Fmt.Style.Fix, 2), 7) &
//                 Fmt.Pad(Fmt.Real(FLOAT(accum)/rStrings, Fmt.Style.Fix, 2), 7)
//               );
//               fputc('\n', wr)
//             ;}
//           ;}
//         ;};
//         assert(accum == nStrings );
//       ;}
//     ;}
//
//   
//   {
//     command_start_message("printing automaton statistics");
//     Print_counts();
//     PrintDegreeDistr();
//     PrintLengthDistr();
//     fflush(st->log);
//     command_finish_message("printing automaton statistics");
//   ;}
//
