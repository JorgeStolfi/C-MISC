#include "stdio.h"
#include "string.h"
#include "ctype.h"
#define MAXLEN 512

/*
   Voynich Transcription Tool 
*/

/* Global variables */
/* These hold the options, defaults to do nothing */
int format=0;   /* Plain Ascii output format */
int high=0;     /* High Ascii left as is */
int kcomm=0;    /* Keep comments */
int hspace=0;   /* Do not touch dot spaces */
int uspace=0;   /* Do not touch comma spaces */
int ast=0;      /* Keep asterisk */
int brack=0;    /* Leave brackets */
int kfoli=0;    /* Leave foliation info */
int gaps=0;     /* Leave - sign */
int para=0;     /* Leave = sign */
int liga=0;     /* Keep ligature indication as is */
int white=0;    /* Leave white space */
int hold=0;     /* Leave placeholders % and ! */
int wrap=0;     /* Maintain line wrapping */
int wwidth;     /* New line wrapping limit */
int infarg= -1; /* Argument of input file name */
int oufarg= -1; /* Argument of output file name */
char auth=' ';  /* Name of transcriber */
char invar[27]; /* List of 'include page' options */
char exvar[27]; /* List of 'exclude page' options */ 
int npgopt = 0; /* Number of include/exclude options */
FILE *fin, *fout; /* File handles */

/* Some file stats */
int nlpart = 0   /* Number of (partial) lines read from stdin */;
int nlread = 0   /* Number of lines  from stdin after unwrapping */;
int nldrop = 0   /* Number of lines dropped */;
int nlhash = 0   /* Number of hash lines suppressed */;
int nlempt = 0   /* Number of empty lines suppressed */;
int nlwrit = 0   /* Number of lines for output */;
int nlwrap = 0   /* Number of wrapped lines written to output*/;

/* These give info about a complete line, set in GetLine/PrepLine */
int  hash=0         /* 1 if line starts with # */;
int  hastext=0      /* 1 if there is more than just comment */;
int  hasfoli=0      /* 1 if there are < > */;
int  newfoli=0      /* 1 if there are < > without a locus (i.e. newpage)*/;
char cwarn=' '      /* Warning character */;
char *folname= "      "  /* Folio name, e.g. f85r3 */;
char *locname= "     "   /* Locus name, e.g. C3, or blank */;
char lineauth= ' '       /* Name of transcriber or blank */;
char pgvar[27]      /* List of page variable settings */;

/* These track the text */
int incomm=0;   /* <>0 if inside { } */
int infoli=0;   /* <>0 if inside < > */
int ind_ligo= -1; /* Position of ligatures ( char */
int ind_ligc= -1; /* Position of ligatures ) char */
int ind_alto= -1; /* Position of alt. readings [ char */
int ind_altb= -1; /* Position of first alt. readings | char */
int ind_altc= -1; /* Position of alt. readings ] char */
int word0 = -1;   /* Position of first char of a word */
int word1 = -1;   /* Position of last char of a word */
int weir0 = -1;   /* Position of 'weirdo' dollar */
int spec0 = -1;   /* Position of 'special' ampsersand */
int spec1 = -1;   /* Position of 'special' semicolon */
int font = 1;     /* New character font for RTF file, 1= EVA-1 */
int cur_font = 5; /* Previous font for RTF file, illegal value */

/* Further global variables for certain options */
int concat = 0;   /* In GetLine: used for judging cr and spaces */
int highasc = 0;  /* In PrepLine: Ascii code of &..;  */
char cue = '.';   /* The 'space' after which wrapping is allowed */

/*-----------------------------------------------------------*/

void shiftl(b, index, nrlost)

/* Shift left by (nrlost) bytes the part of string (b)
   after (index) */

char *b; 
int index, nrlost;

{
  int ii;
  for (ii = index+1; ii<MAXLEN-nrlost; ii++) {
    b[ii] = b[ii+nrlost];
  }
  return;
}

/*-----------------------------------------------------------*/

int FindSpace(buf, width)

/* Locate rightmost hard space character in buf[0 .. width-1].
   Return index, or zero if buf is shorter, or <0 if none
   found */

char *buf; 
int width;

{
  int result = -1, ii = 0, voycode = 1;
  char cb;

  while ((cb = buf[ii]) && ii < width) {
    /* Dedicated check for comments.
       Sets voycode if outside comments */
    if (cb == '<' || cb == '{') voycode = 0;
    else if (cb == '>' || cb == '}') voycode = 1;
    else { 
      if (cb == cue || cue == (char) 0) {
        if (voycode) result = ii;
      }  
    }
    ii += 1;
  }
  if (cb == 0) result = 0;
  return result;
}

/*-----------------------------------------------------------*/
        
void showvar(code, vars, head)
char code, *vars;
int head;

/* Print all page variables */

{
  int i;
  if (head) { /* Print header legend */
    fprintf(stderr, "%c ", code);
    for (i=0; i<=26; i++) fprintf(stderr, "%c", i+64);
    fprintf(stderr, "\n");
  }
  fprintf(stderr, "%c ", code);
  for (i=0; i<=26; i++) fprintf(stderr, "%c", vars[i]);
  fprintf(stderr, "\n");
}

/*-----------------------------------------------------------*/

void clearvar()

/* Reset all page variables */

{
  int i;
  for (i=0; i<=26; i++) pgvar[i]=' ';
  /*
    fprintf(stderr,"clearing variables...\n");
    showvar('P',pgvar,1);
  */
}

/*-----------------------------------------------------------*/

int usepage(iopt)

/* Decide whether to select page based on comparing page
   options in pgvar[] with include and exclude options
   given by user. Page variable <at> abused to select 
   loci, making this a line selection option instead.
   The value of pgvar[0] is set in GetLine. If iopt=0 check
   only pgvar[0]. Otherwise check all others.
*/

int iopt;
{
  int i, select;
  
  select = 1;
  
  if (iopt) {
  
    for (i=1; i<=26; i++) {
      if (invar[i] != ' ') {
        if (invar[i] != pgvar[i]) select=0;
      }
      if (exvar[i] != ' ') {
        if (exvar[i] == pgvar[i]) select=0;
      }
    }
/*  showvar('P',pgvar,1); 
    showvar('I',invar,0);
    showvar('X',exvar,0);
*/ 
  } else {
/*  fprintf(stderr,"loc.type: %c\n",pgvar[0]); */
    if (invar[0] != ' ') {
      if (invar[0] != pgvar[0]) select=0;
    }
    if (exvar[0] != ' ') {
      if (exvar[0] == pgvar[0]) select=0;
    }
  
  }  

/* fprintf (stderr,"select(%d): %d\n",iopt,select); */

  return select;
}

/*-----------------------------------------------------------*/

void trackinit()

/* Initialise tracker */

{
    incomm = 0;    /* >0 if inside { } */
    infoli = 0;    /* >0 if inside < > */
    ind_ligo = -1; /* Position of ligatures ( char */
    ind_ligc = -1; /* Position of ligatures ) char */
    ind_alto = -1; /* Position of alt. readings [ char */
    ind_altb = -1; /* Position of alt. readings | char */
    ind_altc = -1; /* Position of alt. readings ] char */
    word0 = - 1;   /* Position of first char of a word */
    word1 = - 1;   /* Position of last char of a word */
    cur_font = 5;  /* Default is illegal to force first one */
}

/*-----------------------------------------------------------*/

int trackerr(buf, cget)

/* Print track error */

char *buf, cget;
{
  fprintf(stderr, "Offending character: %c\n", cget);
  fprintf(stderr, "Line parsed so far: %s\n", buf);
  return 1;
}

/*-----------------------------------------------------------*/

int Track(cb, index)
/* Keep track of comments, foliation, ligatures, etc */
/* Return 0 if OK, 1 if error */
/* Will only be called for line without # comment */

char cb;
int index;
{

  /* Evolve previous reading of bracket */
  if (incomm == -1) incomm = 1;
  if (incomm == -2) incomm = 0;
  if (infoli == -1) infoli = 1;
  if (infoli == -2) infoli = 0;
  
  /* If a word was just finished, start a new one */
  if (word1 >= 0) {
    word0 = -1; word1 = -1;
  }
  
  /* Check for in-line comments. Error if inside < >,
     warning if inside [ ] or ( ) */

  if (cb == '{') {
    if (infoli != 0 ) {
      fprintf(stderr,"%s\n","Illegal bracket inside <..>");
      return 1;
    }
    if (ind_ligo > ind_ligc || ind_alto > ind_altc) cwarn='{';
    if (incomm != 0) { /* Already open */
      fprintf(stderr,"%s\n","Illegal second open bracket");      
      return 1;
    }
    if (word0 >= 0 && word1 <0) word1 = index-1; /* End word here */ 
    incomm = -1; return 0;
  }
  if (cb == '}') {
    if (incomm < 1) { /* Not open */
      fprintf(stderr,"%s\n","Illegal close bracket");    
      return 1;
    }
    word0 = -1; word1 = -1; /* Prepare for new word */
    incomm = -2; return 0;
  }

  /* Continue checks only outside { } comments */
  if (incomm != 0) {
    if (incomm > 0) incomm += 1;
    return 0;
  }

  /* Check for foliation. It is already guaranteed that nothing
     else is open */
  if (cb == '<') {
    infoli = -1; return 0;
  }
  if (cb == '>') {
    if (infoli < 3) { /* Should really be longer */
      fprintf(stderr, "%s\n", "Foliation field too short");
      return 1;
    }
    infoli = -2; return 0;
  }

  /* Check for alternate reading brackets */
  if (cb == '[') {
    if (infoli != 0) { /* Not allowed inside < > */
      fprintf(stderr, "%s\n", "Illegal bracket inside <..>");
      return 1;
    }
    if (ind_alto > ind_altc) { /* Unclosed bracket already open */
      fprintf(stderr, "%s\n", "Illegal second open bracket");
      return 1;
    }
    ind_alto = index; ind_altb = -1; ind_altc = -1;
  } else if (cb == '|') {
    if (ind_alto < 0) {
      fprintf(stderr, "%s\n", "Illegal bar without [");
      return 1;
    }
    if (ind_altb < 0) ind_altb = index; /* Keep only left-most */
  } else if (cb == ']') {
    if (ind_alto < 0) { /* No bracket open */
      fprintf(stderr, "%s\n", "Illegal close bracket");
      return 1;
    }
    if (ind_altb < 0) { /* Vertical bar missing, issue warning */ 
      cwarn = '|';       
      ind_altb = ind_alto + 2; /* To continue reasonably */
    }
    ind_altc = index;
  }

  /* Check for ligature brackets, not allowed inside < > */

  if (cb == '(') {
    if (infoli != 0) { /* not allowed inside < > */
      fprintf(stderr, "%s\n", "Illegal bracket inside <..>");
      return 1;
    }
    if (ind_ligo > ind_ligc) { /* Unclosed paren already open */
      fprintf(stderr, "%s\n", "Illegal second open bracket");
      return 1;
    }
    ind_ligo = index; ind_ligc = -1;
    return 0;
  }
  if (cb == ')') {
    if (ind_ligo < 0) {
      fprintf(stderr,"%s\n","Illegal close bracket");
      return 1; /* No bracket open */
    }
    ind_ligc = index;
    return 0;
  }

  /* Inside foliation keep track of index */
  if (infoli > 0) infoli += 1;
  
  /* Further checks of word boundaries, not in foliation */
  if (infoli == 0) {
  
    if (cb == ',' || cb == '.' || cb == ' ' ||
        cb == '\n' || cb == '-') { /* Word delimiter */
      if (word0 >= 0) word1 = index - 1;
    } else { /* non-space char */
      if (word0 < 0) word0 = index;
    }

  }
  
  return 0;
}
/*-----------------------------------------------------------*/

void TrackLight(cb, index)
/* Same as track, but do not do error checking */
/* Will only be called for line without # comment */

char cb;
int index;
{

  font = 1;
  /* Evolve previous reading of bracket */
  if (incomm == -1) incomm = 1;
  if (incomm == -2) incomm = 0;
  if (infoli == -1) infoli = 1;
  if (infoli == -2) infoli = 0;
  
  /* If a word was just finished, start a new one */
  if (word1 >= 0) {
    word0 = -1; word1 = -1;
  }
  
  /* Check for in-line comments. */
  if (cb == '{') {
   if (word0 >= 0 && word1 < 0) word1 = index - 1; /* End word here */ 
    incomm = -1; font = 0; return;
  }
  if (cb == '}') {
    word0 = -1; word1 = -1; /* Prepare for new word */
    incomm = -2; font = 0; return;
  }

  /* Continue only outside { } comments */
  if (incomm != 0) {
    if (incomm > 0) incomm += 1;
    font = 0; return;
  }

  /* Check for foliation */
  if (cb == '<') {
    infoli = -1; font = 0; return;
  }
  if (cb == '>') {
    infoli = -2; font = 0; return;
  }

  /* Check for alternate reading brackets */
  if (cb == '[') {
    ind_alto = index; ind_altb = -1; ind_altc = -1;
  } else if (cb == '|') {
    if (ind_altb < 0) ind_altb = index; /* Keep only left-most */
  } else if (cb == ']') {
    ind_altc = index;
  }

  /* Check for ligature brackets */

  if (cb == '(') {
    ind_ligo = index; ind_ligc = -1;
    return;
  }
  if (cb == ')') {
    ind_ligc = index;
    return;
  }

  /* Inside foliation keep track of index */
  if (infoli > 0) {
    infoli += 1; font = 0;
  }
  
  /* Further checks of word boundaries, not in foliation */
  if (infoli == 0) {
  
    if (cb == ',' || cb == '.' || cb == ' ' ||
        cb == '\n' || cb == '-') { /* Word delimiter */
      if (word0 >= 0) word1 = index - 1;
    } else { /* non-space char */
      if (word0 < 0) word0 = index;
    }

  }
  
  return;
}

/*-----------------------------------------------------------*/

void OutChar(cb)

/* Write character cb to Ascii or RTF file using font 'font' */

char cb; 

{
  /* Main processing only required if format == RTF */
  if (format == 1) {
  
    /* Check font and enter font change if necessary */
    /* Font change is never necessary for a space char */
    if (cb != ' ') {
      if (font != cur_font) {
        fprintf (fout, "\\f%1d ", font);
        cur_font = font;
      }
    }
  
    /* Insert \ before some special characters */
    switch (cb) {
      case '\\':
      case '{':
      case '}': fputc ('\\', fout);
        break;
    }
  } 
   
  /* For all formats: print the character itself */
  fputc (cb, fout);
  
  /* If character just written was newline: */
  if (cb == '\n') {
    /* In RTF: add the correct command and force font definition */
    if (format == 1) {
      fprintf (fout, "\\par ");
      cur_font = 5;
    }
    /* Every newline implies an extra 'wrapped line' */
    nlwrap += 1;
  }
  
  return;
}

/*-----------------------------------------------------------*/

void OutString(buf)

/* Write string of characters to Ascii or RTF file */

char *buf; 

{
  int index = 0;
  char cb;
  
  /* RTF: hash comment always in comment font */
  if (hash) font = 0;
  while (cb = buf[index]) {
    /* RTF: if not hash comment: keep track of font */
    if (hash == 0) {
      TrackLight(cb, index);
    } 
    OutChar(buf[index++]);
    /* Complete lines are identified by a newline
       passed via routine OutString */ 
    if (cb == '\n') nlwrit += 1;
  }
  return;
}

/*-----------------------------------------------------------*/

int ParseOpts(argc, argv)
/* Parse command line options. There can be many and each should be
   of one of the following types:
   -pv with lower case p: one of:
     b,c,s,u,f,o,p,w,x where v can be 0-9, plus optionally -tC where 
     C is the transcriber code.
   -Pc or +Pc with upper case P: page variable where c is A-Z or
     0-9. If P=<at> then used for locus type, c=P for normal loci.
   <filename>:
     Maximum two, where first is input file name and second is
     output file name  */
int argc;
char *argv[];
{
  int i;
  char sign, optn, val;

  for (i=0; i<=26; i++) {
    invar[i]=' '; exvar[i]=' ';
  }

  for (i=1; i<argc; i++) {
    sign=argv[i][0];
    if (sign == '-' || sign == '+') {
    
      /* Normal option, not file name */
      
      optn=argv[i][1]; val=argv[i][2];

      /* process each one */

      if (optn >= '@' && optn <= 'Z') {

        /* Upper case: */
        if (sign == '+') {
          invar[optn-64]=val;
        } else {
          exvar[optn-64]=val;
        }
      } else {
      
        /* Lower case: */
        switch (optn) {
           case 'a':
             if (val == '0') high = 0; /* High ascii as is */
             else if (val == '1') high = 1;  /* High ascii to &123; */
             else if (val == '2') high = 2;  /* &123; to high ascii */
             break;
           case 'b':
             if (val == '0') white = 0; /* Keep blanks */
             else if (val == '1') {
               white = 1; hold = 1;     /* Strip them, and also % and ! */
             }
             break;
           case 'c':
             if (val == '0') kcomm=0; /* Keep comments */
             else if (val == '1') kcomm=1; /* Strip # */
             else if (val == '2') kcomm=2; /* Strip { } */
             else if (val == '3') kcomm=3; /* Strip all */
             break;
           case 'f':
             if (val == '0') kfoli=0; /* Keep foliation info */
             else if (val == '1') kfoli=1; /* Strip it */
             break;
           case 'o':
             if (val == '0') format=0; /* Ascii As Is */
             else if (val == '1') { /* RTF */
               format = 1; high = 2; liga = 4;
             }                  
             break;
           case 'p':
             if (val == '0') {
               gaps=0; para=0; /* Keep - and = */
             } else if (val == '1') {
               gaps=1; para=1; /* Gaps and = stripped */
             } else if (val == '2') {
               gaps=2; para=2; /* Gaps to extra space and = to extra NL */
             }
             break;
           case 'l':
             if (val == '0')  liga=0;     /* Leave ligature brackets */
             else if (val == '1') liga=1; /* Strip them */
             else if (val == '2') liga=2; /* Change to square brackets */
             else if (val == '3') liga=3; /* Change to capitalisation rule */
             else if (val == '4') liga=4; /* As 3, but change cth to CTh.. */
             break;
           case 'h':
             if (val == '0') uspace = 0;      /* Keep comma for uncertain space*/
             else if (val == '1') uspace = 1; /* Treat comma as dot */
             else if (val == '2') uspace = 2; /* Strip uncertain spaces */
             else if (val == '3') {
               uspace = 3;                    /* Turn uncertain spaces to * */
               ast = 1;                       /* And turn words with * into ? */
             }
             break;
           case 's':
             if (val == '0') {
               hspace = 0; cue = '.';      /* Keep dot for hard space*/
             } else if (val == '1') {
               hspace = 1; cue = ' ';      /* Turn hard space to blank */
             } else if (val == '2') {
               hspace = 2; cue = (char) 0; /* Strip hard space */
               uspace = 2;        /* And also uncertain spaces of course */
             }
             break;
           case 't':
             auth = val; hold = 1; /* At same time strip placeholders */
             break;
           case 'u':
             if (val == '0') {
               ast=0; brack=0; /* leave * and [] as is */
             } else if (val == '1') {
               ast=0; brack=1; /* Take first of [] */
             } else if (val == '2') {
               ast=0; brack=2; /* Turn [] to * */
             } else if (val == '3') {
               ast=1; brack=2; /* Turn word with * into ? */
             } else if (val == '4') {
               ast=2; brack=2; /* Remove ? or word with * */
             } else if (val == '5') {
               ast=3; brack=2; /* Remove line with ? * or [..] */}
             break;
           case 'w':
             if (val == '0') wrap=0;       /* Maintain line wrapping */
             else if (val == '1') wrap=1;  /* Unwrap all continuation lines */
             else {
               wrap=2;                     /* Re-wrap */
               wwidth = 20 * (val - '0');
             }
             break;
           case 'x':
             if (val == '0') {
               kcomm = 1; gaps = 1; para = 1;
               wrap = 1; white = 1; hold = 1;
             } else if (val == '1') {
               kcomm = 3; kfoli = 1; gaps = 1;
               para = 1; white = 1; hold = 1;
             } else if (val == '2') {
               kcomm = 3; kfoli = 1; gaps = 1;
               para = 1; white = 1; hold = 1; brack = 1;
             } else if (val == '3') {
               hspace = 1; uspace = 1; cue = ' ';
             } else if (val == '4') {
               hspace = 1; uspace = 2; cue = ' ';
             } else if (val == '5') {
               liga = 1; high = 2;
             } else if (val == '6') {
               liga = 4; high = 2;
             } else if (val == '7') {
               kcomm = 3; kfoli = 1; brack = 1; liga = 1; 
               white = 1; hold = 1; gaps = 1; para = 1;
               hspace = 1; cue = ' '; uspace = 1;
             } else if (val == '8') {
               kcomm = 3; kfoli = 1; brack = 1; liga = 1; 
               white = 1; hold = 1; gaps = 1; para = 1;
               hspace = 1; cue = ' '; uspace = 2;
             }
             break;
        }
      }
    } else {
      /* A file name. Check which of two */
      if (infarg < 0) {
        infarg = i;
      } else if (oufarg < 0) {
        oufarg = i;
      } else {
        fprintf (stderr, "%s\n", "Only two file names allowed");
        return 1;
      }
    }
  }
  return 0;
}

/*-----------------------------------------------------------*/

int DumpOpts(argc, argv)
/* Print information on selected options
   and open input and output files if required */
int argc;
char *argv[];
  {
  int i;

  fprintf (stderr,"Voynich Transcription Tool\n\n");

  /* Process each one */
  switch (format) {
     case 0: fprintf (stderr,"%s\n","Ascii output format");
         break;
     case 1: fprintf (stderr,"%s\n","Rich Text Format output");
         break;
  }       
  switch (high) {
     case 0: fprintf (stderr,"%s\n","Leave high ascii as is");
         break;
     case 1: fprintf (stderr,"%s\n","Expand high ascii to &..; notation");
         break;
     case 2: fprintf (stderr,"%s\n","Convert &..; notation to 1 byte");
         break;
  }       
  switch (kcomm) {
     case 0: fprintf (stderr,"%s\n","Keep all comments");
         break;
     case 1: fprintf (stderr,"%s\n","Remove hash comments");
         break;
     case 2: fprintf (stderr,"%s\n","Remove inline comments");
         break;
     case 3: fprintf (stderr,"%s\n","Remove all comments");
         break;
  }
  switch (kfoli) {
     case 0: fprintf (stderr,"%s\n","Keep foliation");
         break;
     case 1: fprintf (stderr,"%s\n","Remove foliation");
         break;
  }
  switch (liga) {
     case 0: fprintf (stderr,"%s\n","Keep ligature brackets");
         break;
     case 1: fprintf (stderr,"%s\n","Remove ligature brackets");
         break;
     case 2: fprintf (stderr,"%s\n","No ligatures, convert to square brackets");
         break;
     case 3: fprintf (stderr,"%s\n","Change ligature brackets to capitalisation");
          break;     
     case 4: fprintf (stderr,"%s\n","Use extended EVA capitalisation");
          break;

  }
  switch (uspace) {
     case 0: fprintf (stderr,"%s\n","Use comma for uncertain spaces");
         break;
     case 1: fprintf (stderr,"%s\n","Treat uncertain spaces as dot spaces");
         break;
     case 2: fprintf (stderr,"%s\n","Remove uncertain spaces");
         break;
     case 3: fprintf (stderr,"%s\n","Turn words next to uncertain spaces to ?");
         break;
  }
  switch (hspace) {
     case 0: fprintf (stderr,"%s\n","Use dot for normal spaces");
         break;
     case 1: fprintf (stderr,"%s\n","Use space for normal spaces");
         break;
     case 2: fprintf (stderr,"%s\n","Remove normal spaces");
         break;
  }
  switch (white) {
     case 0: fprintf (stderr,"%s\n","Keep white space");
         break;
     case 1: fprintf (stderr,"%s\n","Remove white space");
         break;
  }
  switch (hold) {
     case 0: fprintf (stderr,"%s\n","Keep % and ! interlinear placeholders");
         break;
     case 1: fprintf (stderr,"%s\n","Remove % and ! interlinear placeholders");
         break;
  }
  switch (brack) {
     case 0: fprintf (stderr,"%s\n","Keep alternate readings notation");
         break;
     case 1: fprintf (stderr,"%s\n","Take first of alternate readings");
         break;
     case 2: fprintf (stderr,"%s\n","Turn alternate readings into *");
         break;
  }
  switch (ast) {
     case 0: fprintf (stderr,"%s\n","Keep words containing *");
         break;
     case 1: fprintf (stderr,"%s\n","Turn words containing * into ?");
         break;
     case 2: fprintf (stderr,"%s\n","Remove ? and words containing *");
         break;
     case 3: fprintf (stderr,"%s\n","Remove lines containing ? or *");
         break;
  }
  switch (para) {
     case 0: fprintf (stderr,"%s\n","Keep paragraph end sign");
         break;
     case 1: fprintf (stderr,"%s\n","Remove paragraph end sign");
         break;
     case 2: fprintf (stderr,"%s\n","Replace para end sign by extra newline");
         break;
  }
  switch (gaps) {
     case 0: fprintf (stderr,"%s\n","Keep plant gap sign");
         break;
     case 1: fprintf (stderr,"%s\n","Remove plant gap sign");
         break;
     case 2: fprintf (stderr,"%s\n","Change plant gap sign into extra space");
         break;
  }
  switch (wrap) {
     case 0: fprintf (stderr,"%s\n","Maintain line wrapping");
         break;
     case 1: fprintf (stderr,"%s\n","Unwrap continuation lines");
         break;
     default: fprintf (stderr,"(Re)wrap lines at %3d \n", wwidth);
         break;
  }
  if (auth == ' ') {
    fprintf (stderr,"%s\n","Ignore transcriber ID");
  } else {
    fprintf (stderr,"%s %c\n","Use only data from transcriber",auth);
  }

  /* Locus splitter */
  if (invar[0] != ' ') {
    fprintf (stderr,"Include only locus type %c\n", invar[0]);
    npgopt += 1;
  }
  if (exvar[0] != ' ') {
    fprintf (stderr,"Exclude locus type %c\n", exvar[0]);
    npgopt += 1;
  }
  
  /* Page splitter */
  fprintf (stderr,"\nPage selection options (exclude overrules include):\n");

  for (i=1; i<=26; i++) {
    if (invar[i] != ' ') {
      fprintf (stderr, "- Include if variable %c set to %c\n", i+64, invar[i]);
      npgopt += 1;
    }
    if (exvar[i] != ' ') {
      fprintf (stderr, "- Exclude if variable %c set to %c\n", i+64, exvar[i]);
      npgopt += 1;
    }
  }
  if (npgopt == 0) fprintf (stderr, "- Include all.\n");

  /* Files */
  fprintf (stderr,"\nInput file: ");
  if (infarg >= 0) {
    fprintf(stderr, "%s\n",argv[infarg]);
    if ((fin = fopen(argv[infarg], "r")) == NULL) {
      fprintf(stderr, "File does not exist\n");
      return 1;
    }
  } else {
    fin = stdin;
    fprintf (stderr,"<stdin>\n");
  }
  
  fprintf (stderr,"Output file: ");
  if (oufarg >= 0) {
    fprintf(stderr, "%s\n", argv[oufarg]);
    if ((fout = fopen(argv[oufarg], "w")) == NULL) {
      fprintf(stderr, "Cannot open\n");
      return 1;
    }
  } else {
    fout = stdout;
    fprintf (stderr, "<stdout>\n");
  } 
  return 0;
}

/*-----------------------------------------------------------*/

int GetLine(buf)
/* Get line from stdin to buffer */
/* Return 0 if all OK, <0 if EOF, 1 if error */
/* Concatenate lines ending in backslash if wrap option >0 */

char *buf;
{
  int cr = 0, eod = 0, index = 0, iget, blank, ignore;
  char cget;

  /* Set these global parameters */
  hash = 0; hastext = 0; concat = 0;

  while (cr == 0 && eod == 0) {
    ignore = 0;
    iget = fgetc(fin);
    /* fprintf(stderr, "Index %3d,  char %3d\n", index, iget); */

    /* Check for end of file. Only allowed at first read */
    eod = (iget == EOF);
    if (eod) {
      if (index == 0) {
        return -1;       /* Correct EOF */
      } else {
        buf[index] = 0;  /* Add a null character for safety */
        fprintf (stderr, "EOF at record pos. %3d\n", index);        
        fprintf (stderr, "Line read so far: %s\n", buf); 
        return -2;
      }
    }

    cget = (char) iget;
    if (cget == ' ' || cget == '\t' || cget == '\n') {
      blank = 1;
    } else {
      blank = 0;
    }
     
    /* Check for comment line */

    if (cget == '#') {
      /* A hash symbol. Not allowed when unwrapping */
      if (concat > 0) {
        buf[index] = 0;  /* Add a null character for safety */
        fprintf (stderr, "Hash comment after continuation\n");        
        fprintf (stderr, "Line read so far: %s\n", buf); 
        return 1;
      }
      /* Check for valid hash comment line */
      if (index == 0) hash = 1;
    }
    
    /* Unwrap-processing only if this is not a hash comment */
    if (hash == 0 && wrap > 0) {
      if (concat == 0) { /* Search for backslash */
        if (cget == '\\') {
          concat = 1;
          ignore = 1;
        }
      } else if (concat == 1) { /* Search for newline */
        ignore = 1;
        if (cget == '\n') {
          concat = 2;
        } else if (blank == 0) {
          concat = 0;
          ignore = 0;
        }
      } else if (concat == 2) { /* Search for non-blank */
        if (blank == 0) {
          concat = 0;
          nlpart += 1; /* Only here count continuation */
        } else {
          ignore = 1;
        }        
      }
    }

    /* Process this char only if ignore not set */
    if (ignore == 0) {

      /* Add the character to the buffer */
      buf[index++] = cget;
      if (blank == 0) hastext = 1;
     
      /* Check for newline */
      if (cget == '\n') {
        cr = 1;
        buf[index] = 0; /* Add null at end of line */
      }
    } /* End if ignore == 0 */
  }
  return 0;
}

/*-----------------------------------------------------------*/

int PrepLine(buf1, buf2)
/* Preprocess line just read from stdin to buffer */
/* Return 0 if all OK, <0 if EOF, 1 if error */

char *buf1, *buf2;
{
  int ind1 = 0, ind2 = 0, varpt = 0;
  unsigned char cget, varid;

  int inauth = 0;   /* Track transcriber */
  int inlocu = 0;   /* Track locus ID */
  int ignore;       /* Processing white space */
  int ii;           /* Used to convert char(200) to &200; */

  /* Set these global parameters */
  hasfoli = 0; newfoli = 0;
  /* Hash and hastext were already set in GetLine */
  
  trackinit();

  while (cget = buf1[ind1]) {
    
    /* For hash comment just copy: */
    if (hash) {
      buf2[ind2++] = cget;
      /* fprintf(stderr, "Index %3d,  char %3d\n", ind2, cget); */
    } else {
      /* No hash comment: do all the processing */
      
      /* 1. Check for hash later in line and change to = */
      if (cget == '#') {
        cwarn = cget;
        cget = '=';
      }
  
      /* 2. Check for white space, leaving it untouched inside
         inline comments */
      ignore = 0;
      if (cget == ' ' || cget == '\t') {
        if (incomm == 0 && white == 1) ignore = 1;
      }

      /* 3. Process ligature brackets if desired */
      if (liga == 1) {
        if (cget == '(' || cget == ')' ) ignore = 1;
      }
      if (liga == 2) {
        if (cget == '(' ) cget = '[';
        if (cget == ')' ) cget = ']';
      }

      /* Process this char only if ignore not set */
      if (ignore == 0) {

        /* Check for open caret in 1st position */
        if (cget == '<' && incomm == 0) {
          if (ind1 != 0) {
            fprintf(stderr, "%s\n", "Foliation start not in first position");
            return 1;
          } else {
            hasfoli = 1;
          }
        }
  
        /* Track the text */
        if (Track(cget, ind1)) {
          buf2[ind2] = 0;  /* Add a null character for safety */
          return trackerr(buf2, cget);
        }

        /* Process Ascii(128-255) if desired */

        if (incomm == 0 && high == 1) {
          if (cget > 127) {
            buf2[ind2++]= '&';  
            ii = cget / 100; buf2[ind2++] = ii + 48;
            cget -= (100 * ii); ii = cget / 10;  buf2[ind2++] = ii + 48;
            ii = cget - (10 * ii); buf2[ind2++] = ii + 48;
            buf2[ind2++] = ';'; cget = '/'; ignore = 1;
          } 
        }
      
        /* Process &..; if desired */

        if (incomm == 0 && high == 2) {
          if (spec0 >= 0 && spec1 < 0) {
            if (cget >= '0' && cget <= '9') {
              /* Process numbers between & and ; */
              highasc = 10 * highasc + (cget - '0');
              /* fprintf(stderr, "Hiasc(I): %3d\n", highasc); */
            }
          }
          if (cget == '&' || cget == '$') {
            if (infoli != 0 ) {
              fprintf(stderr, "Illegal %c inside <..>\n", cget);
              return 1;
            }
            if (spec0 >= 0) { /* Already open */
              fprintf(stderr, "Illegal %c after ampersand\n", cget);      
              return 1;
            }
            spec1 = -1;
            if (cget == '&') {
              spec0 = ind2; highasc = 0; weir0 = -1;
            } else {
              weir0 = ind2;
            }
          }
          if (cget == ';') {
            if (infoli == 0) { /* Do not get confused by transcriber code */
              if (weir0 >= 0) { /* Just closing a weirdo */
                weir0 = -1;
              } else { /* Now it could be closing a special char. */
                if (spec0 <= 0) { /* Not open */
                  fprintf(stderr, "%s\n", "Illegal semi-colon");    
                  return 1;
                }
                spec1 = ind2; /* This will trigger code at end */
              }
            }
          }
        }

        /* Process locus/author inside foliation brackets */
        if (infoli != 0) {

          /* Evolve pointers */
          if (inlocu == -1) inlocu = 1;
          if (inlocu == -2) inlocu = 0;
          if (inauth == -1) inauth = 1;
          if (inauth == -2) inauth = 0;

          /* Check for < . ; and > */
          switch (cget) {
            case '<':
              folname = " "; locname = " "; 
              lineauth = ' '; newfoli = -1;
              break;
            case '>': /* Terminate (optionally) locus ID */
              if (inlocu > 0) {
                locname[inlocu] = 0; inlocu = -2;
              }
              if (newfoli < 0) {
                /* New folio. Prepare for var. reading */
                newfoli = 1; clearvar();
              }
              break;
            case '.': /* Start of locus */
              inlocu = -1; newfoli = 0;
              break;
            case ';': 
              if (inlocu <= 0) {
                fprintf(stderr, "%s\n", "Start of author ID without locus");
                return 1;
              }
              /* Terminate locus ID */
              locname[inlocu] = 0; inlocu = -2;
              /* Start for author ID */
              inauth = -1;
              break;
            default:
              /* Any other char. Check where it belongs */
              if (inauth > 1) {
                fprintf(stderr, "%s\n", "Illegal author ID length");
                return 1;
              }
              if (inauth == 1) { /* Get author ID */
                lineauth = cget;
              } else {
                if (inlocu > 0) { /* Add to locus */
                  locname[inlocu-1] = cget; 
                } else {
                  /* Add to folio and prepare terminator*/
                  folname[infoli-2] = cget; 
                  folname[infoli-1] = 0;
                }
              }
              break;
          }
          /* Advance pointers */
          if (inlocu > 0) inlocu += 1;
          if (inauth > 0) inauth += 1;
        } /* End if infoli != 0 */

        /* Check for page variables, only if start of new page found */
        if (newfoli == 1) {

          if (varpt == -3) {
            if (cget < 'A' || cget > 'Z') {
              cwarn = '$'; varpt = 0;
            } else {
              varid = cget;
            }
          }

          if (varpt == -1) pgvar[varid-64] = cget;

          /* Advance pointer */
          if (varpt < 0) varpt += 1;

          /* Only check for a new $ when not already in sequence */
          if (varpt == 0) {
            if (cget == '$') {
              varid = ' '; varpt = -3;
            }
          } 
        } /* End if newfoli == 1 */

        /* Add the character to the buffer */
        buf2[ind2++] = cget;
        /* Check for complete &..; code */
        if (high == 2) {
          if (spec1 >= 0) {
            /* fprintf(stderr, "Hiasc(F): %3d\n", highasc); */
            buf2[spec0] = (char) highasc;
            ind2 = spec0 + 1;
            spec0 = spec1 = -1;
          }
        }  
      } /* End if ignore == 0 */
      
    } /* End of switch hash / no hash */
    ind1 += 1;
 
  } /* End of input string */
  buf2[ind2] = 0; /* Complete the output buffer */

  /* Now the entire line was read, do some more checks */

  if (hash == 0) {
    if (incomm) {
      fprintf(stderr, "%s", "Unclosed comment\n");
      return 1;
    }
    if (infoli) {
      fprintf(stderr, "%s", "Unclosed foliation comment\n");
      return 1;
    }
    if (ind_alto > ind_altc) {
      fprintf(stderr, "%s", "Unclosed alternate reading\n");
      return 1;
    }
    if (ind_ligo > ind_ligc) {
      fprintf(stderr, "%s", "Unclosed ligature\n");
      return 1;
    }
    
    /* Set pgvar[0] if a new locus was read */
    if (hasfoli && (newfoli == 0)) {
      pgvar[0] = locname[0];
      if (pgvar[0] >= '0' && pgvar[0] <= '9') pgvar[0] = 'P';
    }
  }
  return 0;
}

/*-----------------------------------------------------------*/

int ProcRead(buf)
/* Process uncertain readings and/or ligature capitalisation rule.
   Read and write to same buffer
   Return 0 if all OK, -1 if line to be deleted,
   1 if error */
   
char *buf;
{
  int index, ii, ret=0, locast, locq;

  char cb;

  /* Check if anything needs to be done at all */
  if (brack == 0 && ast == 0 && liga < 3) return ret;

  /* If it is empty, nothing either */
  if (hastext == 0) return ret;

  /* If it is a hash line, nothing either */
  if (hash != 0) return ret;

  /* First loop over line takes care of ligature
     matters and the [] brackets */
  /* Initialise a few things */
  trackinit(); index = 0;

  while (cb = buf[index]) {

    /* Track the text */
    TrackLight(cb, index);
    
    /* If extended capitalisation is used and the character is h */
    if (liga == 4 && cb == 'h' && index > 0) {
      /* Only if not inside (..) and not inside {..} */
      if (incomm == 0 && ind_ligo <0) {
        switch (buf[index-1]) {
          case 's':
          case 't':
          case 'k':
          case 'p': 
          case 'f':
            buf[index-1] = toupper(buf[index-1]);
        }
      }
    }
     
    /* If a complete set of liga parens found: */
    if (ind_ligc >= 0) {
      if (liga >= 3) {

        /* Capitalise all but last char */
        index = ind_ligc - 2;
        for (ii=ind_ligo; ii<index; ii++) {
          buf[ii] = toupper(buf[ii+1]);
        }
        buf[index] = buf[index+1];
        /* If it ends with a quote, un-capitalise the one before */
        if (buf[index] == '\'') buf[index-1] = tolower(buf[index-1]);

        /* Shift left remainder of line */
        shiftl(buf, index, 2);
      }

      /* reset pointers */
      ind_ligo = -1; ind_ligc = -1;
    }

    /* If a complete set of alt.read brackets found: */
    if (ind_altc >= 0) {
      if (brack == 1) {

        /* Process into best guess: */
        index = ind_altb - 2;
        for (ii=ind_alto; ii<=index; ii++) {
          buf[ii] = buf[ii+1];
        }

      } else if (brack == 2) {

        /* Process into asterisk */
        index = ind_alto;
        buf[index] = '*'; 
      }

      /* Shift left remainder of line */
      shiftl(buf, index, ind_altc-index);
      
      /* reset pointers */
      ind_alto= -1; ind_altb= -1; ind_altc= -1;
    }
    index += 1;
  }
  
  /* If necessary, second loop over line takes care
     of uncertain word readings */
  locast = -1 ; locq= -1; 
  if (ast == 0) return ret;
  
  /* Initialise tracker */
  trackinit(); index = 0;
  /* fprintf(stderr,"%s\n",buf); */
  while (cb = buf[index]) {

    /* Track the text */
    TrackLight(cb, index);

    /* Look for asterisk and ? */
    if (cb == '*') locast = index;
    if (cb == '?') locq = index; 

    if (word1 >= 0) {
      /* A complete word was read: check * and ? */
      /* fprintf(stderr,"word from %d to %d\n",word0,word1); */
      if (locast >= 0 || locq >= 0) {
        /* fprintf(stderr,"ast/q:  %d / %d\n",locast,locq); */
        if (ast == 1) {
          /* Action: change word into ? */
          index = word0;
          buf[index] = '?';
          /* Shift left remainder of line */
          shiftl(buf, index, word1-word0);
          index += 1; word0 = -1; word1 = -1;
        
        } else  if (ast == 2) {
          
          /* Action: drop this word */
          index = word0 - 1;
          shiftl(buf, index, word1 - word0 + 1);
          index += 1; word0 = -1; word1 = -1;
        
        } else  if (ast == 3) {
          
          /* Action: drop this line. */  
          return -1;    
        }  
        /* fprintf(stderr,"%s\n\n",buf); */
        /* Reset pointers */
        locast = -1; locq = -1;  
      }        

    }
    index += 1;
  }
  return 0;
}

/*-----------------------------------------------------------*/

int ProcSpaces(buf1, buf2)
/* Process spaces in buffer. This treats occurrences of
   comma, dot, hyphen, equal, percent and exclam in input text.
   Return 0 if all OK, 1 if error */
char *buf1, *buf2;
{
  int indin=0, indout=0, eol=0;
  int addchar;
  char cb, cbo;

  /* Some standard 'track' initialisations per line */
  trackinit();

  /* If no text, do little */
  if (hastext == 0) {
    buf2[0] = 0;
    return 0;
  }

  while (eol == 0) {
    cb = buf1[indin];
    eol = (cb == (char) 0);

    if (hash == 0) {
      /* Track the text */
      TrackLight(cb, indin);
    }

/*  if (comm == -1) fprintf(stderr,"%s","Start of comment\n");
    if (comm == -2) fprintf(stderr,"%s","End of comment\n");
    if (folio == -1) fprintf(stderr,"%s","Start of foliation\n");
    if (folio == -2) fprintf(stderr,"%s","End of foliation\n");  */

    indin += 1;
    addchar = 1;

    /* Check outside comments only */
    if (hash == 0 && incomm == 0 && infoli == 0) {

      /* Check plant gap hyphen */
      if (cb == '-') {
        if (gaps == 1) addchar = 0;  /* Just remove it */
        else if (gaps == 2) {
          cb = '.';                  /* Convert to extra dot space */
        }
      }

      /* Check real and half spaces */
      cbo = cb;
      if (cb == ',') {
        if (uspace == 1) cb = cbo = '.';  /* Treat like dot space*/
        else if (uspace == 2) addchar=0;  /* Convert to blank */
        else if (uspace == 3) cbo = '*';  /* Convert to unreadable */
      }
      if (cb == '.') {
        if (hspace == 1) cbo = ' ';       /* Just use blank instead */
        else if (hspace == 2) addchar=0;  /* Convert to nothing */
      }
      cb = cbo;
      
      /* Check equal sign (para end) */
      if (cb == '=') {
        if (para == 1) addchar = 0;     /* Convert to nothing */
        else if (para == 2) cb = '\n';  /* Extra CR instead of = */
      }

      /* Here for percent and exclamation mark */
      if (hold == 1) {
        if (cb == '!' || cb == '%')  addchar = 0;   /* Convert to nothing */
      }
    }

    while (addchar--) {
      buf2[indout++] = cb;
    }
  }
  return 0;
}

/*-----------------------------------------------------------*/

int PutLine(buf)
/* Write buffer to output, optionally skipping comments,
   foliation info and a few other things */

/* Return 0 if all OK, 1 if error */

char *buf;
{
  int index=0, indout=0, eol=0, output, ii;
  char cb;
  char wrapbuf[MAXLEN];

  if (hash != 0) {
    /* Block comment: print immediately or exit  */
    if (kcomm == 1 || kcomm == 3) {
      nlhash +=1;
    } else {
      OutString(buf);
    }
    return 0;
  }

  /* If the line is empty, only print if whitespace kept */
  if (white == 1 && hastext == 0) {
    nlempt +=1 ; return 0;
  }
     
  /* For 'wrong author': same thing */
  if (auth != ' ') {
    if (lineauth != ' ' && lineauth != auth) {
      nldrop += 1; return 0;
    }
  }

  /* The usual initialisation for 'track' */
  trackinit();

  while (cb = buf[index]) {
    output = 1;
    TrackLight(cb, index);

    /* Check inline comments */
    if (kcomm > 1) {
      if (incomm != 0) output = 0;
    }

    /* Check foliation comments */
    if (kfoli == 1) {
      if (infoli != 0) output = 0;
    }

    /* Now output if required */
    if (output != 0) {
      wrapbuf[indout++] = cb;
    }
    index += 1;
  }
  
  /* All characters processed and added to wrapbuf */
  if (indout == 0) {
    nlempt += 1;
    return 0;
  } else if (wrapbuf[0] == '\n' && white == 1) {
    nlempt += 1;
    return 0;
  } else {
    wrapbuf[indout] = 0;
  }

  /* Now write wrapbuf whole or in pieces */
  /* Another initialisation for 'track' needed first */
  trackinit();
  
  if (wrap < 2) {
    OutString(wrapbuf);
  } else {
    /* First write all the parts */
    while ((index = FindSpace(wrapbuf, wwidth)) > 0) {
      for (ii=0; ii<=index; ii++) {
        cb = wrapbuf[ii];
        TrackLight(cb, index);
        OutChar(cb);
      }
      OutChar('\\');
      OutChar('\n');
      wrapbuf[0] = ' ';
      shiftl(wrapbuf, 0, index);
    }
    /* Write the last bit or handle the 'no dot found' case */
    if (index == 0) {
      OutString(wrapbuf);
    }
    if (index < 0) {
      fprintf(stderr, "%s", "No space found for line wrapping\n");
      return 1;
    }
  }

  return 0;
}

/*-----------------------------------------------------------*/

void main(argc, argv)
int argc;
char *argv[];
{
  char orig[MAXLEN], buf1[MAXLEN], buf2[MAXLEN];
  int ilin = 0, iprep = 0, iproc = 0, iout = 0;
  int selpage, selloc;

  /* Why not indeed: */
  char *what = "@(#)vtt\t1.6\t1997/09/30\tRZ\n";

  /* Parse command line options */
  if (ParseOpts(argc, argv)) {
    fprintf (stderr, "%s\n", "Error parsing command line");
    return;
  }
 
  /* Dump options and open any files */  
  if (DumpOpts(argc, argv)) {
    fprintf (stderr, "%s\n", "Error opening file(s)");
    return;
  }  
  fprintf (stderr, "\n%s\n", "Starting...");

  clearvar();

  /* Generate output file header (RTF) */
  if (format == 1) {
    fprintf (fout, "%s", "{\\rtf1\\ansi {\\fonttbl");
    fprintf (fout, "%s", "{\\f0\\fmodern\\fcharset0\\fprq1 Courier New;}");
    fprintf (fout, "%s", "{\\f1\\fnil\\fcharset0\\fprq2 EVA Hand 1;}}\n");
    fprintf (fout, "%s", "\\f1\\fs20 ");
  }
  
  /* Initial value of selpage (i.e. prior to any 'new folio')
     depends on whether page selection options were specified */
  if (npgopt == 0) selpage = 1;
  else selpage = 0;
  
  /* Main loop through input file */

  while (ilin == 0) {

    /* Read one line to buffer. 
       This concatenates lines if required but that's all */

    ilin = GetLine(orig);
    if (ilin < 0) {  /* Normal EOF */

      /* Generate output file trailer (RTF) */
      if (format == 1) {
        fprintf (fout, "}\n");
      }

      /* Print statistics */
      if (wrap == 0) {
        fprintf (stderr, "\n%7d lines read in\n", nlread);
      } else {
        fprintf (stderr, "\n%7d lines read in\n", nlpart);
        fprintf (stderr, "%7d lines after unwrapping\n", nlread);
      }
      fprintf (stderr, "%7d lines de-selected\n", nldrop);
      fprintf (stderr, "%7d hash comment lines suppressed\n", nlhash);
      fprintf (stderr, "%7d empty lines suppressed\n", nlempt);
      fprintf (stderr, "%7d lines written to output\n", nlwrit);
      if (wrap > 1) {
        fprintf (stderr, "%7d lines after wrapping\n", nlwrap);
      }
      return;
    }
    else if (ilin > 0) {
      fprintf (stderr, "%s\n", "Error reading line from stdin");
      return;
    }
    nlread += 1; nlpart += 1;
    /* fprintf (stderr, "%2d %s\n", nlread, orig); */
    
    /* Preprocess line */
    /* This also keeps track of foliation and comments, and warns
       about unclosed brackets */
    cwarn = ' ';
    iprep = PrepLine(orig, buf1);
    if (iprep != 0) {
      fprintf (stderr, "%s\n", "Error preprocessing line");
      fprintf (stderr, "Line: %s\n", orig);
      return;
    }
    if (cwarn != ' ') {
      fprintf (stderr, "%c %s\n", cwarn, "warning for line:");
      fprintf (stderr, "%s\n", orig);
      cwarn = ' ';
    }

    /* Check page variables only when new page read */
    if (newfoli == 1) {
      selpage = usepage(1);
      /*
        fprintf (stderr,"Variables on this new folio:\n");
        (void) showvar();
        fprintf (stderr,"Select page: %d\n",selpage);
      */
    }

    /* Check locus type only if page in principle OK */
    if (selloc = selpage) selloc = usepage(0);
  
    /* Only continue if page/locus as requested */
    if (selloc == 0) {
      nldrop += 1;
    } else {

      /* Process spaces */
      if (ProcSpaces(buf1, buf2) != 0) {
        fprintf (stderr, "%s", "Error processing spaces\n");
        fprintf (stderr, "Line: %s\n", orig);
        return;
      }

      /* Process uncertain readings */

      iproc = ProcRead(buf2);
      if (iproc > 0) {
        fprintf (stderr, "%s\n", "Error processing uncertain readings");
        fprintf (stderr, "Line: %s\n", orig);
        return;
      } else {

        /* Further processing and output are skipped if
           locus was not selected */
        if (iproc == 0) {
      
          /* Write line to stdout. Here the optional comment and/or
             foliation removal and/or author selection is handled, 
             as well as any re-wrapping. */

          if (PutLine(buf2) != 0) {
            fprintf (stderr, "Line: %s\n", orig);
            return;
          }
        } else {
          nldrop += 1;
        }
      }
    }  
  }
}
