#include "stdio.h"
#include "string.h"
/* Some global variables */
/* These hold the options, defaults to do nothing */
int kcomm=0     /* Keep comments */;
int hspace=0    /* Do not touch dot spaces */;
int uspace=0    /* Do not touch comma spaces */;
int ast=0       /* Keep asterisk */;
int brack=0     /* Leave brackets */;
int kfoli=0     /* Leave foliation info */;
int gaps=0      /* Leave - sign */;
int para=0      /* Leave = sign */;
int liga=0      /* Keep ligature indication as is */;
int white=0     /* Leave white space */;
int hold=0      /* Leave placeholders % and ! */;
char auth=' '   /* Name of transcriber */;

/* These give info about a complete line, set in getlin */
int  hash=0         /* 1 if line starts with # */;
int  hastext=0      /* 1 if there is more than just comment */;
int  hasfoli=0      /* 1 if there are < > */;
char cwarn=' '      /* Warning character */;
char *folname= "      "  /* Folio name, e.g. f85r3 */;
char *locname= "     "   /* Locus name, e.g. C3, or blank */;
char lineauth= ' '       /* Name of transcriber or blank */;

/* These track the text */
int incomm=0   /* <>0 if inside { } */;
int infoli=0   /* <>0 if inside < > */;
int ind_ligo= -1 /* Position of ligatures ( char */;
int ind_ligc= -1 /* Position of ligatures ) char */;
int ind_alto= -1 /* Position of alt. readings [ char */;
int ind_altb= -1 /* Position of first alt. readings | char */;
int ind_altc= -1 /* Position of alt. readings ] char */;

/*-----------------------------------------------------------*/

int main(argc,argv)
int argc;
char *argv[];
{
  char orig[512], buf1[512], buf2[512];
  char *dmy;
  int ipa, ilin=0, iproc=0, iout=0;

  /* Just for the heck of it: */
  char *what = "@(#)vtt\t1.3\t1996/09/25\n";

  /* Parse command line options */

  ipa = parseopts(argc,argv);
  ipa = dumpopts();
  fprintf (stderr,"\n%s\n","Starting...");

  /* Main loop through input file */

  while (ilin == 0) {

    /* Read one line to buffer. This already keeps track
       of foliation and comments, and warns about unclosed
       brackets */

    cwarn = ' ';
    ilin=getlin(orig);
    if (ilin < 0) return 0; /* Normal EOF */
    else if (ilin > 0) {
      fprintf (stderr,"%s\n","Error reading line:");
      fprintf (stderr,"%s\n",orig);
      return 1;
    }
    if (cwarn != ' ') {
      fprintf (stderr,"%c %s\n",cwarn,"warning for line:");
      fprintf (stderr,"%s\n",orig);
      cwarn=' ';
    }

    /* Process uncertain readings */

    dmy=strcpy(buf1,orig);
    iproc=procread(buf1);
    if (iproc != 0) {
      fprintf (stderr,"%s\n","Error processing uncertain reads");
      fprintf (stderr,"Line: %s\n",orig);
      return 1;
    }

    /* Process spaces */

    iproc=procspc(buf1,buf2);
    if (iproc != 0) {
      fprintf (stderr,"%s","Error processing spaces\n");
      fprintf (stderr,"Line: %s\n",orig);
      return 1;
    }

    /* Write line to stdout. Here the optional comment and/or
       foliation removal and/or author selection is handled */

    iout=putlin(buf2);
  }
}

/*-----------------------------------------------------------*/

int parseopts(argc,argv)
/* Parse command line options. There can be many and each should be
   of the type -pv where p = c,s,u,f,p,w and v can be 0-3 ,
   plus optionally -tC where C is the transcriber code */
int argc;
char *argv[];
{
  int i;
  char optn,val;
  for (i=2; i<=argc; i++) {
    optn=argv[i-1][1]; val=argv[i-1][2];

    /* Now process each one */
    switch (optn) {
      case 'c':
         if (val == '0') kcomm=0; /* Keep comments */
         else if (val == '1') kcomm=1; /* Strip them */
         break;
      case 'f':
         if (val == '0') kfoli=0; /* Keep foliation info */
         else if (val == '1') kfoli=1; /* Strip it */
         break;
      case 'p':
         if (val == '0') {
           gaps=0; para=0; /* Keep - and = as is */
         } else if (val == '1') {
           gaps=1; para=1; /* Gaps to two spaces and = stripepd */
         } else if (val == '2') {
           gaps=1; para=2; /* Gaps to two spaces and = to extra NL */
         }
         break;
      case 'l':
         if (val == '0')  liga=0; /* Leave ligature brackets */
         else if (val == '1') liga=1; /* Strip them */
         else if (val == '2') liga=2; /* Change to square brackets */
         else if (val == '3') liga=3; /* Change to capitalisation rule */
         break;
      case 's':
         if (val == '0') {
           hspace=0; uspace=0; /* Keep dot and comma */
         } else if (val == '1') {
           hspace=1; uspace=1; /* Turn both to space */
         } else if (val == '2') {
           hspace=1; uspace=2; /* Keep hard space but strip uncertain */
         } else if (val == '3') {
           hspace=2; uspace=2; /* Strip all spaces */
         }
         break;
      case 't':
         auth=val; hold=1; /* At same time strip placeholders */
         break;
      case 'u':
         if (val == '0') {
           ast=0; brack=0; /* leave * and [] as is */
         } else if (val == '1') {
           ast=0; brack=1; /* Take first of [] */
         } else if (val == '2') {
           ast=0; brack=2; /* Turn [] to * */
         } else if (val == '3') {
           /* ast=3; brack=2; Not yet supported */
         }
         break;
      case 'w':
         if (val == '0') white=0; /* Keep whitespace */
         else if (val == '1') {
           white=1; hold=1; /* Strip it, including % and ! */
         }
         break;
      case 'x':
         kcomm=1; kfoli=1; white=1; hold=1;
         if (val == '1') {
           hspace=1; uspace=1;
         } else if (val == '2') {
           hspace=1; uspace=2;
         } else if (val == '3') {
           hspace=2; uspace=2;
         } else if (val == '4') {
           ast=0; brack=1;
         } else if (val == '5') {
           ast=0; brack=2;
         }
         break;
    }
  }
  return 0;
}

/*-----------------------------------------------------------*/

int dumpopts()
/* Print information on selected options */
{
  fprintf (stderr,"Voynich Transcription Tool\n\n");

  /* Process each one */
  switch (kcomm) {
     case 0: fprintf (stderr,"%s\n","Leaving comments in");
         break;
     case 1: fprintf (stderr,"%s\n","Stripping comments");
         break;
  }
  switch (kfoli) {
     case 0: fprintf (stderr,"%s\n","Leaving foliation in");
         break;
     case 1: fprintf (stderr,"%s\n","Stripping foliation");
         break;
  }
  switch (liga) {
     case 0: fprintf (stderr,"%s\n","Keep ligature brackets");
         break;
     case 1: fprintf (stderr,"%s\n","Remove ligature brackets");
         break;
     case 2: fprintf (stderr,"%s\n","No ligatures, convert to square brackets");
         break;
     case 3: fprintf (stderr,"%s\n","Ligature brackets changed to capitalisation");
         break;
  }
  switch (uspace) {
     case 0: fprintf (stderr,"%s\n","Use comma for uncertain spaces");
         break;
     case 1: fprintf (stderr,"%s\n","Use space for uncertain spaces");
         break;
     case 2: fprintf (stderr,"%s\n","Strip uncertain spaces");
         break;
  }
  switch (hspace) {
     case 0: fprintf (stderr,"%s\n","Use dot for normal spaces");
         break;
     case 1: fprintf (stderr,"%s\n","Use space for normal spaces");
         break;
     case 2: fprintf (stderr,"%s\n","Strip normal spaces");
         break;
  }
  switch (white) {
     case 0: fprintf (stderr,"%s\n","Keep white space as is");
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
     case 0: fprintf (stderr,"%s\n","Leave alternate readings as is");
         break;
     case 1: fprintf (stderr,"%s\n","Take first of alternate readings");
         break;
     case 2: fprintf (stderr,"%s\n","Turn alternate readings into *");
         break;
  }
  switch (para) {
     case 0: fprintf (stderr,"%s\n","Leave paragraph end sign as is");
         break;
     case 1: fprintf (stderr,"%s\n","Strip paragraph end sign");
         break;
     case 2: fprintf (stderr,"%s\n","Strip paragraph end sign adding blank line");
         break;
  }
  switch (gaps) {
     case 0: fprintf (stderr,"%s\n","Leave plant gap sign as is");
         break;
     case 1: fprintf (stderr,"%s\n","Change plant gap sign to two spaces (as selected)");
         break;
  }
  if (auth == ' ') {
    fprintf (stderr,"%s\n","Ignore transcriber ID");
  } else {
    fprintf (stderr,"%s %c\n","Use only data from transcriber",auth);
  }
  return 0;
}

/*-----------------------------------------------------------*/

int getlin(buf)
/* Get line from stdin to buffer */
/* Return 0 if all OK, -1 if EOF, 1 if error */
char *buf;
{
  int cr=0, eod=0, index=0;
  int iget, itrack;
  char cget, cprev;

  int inauth=0;   /* Track transcriber */
  int inlocu=0;   /* Track locus ID */
  int skipit;     /* processing white space */

  /* Set these global parms */
  hash=0; hastext=0; hasfoli=0;
  (void) trackinit();

  while (cr == 0 && eod == 0) {
    iget = getchar();

    /* Check for end of file. Only allowed at first read */
    eod = (iget == EOF);
    if (eod) {
      if (index != 00) {
        buf[index]=0;  /* Add a null character for safety */
        return 1;
      }
      return -1; /* Here correct EOF */
    }
    cget = (char) iget;

    /* Check for hash */

    if (hash == 0) {
      if (cget == '#') {
        if (index != 00) {
          cwarn=cget;
          cget='='; /* Hash later in line translated to = */
        } else {
          hash=1;
        }
      }
    }
  
    /* Check for white space, leaving it untouched inside comments.
       At this point also process ligature brackets if desired */

    skipit=0;
    if (hash == 0 && incomm == 0) {
      if (white == 1) {
        if (cget == ' ' || cget == '\t') skipit=1;
      }

      if (liga == 1) {
        if (cget == '(' || cget == ')' ) skipit=1;
      }

      if (liga == 2) {
        if (cget == '(' ) cget= '[';
        if (cget == ')' ) cget= ']';
      }
    }

    /* Forget about this char altogether if skipit set */

    if (skipit == 0) {

      /* Some checks only if this is not a hash line */

      if (hash == 0) {

        /* Check for open caret in 1st position */

        if (cget == '<' && incomm == 0) {
          if (index != 00) {
            fprintf(stderr,"%s\n","Foliation start not in first position");
            return 1;
          } else {
            hasfoli=1;
          }
        }
  
        /* Track the text */
        itrack = track(cget,index);
        if (itrack != 0) {
          buf[index]=0;  /* Add a null character for safety */
          (void) trackerr(buf,cget);
          return 1;
        }

        /* Process locus/author inside foliation brackets */
        if (infoli != 0) {

          /* Evolve pointers */
          if (inlocu == -1) inlocu=1;
          if (inlocu == -2) inlocu=0;
          if (inauth == -1) inauth=1;
          if (inauth == -2) inauth=0;

          /* Check for < > . and ;*/
          if (cget == '<') {
            folname=" "; locname=" "; lineauth=' ';
          } else if (cget == '>') {
            /* Terminate (optionally) locus ID */
            if (inlocu >= 0) {
              locname[inlocu]=0; inlocu = -2;
            }
          } else if (cget == '.') {
            /* Start of locus */
            inlocu = -1;
          } else if (cget == ';') {
            if (inlocu <= 0) {
              fprintf(stderr,"%s\n","Start of author ID without locus");
              return 1;
            }
            locname[inlocu]=0; inlocu = -2; /* Terminate locus ID */
            inauth= -1; /* Start for author ID */
          } else {

            /* Any other char. Check where it belongs */
            if (inauth > 1) {
              fprintf(stderr,"%s\n","Illegal author ID length");
              return 1;
            } else if (inauth == 1) {
              lineauth = cget;        /* Get author ID */
            } else if (inlocu >= 0) {
              locname[inlocu-1] = cget;  /* Add to locus */
            } else {
              folname[infoli-2] = cget;   /* Add to folio */
              folname[infoli-1] = 0;      /* And prepare terminator */
            }
          }
          if (inlocu > 0) inlocu += 1;   /* Advance pointers */
          if (inauth > 0) inauth += 1;
        }
      }

      /* Add the character to the buffer */
      buf[index++] = cget;
      if (cget != ' ' && cget != '\t' && cget != '\n') hastext=1;

      /* Check for CR */
      cr = (cget == '\n');
      if (cr) {
        buf[index] = 0;
      }
    }
  }

  /* Now the entire line was read, do some more checks */

  if (hash == 0) {
    if (incomm != 0) {
      fprintf(stderr,"%s","Unclosed comment\n");
      return 1;
    }
    if (infoli != 0) {
      fprintf(stderr,"%s","Unclosed foliation comment\n");
      return 1;
    }
    if (ind_alto > ind_altc) {
      fprintf(stderr,"%s","Unclosed alternate reading\n");
      return 1;
    }
    if (ind_ligo > ind_ligc) {
      fprintf(stderr,"%s","Unclosed ligature\n");
      return 1;
    }
  }
  return 0;
}

/*-----------------------------------------------------------*/

int procread(buf)
/* Process uncertain readings and/or ligature capitalisation rule.
   Read and write to same buffer */

/* Return 0 if all OK, 1 if error */
char *buf;
{
  int index=0, itrack, ii, nrlost;
  char cb;

  /* Check if anything needs to be done at all */
  if (brack == 0 && liga != 3) return 0;

  /* If it is empty, nothing either */
  if (hastext == 0) return 0;

  /* If it is a hash line, nothing either */
  if (hash == 1) return 0;

  /* Initialise a few things */
  (void) trackinit();

  while (cb = buf[index]) {

    /* Track the text */
    itrack = track(cb,index);
    if (itrack != 0) {
      (void) trackerr(buf,cb);
      return 1;
    }

    /* If a complete set of parentheses found: */
    if (ind_ligc >= 0) {
      if (liga == 3) {

        /* Capitalise all but last char: */
        index = ind_ligc - 2;
        for (ii=ind_ligo;ii<index;ii++) {
          buf[ii]=toupper(buf[ii+1]);
        }
        buf[index]=buf[index+1];

      }

      /* Shift left remainder of line */
      nrlost = 2;
      for (ii=index+1;ii<512-nrlost;ii++) {
        buf[ii]=buf[ii+nrlost];
      }

      /* reset pointers */
      ind_ligo= -1; ind_ligc= -1;
    }

    /* If a complete set of brackets found: */
    if (ind_altc >= 0) {
      if (brack == 1) {

        /* Process into best guess: */
        index = ind_altb - 2;
        for (ii=ind_alto;ii<=index;ii++) {
          buf[ii]=buf[ii+1];
        }

      } else if (brack == 2) {

        /* Process into asterisk */
        index = ind_alto;
        buf[index] = '*'; 
      }

      /* Shift left remainder of line */
      nrlost = ind_altc - index;
      for (ii=index+1;ii<512-nrlost;ii++) {
        buf[ii]=buf[ii+nrlost];
      }

      /* reset pointers */
      ind_alto= -1; ind_altb= -1; ind_altc= -1;
    }
    index += 1;
  }
  return 0;
}

/*-----------------------------------------------------------*/

int procspc(buf1,buf2)
/* Process spaces in buffer. This treats occurrences of
   comma, dot, hyphen, equal, percent and exclam in input text */
/* Return 0 if all OK, 1 if error */
char *buf1, *buf2;
{
  int indin=0, indout=0, eol=0;
  int itrack=0, addchar;
  char cb;

  /* Some standard 'track' initialisations per line */
  (void) trackinit();

  /* If no text, do little */
  if (hastext == 0) {
    buf2[0]=0;
    return 0;
  }

  while (eol == 0) {
    cb=buf1[indin];
    eol = (cb == (char) 0);

    if (hash == 0) {

      /* Track the text */
      itrack = track(cb,indin);
      if (itrack != 0) {
        (void) trackerr(buf2,cb);
        return 1;
      }
    }
/*  if (comm == -1) fprintf(stderr,"%s","Start of comment\n");
    if (comm == -2) fprintf(stderr,"%s","End of comment\n");
    if (folio == -1) fprintf(stderr,"%s","Start of foliation\n");
    if (folio == -2) fprintf(stderr,"%s","End of foliation\n");  */

    indin += 1;
    addchar = 1;

    /* Check outside comments only */
    if (hash == 0 && incomm == 0 && infoli == 0) {

      /* Here for real and half spaces */
      if (cb == ',') {
        if (uspace == 1) cb = ' ';        /* Just use space instead */
        else if (uspace == 2) addchar=0;  /* Convert to nothing */
      }
      if (cb == '.') {
        if (hspace == 1) cb = ' ';        /* Just use space instead */
        else if (hspace == 2) addchar=0;  /* Convert to nothing */
      }

      /* Here for hyphen and equal */
      if (cb == '-') {
        if (gaps != 0) {
          cb = ' '; addchar = 2;        /* Convert to 2 chars */
        }
      }
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

int putlin(buf)
/* Write buffer to stdout, optionally skipping comments,
   foliation info and a few other things */

/* Return 0 if all OK, 1 if error */

char *buf;
{
  int index=0, eol=0, itrack, output, printed=0;
  char cb;

  /* For a block comment, exit possibly now */
  if (kcomm == 1 && hash == 1) return 0;

  /* If the line is empty, only print if whitespace kept */
  if (white == 1 && hastext == 0) return 0;
     
  /* For 'wrong author': same thing */
  if (auth != ' ') {
    if (lineauth != ' ' && lineauth != auth) return 0;
  }

  /* The usual initialisation for 'track' */
  (void) trackinit();

  while (cb = buf[index]) {
    output=1;
    itrack=track(cb,index);

    /* Check inline comments */
    if (kcomm == 1) {
      if (incomm != 0) output = 0;
    }

    /* Check foliation comments */
    if (kfoli == 1) {
      if (infoli != 0) output = 0;
    }

    /* Now output if required */
    if (output !=0) {

      /* Do not print entirely empty lines */
      if (printed != 0 || cb != '\n') {
        (void) putchar(cb); printed += 1;
      }
    }

    index += 1;
  }
  return 0;
}

/*-----------------------------------------------------------*/

int track(cb,index)
/* Keep track of comments, foliation, ligatures, etc */
/* Return 0 if OK, 1 if error */
/* Will only be called for line without # comment */

char cb;
int index;
{
  /* Evolve previous reading of bracket */
  if (incomm == -1) incomm=1;
  if (incomm == -2) incomm=0;
  if (infoli == -1) infoli=1;
  if (infoli == -2) infoli=0;

  /* Check for in-line comments. Error if inside < >,
     warning if inside [ ] or ( ) */

  if (cb == '{') {
    if (infoli !=0 ) return 1;
    if (ind_ligo > ind_ligc || ind_alto > ind_altc) cwarn='{';
    if (incomm != 0) return 1; /* Already open */
    incomm=-1; return 0;
  }
  if (cb == '}') {
    if (incomm < 1) return 1; /* Not open */
    incomm=-2; return 0;
  }

  /* Continue checks only outside { } comments */
  if (incomm !=0) {
    if (incomm > 0) incomm += 1;
    return 0;
  }

  /* Check for foliation. It is already guaranteed that nothing
     else is open */
  if (cb == '<') {
    infoli=-1; return 0;
  }
  if (cb == '>') {
    if (infoli < 3) return 1; /* Should really be longer */
    infoli=-2; return 0;
  }

  /* Check for alternate reading brackets, not allowed inside < > */
  if (cb == '[') {
    if (infoli != 0) return 1;
    if (ind_alto > ind_altc) {
      return 1; /* Unclosed bracket already open */
    }
    ind_alto = index; ind_altb = -1; ind_altc = -1;
  } else if (cb == '|') {
    if (infoli != 0) {
      return 1;
    }
    if (ind_alto < 0) return 1; /* No bracket open */
    if (ind_altb < 0) ind_altb = index; /* Keep only left-most */
  } else if (cb == ']') {
    if (infoli != 0) {
      return 1;
    }
    if (ind_alto < 0) {
      return 1; /* No bracket open */
    }
    if (ind_altb < 0) { 
      cwarn='|';         /* Vertical bar missing, issue warning */
      ind_altb=ind_alto+2; /* So that it can continue reasonably */
    }
    ind_altc = index;
  }

  /* Check for ligature brackets, not allowed inside < > */

  if (cb == '(') {
    if (infoli != 0) {
      return 1;
    }
    if (ind_ligo > ind_ligc) {
      return 1; /* Unclosed paren already open */
    }
    ind_ligo = index; ind_ligc = -1;
    return 0;
  }
  if (cb == ')') {
    if (infoli != 0) {
      return 1;
    }
    if (ind_ligo < 0) {
      return 1; /* No bracket open */
    }
    ind_ligc = index;
    return 0;
  }

  /* Inside foliation keep track of index */
  if (infoli > 0) infoli += 1;

  return 0;
}

/*-----------------------------------------------------------*/

int trackinit()

/* Initialise tracker */

{
    incomm=0   /* 1 if inside { } */;
    infoli=0   /* 1 if inside < > */;
    ind_ligo= -1 /* Position of ligatures ( char */;
    ind_ligc= -1 /* Position of ligatures ) char */;
    ind_alto= -1 /* Position of alt. readings [ char */;
    ind_altb= -1 /* Position of alt. readings | char */;
    ind_altc= -1 /* Position of alt. readings ] char */;
}

/*-----------------------------------------------------------*/

int trackerr(buf,cget)

/* Print track error */

char *buf, cget;
{
  fprintf(stderr,"%s %c\n","Track finds illegal character:",cget);
  fprintf(stderr,"%s %s\n","Line parsed so far:",buf);
}
